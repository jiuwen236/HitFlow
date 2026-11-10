# rules/ngs_fastp.smk
# Snakemake rules for fastp quality control and adapter trimming
# This rule runs before alignment to clean and filter raw FASTQ files

ENV_FASTP = "../envs/ngs/fastp.yaml"

# Container configuration
CONTAINER_FASTP = config.get("containers", {}).get("fastp", "")

# Get directory configurations from main Snakefile
NGS_PREPROCESSING_DIR = config["output"]["ngs"]["preprocessing"]
NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
LOGS_DIR = config["output"]["logs"]

# S03 Fastp Quality Control and Adapter Trimming
rule fastp_qc:
    input:
        r1 = lambda wildcards: get_fastq_path(wildcards, "1"),
        r2 = lambda wildcards: get_fastq_path(wildcards, "2")
    output:
        r1_clean = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastp/{{sample_id}}.clean.R1.fq.gz",
        r2_clean = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastp/{{sample_id}}.clean.R2.fq.gz",
        html = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastp/{{sample_id}}.fastp.html",
        json = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastp/{{sample_id}}.fastp.json"
    params:
        opts = config["fastp"]["options"],
        # Convert list to space-separated string for input files
        r1_input = lambda wildcards, input: " ".join(input.r1) if isinstance(input.r1, list) else input.r1,
        r2_input = lambda wildcards, input: " ".join(input.r2) if isinstance(input.r2, list) else input.r2
    threads: config["fastp"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["fastp"]["mem"].rstrip("g")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["fastp_module"]
    conda:
        ENV_FASTP
    container:
        CONTAINER_FASTP
    log:
        f"{LOGS_DIR}/fastp/{{sample_id}}.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/fastp/{{sample_id}}.txt"
    shell:
        """
        set -euo pipefail

        tmp_merge=""
        cleanup() {{
            if [ -n "$tmp_merge" ] && [ -d "$tmp_merge" ]; then
                rm -rf "$tmp_merge"
            fi
        }}
        trap cleanup EXIT

        if [ $(echo {params.r1_input} | wc -w) -gt 1 ]; then
            echo "Merging multiple R1 files: {params.r1_input}" >> {log}
            echo "Merging multiple R2 files: {params.r2_input}" >> {log}

            tmp_merge=$(mktemp -d)
            r1_concat="$tmp_merge/{wildcards.sample_id}.merged.R1.fastq.gz"
            r2_concat="$tmp_merge/{wildcards.sample_id}.merged.R2.fastq.gz"

            cat {params.r1_input} > "$r1_concat"
            cat {params.r2_input} > "$r2_concat"

            in_r1="$r1_concat"
            in_r2="$r2_concat"
        else
            in_r1="{params.r1_input}"
            in_r2="{params.r2_input}"
        fi

        fastp \
          -i "$in_r1" \
          -I "$in_r2" \
          -o {output.r1_clean} \
          -O {output.r2_clean} \
          -w {threads} \
          -h {output.html} \
          -j {output.json} \
          {params.opts} \
          2>> {log}
        """
