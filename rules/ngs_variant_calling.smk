# rules/ngs_variant_calling.smk
# Snakemake rules for GATK4 HaplotypeCaller variant calling

ENV_GATK4 = "../envs/ngs/gatk4.yaml"

# Get directory configurations from main Snakefile
NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
NGS_PREPROCESSING_DIR = config["output"]["ngs"]["preprocessing"]
NGS_VARIANTS_SINGLE_DIR = config["output"]["ngs"]["variants"]["single_sample"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]
DBSNP_RESOURCE = config["resolved_paths"].get("known_sites_vcf", "")

# VerifyBamID2 enabled flag
VERIFYBAMID2_ENABLED = config.get("verifybamid2", {}).get("enabled", True)

# Helper function to get optional verifybamid dependency
def get_verifybamid_input(wildcards):
    if VERIFYBAMID2_ENABLED:
        return f"{NGS_ALIGNMENT_DIR}/{wildcards.sample_id}/qc/verifybamid.success"
    return []

# S19 GATK4 HaplotypeCaller: Call variants per chromosome and create GVCFs
rule gatk4_haplotype_caller:
    input:
        bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.recal.bam",
        bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.recal.bai",
        ref = config["resolved_paths"]["reference_genome"],
        fai = config["resolved_paths"]["reference_genome"] + ".fai",
        dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict",
        dbsnp = lambda wildcards: DBSNP_RESOURCE if DBSNP_RESOURCE else [],
        verifybamid_success = get_verifybamid_input
    output:
        gvcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/gvcf/{{sample_id}}.{{chromosome}}.hc.g.vcf.bgz"
    params:
        java_opts = config["gatk_paras"]["haplotypecaller"]["java_opts"],
        contamination = config["gatk_paras"]["haplotypecaller"]["contamination"],
        gqb_params = lambda wildcards: " ".join([f"-GQB {i}" for i in config["gatk_paras"]["haplotypecaller"]["gqb_values"]]),
        annotation_params = lambda wildcards: " ".join([f"-G {a}" for a in config["gatk_paras"]["haplotypecaller"]["annotations"]]),
        wes_bed = config.get("resolved_paths", {}).get("wes", {}).get("targets_bed", ""),
        wes_pad = config.get("resolved_paths", {}).get("wes", {}).get("interval_padding", 0),
        dbsnp_arg = lambda wildcards, input: (f"--dbsnp {input.dbsnp}" if input.dbsnp else "")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(config["gatk_paras"]["haplotypecaller"]["mem"].rstrip("GgMmKk")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["java_module"],
        config["tool_modules"]["gatk_module"]
    conda:
        ENV_GATK4
    container:
        config["containers"]["gatk4"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/gatk_haplotypecaller/{{sample_id}}.{{chromosome}}.txt"
    log:
        f"{LOGS_DIR}/gatk_haplotypecaller/{{sample_id}}.{{chromosome}}.haplotypecaller.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/gatk_hc_{wildcards.sample_id}_{wildcards.chromosome}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.gvcf})

        # Build interval arguments; when WES BED is provided, use INTERSECTION to
        # restrict calls to the overlap of chromosome and BED targets
        INTERVALS="--intervals {wildcards.chromosome}"
        INTERVAL_SET_RULE=""
        if [ -n "{params.wes_bed}" ]; then
            INTERVALS="$INTERVALS -L {params.wes_bed} --interval-padding {params.wes_pad}"
            INTERVAL_SET_RULE="--interval-set-rule INTERSECTION"
        fi

        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" HaplotypeCaller \
            --reference {input.ref} \
            {params.dbsnp_arg} \
            --input {input.bam} \
            --contamination {params.contamination} {params.annotation_params} {params.gqb_params} \
            --output {output.gvcf} \
            $INTERVALS \
            $INTERVAL_SET_RULE \
            --emit-ref-confidence GVCF \
            > {log} 2>&1 || true

        # When INTERSECTION is empty (e.g., chromosome not covered in WES BED),
        # GATK returns non-zero. Create an empty GVCF so downstream rules can proceed.
        # Use Python gzip since bgzip may not be available in the GATK container.
        if [ ! -s {output.gvcf} ]; then
            python3 -c "
import gzip
header = b'##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t{wildcards.sample_id}\\n'
with gzip.open('{output.gvcf}', 'wb') as f:
    f.write(header)
"
        fi
        
        EXIT_CODE=$?

        # If successful, append/update samples_vcf.csv with rawHaplotypeCaller directory
        if [ $EXIT_CODE -eq 0 ]; then
            SAMPLES_VCF_CSV="{config[samples_vcf_csv]}"
            if [ -n "$SAMPLES_VCF_CSV" ]; then
                python - << PY
import csv, os

csv_path = "$SAMPLES_VCF_CSV"
sample = "{wildcards.sample_id}"
raw_dir = os.path.dirname("{output.gvcf}")

fieldnames = ["sample", "raw_hc_dir"]
rows = []
found = False

if os.path.exists(csv_path):
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            if r["sample"] == sample:
                r["raw_hc_dir"] = raw_dir
                found = True
            rows.append(r)

if not found:
    rows.append({{"sample": sample, "raw_hc_dir": raw_dir}})

with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
PY
            fi
        fi

        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """ 