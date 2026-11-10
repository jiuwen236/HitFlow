# rules/workflow3 qc.smk
# Snakemake rules for QC steps: qualimap_bamqc, qualimap_filter, verify_bam_id, verify_bam_id_filter
ENV_QUALIMAP = "../envs/ngs/qualimap.yaml"
ENV_PERL = "../envs/basic.yaml"
ENV_VERIFYBAMID = "../envs/ngs/verifybamid.yaml"

# Container configurations
CONTAINER_QUALIMAP = config.get("containers", {}).get("qualimap", "")
CONTAINER_VERIFYBAMID = config.get("containers", {}).get("verifybamid2", "")
CONTAINER_PERL = config.get("containers", {}).get("perl_json", "")

# Get directory configurations from main Snakefile
NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
QC_FASTQC_DIR = config["output"]["ngs"]["preprocessing"]
LOGS_DIR = config["output"]["logs"]

# VerifyBamID2 enabled flag
VERIFYBAMID2_ENABLED = config.get("verifybamid2", {}).get("enabled", True)

# S14 Qualimap BamQC
rule qualimap_bamqc:
    input:
        bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bam"
    output:
        genome_results = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/qualimap/genome_results.txt",
        qualimap_report = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/qualimap/qualimapReport.html"
    params:
        outdir = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/qualimap",
        chroms = ",".join(config["chromosomes"]),
        ngs_scripts = config["ngs_scripts"]
    threads: 8
    resources:
        mem_mb = 45 * 1024,  # 45g
        time = "200:00:00"
    envmodules:
        config["tool_modules"]["qualimap_module"],
        config["tool_modules"]["java_module"]
    conda:
        ENV_QUALIMAP
    container:
        CONTAINER_QUALIMAP
    benchmark:
        f"{LOGS_DIR}/benchmarks/qualimap/{{sample_id}}.bamqc.txt"
    log:
        f"{LOGS_DIR}/qualimap/{{sample_id}}.bamqc.log"
    shell:
        """
        mkdir -p {params.outdir}
        
        qualimap bamqc --java-mem-size=45g -nt {threads} \
            -bam {input.bam} -outdir {params.outdir} > {log} 2>&1
        """
        

# S15 Qualimap Filter
def get_qualimap_filter_inputs(wildcards):
    """Get inputs for qualimap_filter. Skip FastQC for BAM-only samples."""
    inputs = {
        "qualimap_genome_results": f"{NGS_ALIGNMENT_DIR}/{wildcards.sample_id}/qc/qualimap/genome_results.txt",
        "qualimap_report": f"{NGS_ALIGNMENT_DIR}/{wildcards.sample_id}/qc/qualimap/qualimapReport.html"
    }
    # Only add FastQC dependency for samples with FASTQ files
    if not sample_has_bam(wildcards.sample_id):
        inputs["fastqc_report"] = f"{NGS_PREPROCESSING_DIR}/{wildcards.sample_id}/fastqc/{wildcards.sample_id}.fastqc.report"
    return inputs

rule qualimap_filter:
    input:
        unpack(get_qualimap_filter_inputs)
    output:
        report = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/{{sample_id}}.qualimap.report"
    params:
        flag_success = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/qualimap.success",
        flag_failed = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/qualimap.failed",
        fastqc_dir = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc",
        qualimap_dir = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/qualimap",
        ngs_scripts = config["ngs_scripts"],
        configure_path = config["configure_path"],
        qc_results_dir = config["output"]["ngs"]["qc_reports"],
        failed_samples_dir = config["output"]["ngs"]["qc_reports"] + "/failed_samples",
        fail_log = config["output"]["ngs"]["qc_reports"] + "/failed_samples.txt"
    resources:
        mem_mb = 1024,
        time = "144:00:00"
    envmodules:
        config["tool_modules"]["perl_module"]
    conda:
        ENV_PERL
    container:
        CONTAINER_PERL
    benchmark:
        f"{LOGS_DIR}/benchmarks/qualimap/{{sample_id}}.filter.txt"
    log:
        f"{LOGS_DIR}/qualimap/{{sample_id}}.filter.log"
    shell:
        """
        mkdir -p {params.qualimap_dir}
        
        # Run the qualimap filter script (allow non-zero exit for QC failures)
        set +e
        perl {params.ngs_scripts}/qualimap_verifybamid_filter.pl \
            -qualimap1 {input.qualimap_genome_results} \
            -qualimap2 {input.qualimap_report} \
            -fastqcdir {params.fastqc_dir} \
            -tool qualimap \
            -conf {params.configure_path} \
            -out {output.report} > {log} 2>&1
            
        EXIT_CODE=$?
        set -e
        sleep 1
        
        if [ $EXIT_CODE -eq 0 ]; then
            echo "success" >> {log}
            mkdir -p {params.qc_results_dir}
            cp -af {params.qualimap_dir} {params.qc_results_dir}/
            touch {params.flag_success}
        else
            echo "{wildcards.sample_id} qualimap '$(date)'" >> {params.fail_log}
            mkdir -p {params.failed_samples_dir}
            cp -af {params.qualimap_dir} {params.failed_samples_dir}/
            touch {params.flag_failed}
            echo "WARNING: Qualimap QC failed for {wildcards.sample_id}; results kept for review" >> {log}
        fi
        """

# S16 Verify BAM ID (only if enabled)
if VERIFYBAMID2_ENABLED:
    rule verify_bam_id:
        input:
            bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.recal.bam",
            ref = config["resolved_paths"]["reference_genome"]
        output:
            selfSM = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/verifybamid/{{sample_id}}.selfSM",
            json = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/verifybamid/{{sample_id}}.json"
        params:
            outdir = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/verifybamid",
            outprefix = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/verifybamid/{{sample_id}}",
            ngs_scripts = config["ngs_scripts"],
            svd_prefix = config["resolved_paths"]["svd_prefix"]
        threads: 4
        resources:
            mem_mb = 20 * 1024,  # 20g
            time = "48:00:00"
        envmodules:
            config["tool_modules"]["verifybamid_module"]
        conda:
            ENV_VERIFYBAMID
        container:
            CONTAINER_VERIFYBAMID
        benchmark:
            f"{LOGS_DIR}/benchmarks/verify_bam_id/{{sample_id}}.txt"
        log:
            f"{LOGS_DIR}/verify_bam_id/{{sample_id}}.log"
        shell:
            """
            mkdir -p {params.outdir}
            
            # Detect verifybamid2 or verifyBamID
            if command -v verifybamid2 >/dev/null 2>&1; then
                VERIFYBAMID=verifybamid2
            else
                VERIFYBAMID=verifyBamID
            fi
            
            $VERIFYBAMID --Reference {input.ref} \
                --BamFile {input.bam} \
                --Output {params.outprefix} \
                --NumThread {threads} \
                --SVDPrefix {params.svd_prefix} > {log} 2>&1
            
            if command -v python3 >/dev/null 2>&1; then
                PYTHON_BIN=python3
            elif command -v python >/dev/null 2>&1; then
                PYTHON_BIN=python
            else
                PYTHON_BIN=""
            fi

            if [ -n "$PYTHON_BIN" ]; then
                echo "INFO: verify_bam_id will use $PYTHON_BIN for JSON conversion" >> {log}
                $PYTHON_BIN {params.ngs_scripts}/verifybamid_json.py \
                    {output.selfSM} {output.json} >> {log} 2>&1
            else
                echo "WARNING: python not found, using sed/cut fallback for JSON conversion" >> {log}
                # Use sed/cut instead of awk for better compatibility with minimal containers
                SAMPLE_ID=$(grep -v "^#" {output.selfSM} | head -1 | cut -f1)
                FREEMIX=$(grep -v "^#" {output.selfSM} | head -1 | cut -f7)
                printf '%s' '{{"'"$SAMPLE_ID"'":{{"contaminate":"'"$FREEMIX"'"}}}}' > {output.json}
                echo "INFO: sed/cut fallback finished writing {output.json}" >> {log}
            fi
            """

    # S17 Verify BAM ID Filter
    rule verify_bam_id_filter:
        input:
            selfSM = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/verifybamid/{{sample_id}}.selfSM"
        output:
            report = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/{{sample_id}}.verifybamid.report",
            flag_success = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/verifybamid.success"
        params:
            flag_failed = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/verifybamid.failed",
            outdir = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/qc/verifybamid",
            ngs_scripts = config["ngs_scripts"],
            configure_path = config["configure_path"],
            qc_results_dir = config["output"]["ngs"]["qc_reports"],
            failed_samples_dir = config["output"]["ngs"]["qc_reports"] + "/failed_samples",
            fail_log = config["output"]["ngs"]["qc_reports"] + "/failed_samples.txt"
        resources:
            mem_mb = 1024,
            time = "144:00:00"
        envmodules:
            config["tool_modules"]["perl_module"]
        conda:
            ENV_PERL
        container:
            CONTAINER_PERL
        benchmark:
            f"{LOGS_DIR}/benchmarks/verify_bam_id/{{sample_id}}.filter.txt"
        log:
            f"{LOGS_DIR}/verify_bam_id/{{sample_id}}.filter.log"
        shell:
            """
            mkdir -p {params.outdir}
            
            perl {params.ngs_scripts}/qualimap_verifybamid_filter.pl \
                -verifybamid {input.selfSM} \
                -tool verifybamid \
                -conf {params.configure_path} \
                -out {output.report} > {log} 2>&1
                
            EXIT_CODE=$?
            sleep 1
            
            if [ $EXIT_CODE -eq 0 ]; then
                echo "success" >> {log}
                mkdir -p {params.qc_results_dir}
                cp -a {params.outdir} {params.qc_results_dir}/
                touch {output.flag_success}
            else
                echo "{wildcards.sample_id} verifybamid '$(date)'" >> {params.fail_log}
                mkdir -p {params.failed_samples_dir}
                cp -a {params.outdir} {params.failed_samples_dir}/
                touch {params.flag_failed}
                exit 1
            fi
            """ 