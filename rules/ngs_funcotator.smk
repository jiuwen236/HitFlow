# rules/ngs_funcotator.smk
# GATK Funcotator annotation for single-sample variants

ENV_FUNCOTATOR = "../envs/ngs/funcotator.yaml"

# Container configurations
CONTAINER_GATK4 = config.get("containers", {}).get("gatk4", "")
CONTAINER_BCFTOOLS = config.get("containers", {}).get("bcftools", "")

import os

# Get directory configurations
NGS_VARIANTS_SINGLE_DIR = config["output"]["ngs"]["variants"]["single_sample"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# Get reference paths
REF_GENOME = config["resolved_paths"]["reference_genome"]
REF_FAI = REF_GENOME + ".fai"
REF_DICT = os.path.splitext(REF_GENOME)[0] + ".dict"

# Get Funcotator configuration
FUNCOTATOR_CFG = config.get("ngs_single_sample_variants", {}).get("annotation", {}).get("funcotator", {})


# S01: Funcotator annotation
# S01a: Funcotator annotation (output uncompressed VCF)
rule funcotate_variants_gatk:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.vcf.gz.tbi",
        ref = REF_GENOME,
        fai = REF_FAI,
        dict = REF_DICT
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.funcotator.vcf"
    params:
        java_opts = FUNCOTATOR_CFG.get("java_opts", "-Xmx10G -XX:ParallelGCThreads=1"),
        data_sources = FUNCOTATOR_CFG.get("data_sources", "ref_data/funcotator_dataSources"),
        ref_version = FUNCOTATOR_CFG.get("ref_version", "hg38"),
        transcript_selection = FUNCOTATOR_CFG.get("transcript_selection_mode", "CANONICAL"),
        annotation_defaults = FUNCOTATOR_CFG.get("annotation_defaults", []),
        annotation_overrides = FUNCOTATOR_CFG.get("annotation_overrides", []),
        exclude_fields = FUNCOTATOR_CFG.get("exclude_fields", []),
        output_format = "VCF"  # Force VCF output
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(FUNCOTATOR_CFG.get("mem", "12G").rstrip("GgMmKk")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["java_module"],
        config["tool_modules"]["gatk_module"]
    conda:
        ENV_FUNCOTATOR
    container:
        CONTAINER_GATK4
    benchmark:
        f"{LOGS_DIR}/benchmarks/funcotator/{{sample_id}}.funcotate_gatk.tsv"
    log:
        f"{LOGS_DIR}/funcotator/{{sample_id}}.funcotate_gatk.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/funcotator_{wildcards.sample_id}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.vcf})
        
        # Build annotation defaults arguments
        ANNO_DEFAULTS=""
        if [ -n "{params.annotation_defaults}" ]; then
            for anno in {params.annotation_defaults}; do
                ANNO_DEFAULTS="$ANNO_DEFAULTS --annotation-default $anno"
            done
        fi
        
        # Build annotation overrides arguments
        ANNO_OVERRIDES=""
        if [ -n "{params.annotation_overrides}" ]; then
            for anno in {params.annotation_overrides}; do
                ANNO_OVERRIDES="$ANNO_OVERRIDES --annotation-override $anno"
            done
        fi
        
        # Build exclude fields arguments
        EXCLUDE_FIELDS=""
        if [ -n "{params.exclude_fields}" ]; then
            for field in {params.exclude_fields}; do
                EXCLUDE_FIELDS="$EXCLUDE_FIELDS --exclude-field $field"
            done
        fi
        
        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" Funcotator \\
            --variant {input.vcf} \\
            --reference {input.ref} \\
            --ref-version {params.ref_version} \\
            --data-sources-path {params.data_sources} \\
            --output {output.vcf} \\
            --output-file-format {params.output_format} \\
            --transcript-selection-mode {params.transcript_selection} \\
            $ANNO_DEFAULTS \\
            $ANNO_OVERRIDES \\
            $EXCLUDE_FIELDS \\
            > {log} 2>&1
        
        rm -rf ${{JOB_TMP_DIR}}
        """

# S01b: Compress and index Funcotator output with bcftools
rule funcotate_variants_compress:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.funcotator.vcf"
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.funcotator.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.funcotator.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_FUNCOTATOR
    container:
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/funcotator/{{sample_id}}.funcotate_compress.tsv"
    log:
        f"{LOGS_DIR}/funcotator/{{sample_id}}.funcotate_compress.log"
    shell:
        """
        # Compress the VCF with bcftools
        bcftools view -O z -o {output.vcf} {input.vcf} > {log} 2>&1
        
        # Create index
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """


# S02: Extract Funcotator annotations to TSV table
rule extract_funcotator_table:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.funcotator.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.funcotator.vcf.gz.tbi"
    output:
        tsv = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.funcotator.tsv"
    params:
        # Key fields to extract from Funcotator annotations
        fields = "CHROM,POS,ID,REF,ALT,QUAL,FILTER,FUNCOTATION"
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_FUNCOTATOR
    container:
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/funcotator/{{sample_id}}.extract_table.tsv"
    log:
        f"{LOGS_DIR}/funcotator/{{sample_id}}.extract_table.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        
        # Extract annotations and parse Funcotator field
        bcftools query -H -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/FUNCOTATION\\n' \\
            {input.vcf} > {output.tsv} 2> {log}
        
        # If bcftools query fails or produces empty output, create minimal output
        if [ ! -s {output.tsv} ]; then
            echo "Warning: No variants with Funcotator annotations found" >> {log}
            echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tFUNCOTATION" > {output.tsv}
        fi
        """
