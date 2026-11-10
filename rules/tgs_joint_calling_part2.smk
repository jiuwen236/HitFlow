# rules/tgs_joint_calling_part2.smk
# TGS Joint Calling - Part 2: Normalization and Filtering
# Note: VQSR steps removed as DeepVariant/Clair3 provide quality scores

ENV_BCFTOOLS = "../envs/basic.yaml"
ENV_VCFTOOLS = "../envs/ngs/vcftools.yaml"

import os

TGS_VARIANTS_SINGLE_DIR = config["output"]["tgs"]["variants"]["single_sample"]
TGS_VARIANTS_COHORT_DIR = config["output"]["tgs"]["variants"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]


# TGS-04: VCF normalization with bcftools norm
rule tgs_normalize_vcf:
    input:
        vcf = f"{TGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz",
        idx = f"{TGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz.tbi",
        ref = config["resolved_paths"]["reference_genome"],
        fai = config["resolved_paths"]["reference_genome"] + ".fai"
    output:
        vcf = f"{TGS_VARIANTS_COHORT_DIR}/normalized/cohort.normalized.vcf.gz",
        idx = f"{TGS_VARIANTS_COHORT_DIR}/normalized/cohort.normalized.vcf.gz.tbi"
    params:
        options = config["tgs_joint_calling"]["vcf_processing"]["normalize"]["options"]
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BCFTOOLS
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/tgs_joint_calling/normalize_vcf.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_joint_calling/normalize_vcf.benchmark.txt"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        bcftools norm \\
            --fasta-ref {input.ref} \\
            {params.options} \\
            --output-type z \\
            --output {output.vcf} \\
            {input.vcf} \\
            > {log} 2>&1

        # Index the output VCF
        bcftools index --tbi --force {output.vcf}
        """


# TGS-06: Quality filtering with vcftools
rule tgs_filter_variants:
    input:
        vcf = f"{TGS_VARIANTS_COHORT_DIR}/normalized/cohort.normalized.vcf.gz",
        idx = f"{TGS_VARIANTS_COHORT_DIR}/normalized/cohort.normalized.vcf.gz.tbi"
    output:
        vcf = f"{TGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf"
    params:
        max_missing = config["tgs_joint_calling"]["vcf_processing"]["filtering"]["max_missing"],
        hwe_threshold = config["tgs_joint_calling"]["vcf_processing"]["filtering"]["hwe_threshold"],
        out_prefix = f"{TGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered"
    threads: 1
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["vcftools_module"]
    conda:
        ENV_VCFTOOLS
    container:
        config["containers"]["vcftools"]
    log:
        f"{LOGS_DIR}/tgs_joint_calling/filter_variants.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_joint_calling/filter_variants.benchmark.txt"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        # Run vcftools filtering
        vcftools \\
            --recode \\
            --recode-INFO-all \\
            --gzvcf {input.vcf} \\
            --max-missing {params.max_missing} \\
            --hwe {params.hwe_threshold} \\
            --out {params.out_prefix} \\
            > {log} 2>&1
        """


# TGS-06b: Compress and index filtered VCF with bcftools
rule tgs_compress_filtered_vcf:
    input:
        vcf = f"{TGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf"
    output:
        vcf = f"{TGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz",
        idx = f"{TGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BCFTOOLS
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/tgs_joint_calling/compress_filtered_vcf.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_joint_calling/compress_filtered_vcf.benchmark.txt"
    shell:
        """
        # Compress the output VCF with bcftools and create index
        bcftools view -O z -o {output.vcf} {input.vcf} > {log} 2>&1
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """


def get_tgs_joint_calling_outputs():
    """Get TGS joint calling outputs based on configuration and sample count."""
    if not config["workflows"].get("enable_tgs_joint_calling", False):
        return []
    if not config["tgs_joint_calling"].get("enabled", False):
        return []
    
    try:
        samples = get_tgs_joint_samples()
    except (NameError, KeyError):
        return []
    
    if len(samples) < 2:  # Require at least 2 samples
        return []
    
    outputs = [f"{TGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz"]
    # Per-chromosome genotyped VCFs only exist in GATK backend
    if config.get("tgs_joint_calling", {}).get("tool", "glnexus").strip().lower() == "gatk":
        outputs.extend([f"{TGS_VARIANTS_COHORT_DIR}/genotyped/{chromosome}.genotyped.vcf.gz"
                       for chromosome in config["chromosomes"]])

    phasing_cfg = config.get("tgs_joint_calling", {}).get("phasing", {})
    if config["workflows"].get("enable_tgs_phasing", False) and phasing_cfg.get("enabled", False):
        skip_chroms = config.get("phasing", {}).get("skip_chromosomes", [])
        outputs.extend([
            f"{TGS_VARIANTS_COHORT_DIR}/phased/{chromosome}.phased.vcf.gz"
            for chromosome in config["chromosomes"]
            if chromosome not in skip_chroms
        ])
    
    return outputs

