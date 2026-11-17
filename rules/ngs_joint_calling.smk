# rules/ngs_joint_calling.smk
# Snakemake rules for GATK joint calling, VQSR, and population-level variant analysis

ENV_GATK4 = "../envs/ngs/gatk4.yaml"
ENV_BCFTOOLS = "../envs/basic.yaml"
ENV_VCFTOOLS = "../envs/ngs/vcftools.yaml"

# Container configurations
CONTAINER_GATK4 = config.get("containers", {}).get("gatk4", "")
CONTAINER_BCFTOOLS = config.get("containers", {}).get("bcftools", "")
CONTAINER_VCFTOOLS = config.get("containers", {}).get("vcftools", "")
CONTAINER_SHAPEIT5 = config.get("containers", {}).get("shapeit5", "")

import os
from pathlib import Path
import csv
import glob

# Get directory configurations from main Snakefile
NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
NGS_VARIANTS_SINGLE_DIR = config["output"]["ngs"]["variants"]["single_sample"]
NGS_VARIANTS_COHORT_DIR = config["output"]["ngs"]["variants"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# File path resolution function (imported from main Snakefile logic)
def resolve_input_file_path(file_path, base_dirs=None):
    """
    Resolve input file path by checking if it exists, and if not, try with base directories.
    
    Args:
        file_path (str): The file path to resolve
        base_dirs (list): List of base directories to try if file doesn't exist
    
    Returns:
        str: The resolved file path
    """
    if base_dirs is None:
        base_dirs = config["input_data"]["base_dir"]
    
    # If it's an absolute path, return as is
    if os.path.isabs(file_path):
        return file_path
    
    # If the relative path exists, return as is
    if os.path.exists(file_path):
        return file_path
    
    # Try with each base directory
    for base_dir in base_dirs:
        full_path = os.path.join(base_dir, file_path)
        if os.path.exists(full_path):
            return full_path
    
    # If file not found, return original path (let snakemake handle the error)
    return file_path

# ---------- Joint sample discovery (FASTQ + VCF + scan_ngs_dirs) ----------
BASE_DIRS = config["input_data"]["base_dir"]
SCAN_DIRS = config.get("scan_ngs_dirs", [])
SAMPLES_VCF_CSV = resolve_input_file_path(config.get("samples_vcf_csv", "samples_vcf.csv"))

def _load_vcf_samples(csv_path):
    dirs = {}
    if os.path.isfile(csv_path):
        with open(csv_path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                dirs[row["sample"]] = row["raw_hc_dir"]
    return dirs

VCF_SAMPLE_DIRS = _load_vcf_samples(SAMPLES_VCF_CSV)


def get_joint_samples():
    """Union of samples defined via FASTQ pipeline, VCF CSV file and scan_ngs_dirs."""
    # Get samples from the main workflow's ALL_SAMPLES variable
    try:
        main_samples = set(ALL_SAMPLES.keys()) if 'ALL_SAMPLES' in globals() else set()
    except NameError:
        main_samples = set()
    
    samples = main_samples | set(VCF_SAMPLE_DIRS.keys())
    
    # Filter out empty strings to avoid scanning current directory
    scan_dirs = [d for d in SCAN_DIRS if d and d.strip()]
    
    for d in scan_dirs:
        for vcf_file in glob.glob(os.path.join(d, "*.hc.g.vcf.bgz")):
            name = os.path.basename(vcf_file).split(".")[0]
            samples.add(name)
    return sorted(samples)


def get_gvcf_path(sample, chrom):
    """Get GVCF path for a sample and chromosome, checking VCF CSV first, then new directory structure"""
    if sample in VCF_SAMPLE_DIRS:
        return os.path.join(VCF_SAMPLE_DIRS[sample], f"{sample}.{chrom}.hc.g.vcf.bgz")
    return f"{NGS_VARIANTS_SINGLE_DIR}/{sample}/gvcf/{sample}.{chrom}.hc.g.vcf.bgz"


def get_all_sample_gvcfs_for_chrom(wildcards):
    """Return list of GVCF paths for the requested chromosome."""
    return [get_gvcf_path(s, wildcards.chromosome) for s in get_joint_samples()]

# S01: GATK GenomicsDBImport - Import GVCFs into genomics database per chromosome
rule genomics_db_import:
    input:
        gvcfs = get_all_sample_gvcfs_for_chrom,
        ref = config["resolved_paths"]["reference_genome"],
        fai = config["resolved_paths"]["reference_genome"] + ".fai",
        dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict"
    output:
        db_dir = directory(f"{NGS_VARIANTS_COHORT_DIR}/genomicsdb/{{chromosome}}"),
        sample_map = f"{NGS_VARIANTS_COHORT_DIR}/sample_maps/{{chromosome}}.sample_map"
    params:
        java_opts = config["joint_calling"]["genomicsdbimport"]["java_opts"],
        batch_size = config["joint_calling"]["genomicsdbimport"]["batch_size"],
        reader_threads = config["joint_calling"]["genomicsdbimport"]["reader_threads"],
        sample_ids = lambda wildcards: " ".join(get_joint_samples()),
        wes_bed = config.get("resolved_paths", {}).get("wes", {}).get("targets_bed", ""),
        wes_pad = config.get("resolved_paths", {}).get("wes", {}).get("interval_padding", 0)
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(config["joint_calling"]["genomicsdbimport"]["mem"].rstrip("GgMmKk")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["java_module"],
        config["tool_modules"]["gatk_module"]
    conda:
        ENV_GATK4
    container:
        CONTAINER_GATK4
    benchmark:
        f"{LOGS_DIR}/benchmarks/joint_calling/genomicsdb_import.{{chromosome}}.txt"
    log:
        f"{LOGS_DIR}/joint_calling/genomicsdb_import.{{chromosome}}.log"
    shell:
        """
        # Create sample map file
        mkdir -p $(dirname {output.sample_map})
        
        # Generate sample map
        python - << PY
import os, csv
chrom = "{wildcards.chromosome}"
out_path = "{output.sample_map}"
vcf_dirs = {{}}
csv_path = "{config[samples_vcf_csv]}"
variants_single_dir = "{NGS_VARIANTS_SINGLE_DIR}"

if os.path.isfile(csv_path):
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            vcf_dirs[r["sample"]] = r["raw_hc_dir"]

samples = "{params.sample_ids}".split()
with open(out_path, "w") as out_f:
    for s in samples:
        if s in vcf_dirs:
            path = os.path.join(vcf_dirs[s], f"{{s}}.{{chrom}}.hc.g.vcf.bgz")
        else:
            path = f"{{variants_single_dir}}/{{s}}/gvcf/{{s}}.{{chrom}}.hc.g.vcf.bgz"
        out_f.write(f"{{s}}\\t{{path}}\\n")
PY
        
        JOB_TMP_DIR="{TMP_DIR}/genomicsdb_{wildcards.chromosome}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.db_dir})

        INTERVALS="-L {wildcards.chromosome}"
        if [ -n "{params.wes_bed}" ]; then
            INTERVALS="$INTERVALS -L {params.wes_bed} --interval-padding {params.wes_pad}"
        fi

        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" GenomicsDBImport \\
            -R {input.ref} \\
            --sample-name-map {output.sample_map} \\
            $INTERVALS \\
            --batch-size {params.batch_size} \\
            --reader-threads {params.reader_threads} \\
            --merge-input-intervals \\
            --consolidate \\
            --genomicsdb-shared-posixfs-optimizations true \\
            --genomicsdb-workspace-path {output.db_dir} \\
            --tmp-dir ${{JOB_TMP_DIR}} \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """

# S02: GATK GenotypeGVCFs - Perform joint genotyping on genomics database
rule genotype_gvcfs:
    input:
        db_dir = f"{NGS_VARIANTS_COHORT_DIR}/genomicsdb/{{chromosome}}",
        ref = config["resolved_paths"]["reference_genome"],
        fai = config["resolved_paths"]["reference_genome"] + ".fai",
        dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict"
    output:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/genotyped/{{chromosome}}.genotyped.vcf.gz"
    params:
        java_opts = config["joint_calling"]["genotypegvcfs"]["java_opts"],
        wes_bed = config.get("resolved_paths", {}).get("wes", {}).get("targets_bed", ""),
        wes_pad = config.get("resolved_paths", {}).get("wes", {}).get("interval_padding", 0)
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(config["joint_calling"]["genotypegvcfs"]["mem"].rstrip("GgMmKk")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["java_module"],
        config["tool_modules"]["gatk_module"]
    conda:
        ENV_GATK4
    container:
        CONTAINER_GATK4
    benchmark:
        f"{LOGS_DIR}/benchmarks/joint_calling/genotype_gvcfs.{{chromosome}}.txt"
    log:
        f"{LOGS_DIR}/joint_calling/genotype_gvcfs.{{chromosome}}.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/genotype_{wildcards.chromosome}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.vcf})

        if [ -n "${{CONDA_PREFIX:-}}" ]; then
            export PATH="${{CONDA_PREFIX}}/bin:$PATH"
            export JAVA_HOME="${{CONDA_PREFIX}}"
        fi 

        INTERVALS="--intervals {wildcards.chromosome}"
        if [ -n "{params.wes_bed}" ]; then
            INTERVALS="$INTERVALS -L {params.wes_bed} --interval-padding {params.wes_pad}"
        fi

        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" GenotypeGVCFs \\
            --reference {input.ref} \\
            $INTERVALS \\
            --variant gendb://{input.db_dir} \\
            --output {output.vcf} \\
            --annotation AS_BaseQualityRankSumTest \\
            --annotation AS_FisherStrand \\
            --annotation AS_InbreedingCoeff \\
            --annotation AS_MappingQualityRankSumTest \\
            --annotation AS_QualByDepth \\
            --annotation AS_ReadPosRankSumTest \\
            --annotation AS_RMSMappingQuality \\
            --annotation AS_StrandOddsRatio \\
            --annotation Coverage \\
            --annotation ExcessHet \\
            --annotation FisherStrand \\
            --annotation MappingQualityRankSumTest \\
            --annotation QualByDepth \\
            --annotation ReadPosRankSumTest \\
            --annotation RMSMappingQuality \\
            --annotation StrandOddsRatio \\
            --allow-old-rms-mapping-quality-annotation-data \\
            --merge-input-intervals \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """

# S03a: GATK GatherVcfs - Gather per-chromosome VCFs into single cohort VCF  
rule gather_cohort_vcfs_gatk:
    input:
        vcfs = expand(f"{NGS_VARIANTS_COHORT_DIR}/genotyped/{{chromosome}}.genotyped.vcf.gz", 
                     chromosome=config["chromosomes"])
    output:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.unindexed.vcf.gz"
    params:
        java_opts = config["joint_calling"]["gathervcfs"]["java_opts"],
        input_args = lambda wildcards, input: " ".join([f"--INPUT {vcf}" for vcf in input.vcfs])
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(config["joint_calling"]["gathervcfs"]["mem"].rstrip("GgMmKk")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["java_module"],
        config["tool_modules"]["gatk_module"]
    conda:
        ENV_GATK4
    container:
        CONTAINER_GATK4
    benchmark:
        f"{LOGS_DIR}/benchmarks/joint_calling/gather_vcfs_gatk.txt"
    log:
        f"{LOGS_DIR}/joint_calling/gather_vcfs_gatk.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/gather_vcfs_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.vcf})

        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" GatherVcfs \\
            {params.input_args} \\
            --OUTPUT {output.vcf} \\
            > {log} 2>&1
        
        rm -rf ${{JOB_TMP_DIR}}
        """

# S03b: Index gathered VCF with bcftools
rule gather_cohort_vcfs_index:
    input:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.unindexed.vcf.gz"
    output:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz",
        idx = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz.tbi"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/joint_calling/gather_vcfs_index.txt"
    log:
        f"{LOGS_DIR}/joint_calling/gather_vcfs_index.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        ln -f {input.vcf} {output.vcf}

        # Create index for the gathered VCF file
        bcftools index --tbi --force {output.vcf} > {log} 2>&1
        """

# Conditional VQSR rules - only included if VQSR is enabled
if config["vqsr"]["enabled"]:
    # S04a: GATK VariantRecalibrator for SNPs (optional, controlled by config)
    rule variant_recalibrator_snps:
        input:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz",
            idx = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz.tbi",
            ref = config["resolved_paths"]["reference_genome"],
            fai = config["resolved_paths"]["reference_genome"] + ".fai",
            dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict",
            hapmap = config["resolved_paths"]["vqsr_resources"]["hapmap"],
            omni = config["resolved_paths"]["vqsr_resources"]["omni"],
            snps_1000g = config["resolved_paths"]["vqsr_resources"]["1000g_snps"],
            dbsnp = config["resolved_paths"]["vqsr_resources"]["dbsnp"]
        output:
            recal = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/snps.recal",
            tranches = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/snps.tranches"
        params:
            java_opts = config["vqsr_paras"]["snp_recalibration"]["java_opts"],
            tranches_args = lambda wildcards: " ".join([f"--tranche {t}" for t in config["vqsr_paras"]["snp_recalibration"]["tranches"]]),
            annotations_args = lambda wildcards: " ".join([f"-an {a}" for a in config["vqsr_paras"]["snp_recalibration"]["annotations"]]),
            max_gaussians = config["vqsr_paras"]["snp_recalibration"]["max_gaussians"]
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt: int(config["vqsr_paras"]["snp_recalibration"]["mem"].rstrip("GgMmKk")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["java_module"],
            config["tool_modules"]["gatk_module"]
        conda:
            ENV_GATK4
        container:
            CONTAINER_GATK4
        log:
            f"{LOGS_DIR}/joint_calling/variant_recalibrator_snps.log"
        shell:
            """
            JOB_TMP_DIR="{TMP_DIR}/vqsr_snps_$$"
            mkdir -p ${{JOB_TMP_DIR}}
            mkdir -p $(dirname {output.recal})

            gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" VariantRecalibrator \\
                --reference {input.ref} \\
                --variant {input.vcf} \\
                --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \\
                --resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} \\
                --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.snps_1000g} \\
                --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 {input.dbsnp} \\
                --trust-all-polymorphic \\
                {params.tranches_args} \\
                {params.annotations_args} \\
                --mode SNP \\
                --max-gaussians {params.max_gaussians} \\
                --output {output.recal} \\
                --tranches-file {output.tranches} \\
                > {log} 2>&1
            
            EXIT_CODE=$?
            rm -rf ${{JOB_TMP_DIR}}
            exit $EXIT_CODE
            """

    # S04b: GATK VariantRecalibrator for INDELs (optional, controlled by config)
    rule variant_recalibrator_indels:
        input:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz",
            idx = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz.tbi",
            ref = config["resolved_paths"]["reference_genome"],
            fai = config["resolved_paths"]["reference_genome"] + ".fai",
            dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict",
            mills = config["resolved_paths"]["vqsr_resources"]["mills_indels"],
            dbsnp = config["resolved_paths"]["vqsr_resources"]["dbsnp"]
        output:
            recal = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/indels.recal",
            tranches = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/indels.tranches"
        params:
            java_opts = config["vqsr_paras"]["indel_recalibration"]["java_opts"],
            tranches_args = lambda wildcards: " ".join([f"--tranche {t}" for t in config["vqsr_paras"]["indel_recalibration"]["tranches"]]),
            annotations_args = lambda wildcards: " ".join([f"-an {a}" for a in config["vqsr_paras"]["indel_recalibration"]["annotations"]]),
            max_gaussians = config["vqsr_paras"]["indel_recalibration"]["max_gaussians"]
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt: int(config["vqsr_paras"]["indel_recalibration"]["mem"].rstrip("GgMmKk")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["java_module"],
            config["tool_modules"]["gatk_module"]
        conda:
            ENV_GATK4
        container:
            CONTAINER_GATK4
        log:
            f"{LOGS_DIR}/joint_calling/variant_recalibrator_indels.log"
        shell:
            """
            JOB_TMP_DIR="{TMP_DIR}/vqsr_indels_$$"
            mkdir -p ${{JOB_TMP_DIR}}
            mkdir -p $(dirname {output.recal})

            gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" VariantRecalibrator \\
                --reference {input.ref} \\
                --variant {input.vcf} \\
                --resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} \\
                --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \\
                --trust-all-polymorphic \\
                {params.tranches_args} \\
                {params.annotations_args} \\
                --mode INDEL \\
                --max-gaussians {params.max_gaussians} \\
                --output {output.recal} \\
                --tranches-file {output.tranches} \\
                > {log} 2>&1
            
            EXIT_CODE=$?
            rm -rf ${{JOB_TMP_DIR}}
            exit $EXIT_CODE
            """

    # S04c-1: GATK ApplyVQSR for SNPs (optional, controlled by config)
    rule apply_vqsr_snps_gatk:
        input:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz",
            idx = f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz.tbi",
            ref = config["resolved_paths"]["reference_genome"],
            fai = config["resolved_paths"]["reference_genome"] + ".fai",
            dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict",
            recal = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/snps.recal",
            tranches = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/snps.tranches"
        output:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.snps_vqsr.unindexed.vcf.gz"
        params:
            java_opts = config["vqsr_paras"]["applyvqsr"]["java_opts"],
            ts_filter_level = config["vqsr_paras"]["snp_recalibration"]["ts_filter_level"]
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt: int(config["vqsr_paras"]["applyvqsr"]["mem"].rstrip("GgMmKk")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["java_module"],
            config["tool_modules"]["gatk_module"]
        conda:
            ENV_GATK4
        container:
            CONTAINER_GATK4
        log:
            f"{LOGS_DIR}/joint_calling/apply_vqsr_snps_gatk.log"
        shell:
            """
            JOB_TMP_DIR="{TMP_DIR}/apply_vqsr_snps_$$"
            mkdir -p ${{JOB_TMP_DIR}}
            mkdir -p $(dirname {output.vcf})

            gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" ApplyVQSR \\
                --reference {input.ref} \\
                --variant {input.vcf} \\
                --output {output.vcf} \\
                --tranches-file {input.tranches} \\
                --recal-file {input.recal} \\
                --ts-filter-level {params.ts_filter_level} \\
                --mode SNP \\
                > {log} 2>&1
            
            rm -rf ${{JOB_TMP_DIR}}
            """

    # S04c-2: Index SNP VQSR output with bcftools
    rule apply_vqsr_snps_index:
        input:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.snps_vqsr.unindexed.vcf.gz"
        output:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.snps_vqsr.vcf.gz",
            idx = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.snps_vqsr.vcf.gz.tbi"
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
            CONTAINER_BCFTOOLS
        log:
            f"{LOGS_DIR}/joint_calling/apply_vqsr_snps_index.log"
        shell:
            """
            # Move the unindexed file to final location
            mv {input.vcf} {output.vcf}
            
            # Create index for the VQSR output VCF file
            bcftools index --tbi --force {output.vcf} > {log} 2>&1
            """

    # S04d: GATK ApplyVQSR for INDELs (optional, controlled by config)
    rule apply_vqsr_indels:
        input:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.snps_vqsr.vcf.gz",
            idx = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.snps_vqsr.vcf.gz.tbi",
            ref = config["resolved_paths"]["reference_genome"],
            fai = config["resolved_paths"]["reference_genome"] + ".fai",
            dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict",
            recal = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/indels.recal",
            tranches = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/indels.tranches"
        output:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.vqsr.unindexed.vcf.gz"
        params:
            java_opts = config["vqsr_paras"]["applyvqsr"]["java_opts"],
            ts_filter_level = config["vqsr_paras"]["indel_recalibration"]["ts_filter_level"]
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt: int(config["vqsr_paras"]["applyvqsr"]["mem"].rstrip("GgMmKk")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["java_module"],
            config["tool_modules"]["gatk_module"]
        conda:
            ENV_GATK4
        container:
            CONTAINER_GATK4
        log:
            f"{LOGS_DIR}/joint_calling/apply_vqsr_indels_gatk.log"
        shell:
            """
            JOB_TMP_DIR="{TMP_DIR}/apply_vqsr_indels_$$"
            mkdir -p ${{JOB_TMP_DIR}}
            mkdir -p $(dirname {output.vcf})

            gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" ApplyVQSR \\
                --reference {input.ref} \\
                --variant {input.vcf} \\
                --output {output.vcf} \\
                --tranches-file {input.tranches} \\
                --recal-file {input.recal} \\
                --ts-filter-level {params.ts_filter_level} \\
                --mode INDEL \\
                > {log} 2>&1
            
            rm -rf ${{JOB_TMP_DIR}}
            """

    # S04d-2: Index INDEL VQSR output with bcftools
    rule apply_vqsr_indels_index:
        input:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.vqsr.unindexed.vcf.gz"
        output:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.vqsr.vcf.gz",
            idx = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.vqsr.vcf.gz.tbi"
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
            CONTAINER_BCFTOOLS
        log:
            f"{LOGS_DIR}/joint_calling/apply_vqsr_indels_index.log"
        shell:
            """
            # Move the unindexed file to final location
            mv {input.vcf} {output.vcf}
            
            # Create index for the VQSR output VCF file
            bcftools index --tbi --force {output.vcf} > {log} 2>&1
            """

# S05: VCF normalization with bcftools norm
rule normalize_vcf:
    input:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.vqsr.vcf.gz" if config["vqsr"]["enabled"] else f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz",
        idx = f"{NGS_VARIANTS_COHORT_DIR}/vqsr/cohort.vqsr.vcf.gz.tbi" if config["vqsr"]["enabled"] else f"{NGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz.tbi",
        ref = config["resolved_paths"]["reference_genome"],
        fai = config["resolved_paths"]["reference_genome"] + ".fai"
    output:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/normalized/cohort.normalized.vcf.gz"
    params:
        options = config["vcf_processing"]["normalize"]["options"]
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/joint_calling/normalize_vcf.txt"
    log:
        f"{LOGS_DIR}/joint_calling/normalize_vcf.log"
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

# S06a: Quality filtering with vcftools (output uncompressed VCF)
rule filter_variants_vcftools:
    input:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/normalized/cohort.normalized.vcf.gz"
    output:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf"
    params:
        max_missing = config["vcf_processing"]["filtering"]["max_missing"],
        hwe_threshold = config["vcf_processing"]["filtering"]["hwe_threshold"],
        out_prefix = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered"
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
        CONTAINER_VCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/joint_calling/filter_variants_vcftools.txt"
    log:
        f"{LOGS_DIR}/joint_calling/filter_variants_vcftools.log"
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

# S06b: Compress and index filtered VCF with bcftools
rule compress_filtered_vcf:
    input:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf"
    output:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz",
        idx = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz.tbi"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/joint_calling/compress_filtered_vcf.txt"
    log:
        f"{LOGS_DIR}/joint_calling/compress_filtered_vcf.log"
    shell:
        """
        # Compress the VCF with bcftools
        bcftools view -O z -o {output.vcf} {input.vcf} > {log} 2>&1
        
        # Create index
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """

# S06c: Cohort and per-sample (-s -) summary statistics for MultiQC
rule ngs_cohort_vcf_stats:
    input:
        vcf = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz",
        idx = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz.tbi"
    output:
        stats = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.stats.txt"
    threads: 2
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BCFTOOLS
    container:
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/joint_calling/cohort_vcf_stats.txt"
    log:
        f"{LOGS_DIR}/joint_calling/cohort_vcf_stats.log"
    shell:
        """
        bcftools stats --threads {threads} -s - {input.vcf} > {output.stats} 2> {log}
        """

# Conditional phasing rules - only included if phasing is enabled
if config["phasing"]["enabled"]:
    # S07a: Phasing with SHAPEIT5 per chromosome (outputs BCF)
    rule phase_variants:
        input:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz",
            idx = f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz.tbi"
        output:
            bcf = f"{NGS_VARIANTS_COHORT_DIR}/phased/{{chromosome}}.phased.bcf"
        params:
            threads = config["phasing_paras"]["shapeit"]["threads"],
            pbwt_mac = config["phasing_paras"]["shapeit"]["pbwt_mac"],
            use_map = config["phasing"]["use_map"],
            genetic_map = lambda wildcards: resolve_input_file_path(config["phasing"]["shapeit"]["genetic_map"].format(chrom=wildcards.chromosome.replace("chr", "")), config["input_data"]["base_dir"]),
            use_reference = config["phasing"].get("use_reference", False),
            reference = lambda wildcards: resolve_input_file_path(config["phasing"]["shapeit"].get("reference", "").format(chrom=wildcards.chromosome.replace("chr", "")), config["input_data"]["base_dir"]) if config["phasing"].get("use_reference", False) else "",
        threads: lambda wildcards: config["phasing_paras"]["shapeit"]["threads"]
        resources:
            mem_mb = 16384,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        container:
            CONTAINER_SHAPEIT5
        log:
            f"{LOGS_DIR}/joint_calling/phase_variants.{{chromosome}}.log"
        shell:
            """
            mkdir -p $(dirname {output.bcf})

            # Resolve phase_common command: container uses SHAPEIT5_phase_common
            if [ -n "$SINGULARITY_CONTAINER" ] || [ -n "$container" ]; then
                PHASE_CMD="SHAPEIT5_phase_common"
            else
                source {IMPORT_SHAPEIT5} 2>/dev/null || true
                PHASE_CMD="phase_common"
            fi

            MAP_ARG=""
            if [ "{params.use_map}" = "True" ]; then
                MAP_ARG="--map {params.genetic_map}"
            fi

            PBWT_ARG=""
            if [ "{params.pbwt_mac}" != "" ]; then
                PBWT_ARG="--pbwt-mac {params.pbwt_mac}"
            fi

            REF_ARG=""
            if [ "{params.use_reference}" = "True" ] && [ -n "{params.reference}" ]; then
                REF_ARG="--reference {params.reference}"
            fi

            $PHASE_CMD \\
            --input  {input.vcf} \\
            --region {wildcards.chromosome} \\
            --output {output.bcf} \\
            --thread {params.threads} \\
            $MAP_ARG \\
            $PBWT_ARG \\
            $REF_ARG \\
            > {log} 2>&1
            """

    # S07b: Convert phased BCF to compressed VCF and create index
    rule convert_phased_bcf:
        input:
            bcf = f"{NGS_VARIANTS_COHORT_DIR}/phased/{{chromosome}}.phased.bcf"
        output:
            vcf = f"{NGS_VARIANTS_COHORT_DIR}/phased/{{chromosome}}.phased.vcf.gz",
            idx = f"{NGS_VARIANTS_COHORT_DIR}/phased/{{chromosome}}.phased.vcf.gz.tbi"
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
            CONTAINER_BCFTOOLS
        log:
            f"{LOGS_DIR}/joint_calling/convert_phased.{{chromosome}}.log"
        shell:
            """
            bcftools view -O z -o {output.vcf} {input.bcf} > {log} 2>&1
            bcftools index --tbi --force {output.vcf} >> {log} 2>&1
            """
