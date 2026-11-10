# rules/tgs_joint_calling.smk
# Snakemake rules for TGS joint calling and population-level variant analysis
# Supports Clair3 and DeepVariant GVCF inputs
# Backends: GLnexus (recommended, fast) or GATK GenomicsDBImport+GenotypeGVCFs (legacy)
# Note: VQSR is not included as DeepVariant/Clair3 already provide quality scores

ENV_GATK4 = "../envs/ngs/gatk4.yaml"
ENV_BCFTOOLS = "../envs/basic.yaml"

import os
from pathlib import Path

# Get directory configurations from main Snakefile
TGS_VARIANTS_SINGLE_DIR = config["output"]["tgs"]["variants"]["single_sample"]
TGS_VARIANTS_COHORT_DIR = config["output"]["tgs"]["variants"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# Determine which variant calling tool is being used
# Supports "Default" which auto-selects based on platform (ONT->Clair3, PACBIO_HIFI->DeepVariant)
VARIANT_CALLING_TOOL = get_tgs_variant_calling_tool()

# Joint calling backend: "glnexus" or "gatk"
TGS_JOINT_CALLING_TOOL = config.get("tgs_joint_calling", {}).get("tool", "glnexus").strip().lower()

# Helper functions for TGS joint calling
def get_tgs_joint_samples():
    """Get list of TGS samples that have GVCF files (requires variant calling to be run first)."""
    try:
        return get_ont_samples()
    except NameError:
        return []


def get_tgs_gvcf_path(sample, chrom):
    """Get GVCF path for a TGS sample based on the variant calling tool used."""
    if VARIANT_CALLING_TOOL == "clair3":
        # Clair3 produces chromosome-specific GVCFs - need to extract per chromosome
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.gvcf.gz"
    elif VARIANT_CALLING_TOOL == "deepvariant":
        # DeepVariant GVCFs are sanitized for GATK compatibility
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.gatk.gvcf.gz"
    else:
        raise ValueError(f"Unsupported TGS variant calling tool: {VARIANT_CALLING_TOOL}")


def get_all_tgs_gvcfs_for_chrom(wildcards):
    """Return list of GVCF paths for the requested chromosome."""
    samples = get_tgs_joint_samples()
    return [get_tgs_gvcf_path(s, wildcards.chromosome) for s in samples]


def get_all_tgs_gvcfs(wildcards=None):
    """Return list of all GVCF paths for all TGS samples (used by GLnexus)."""
    samples = get_tgs_joint_samples()
    return [get_tgs_gvcf_path(s, None) for s in samples]


def _get_glnexus_preset():
    """Resolve GLnexus preset based on config or auto-detect from variant calling tool."""
    preset = config.get("tgs_joint_calling", {}).get("glnexus", {}).get("preset", "").strip()
    if preset:
        return preset
    # Auto-select based on variant calling tool
    if VARIANT_CALLING_TOOL == "deepvariant":
        return "DeepVariantWGS"
    else:
        # Clair3 GVCFs work well with DeepVariant_unfiltered preset
        return "DeepVariant_unfiltered"


# =============================================================================
# GLnexus backend: single rule replaces GenomicsDBImport + GenotypeGVCFs + GatherVcfs
# =============================================================================

if TGS_JOINT_CALLING_TOOL == "glnexus":

    # TGS-01-GLnexus: Joint genotyping with GLnexus
    rule tgs_glnexus_joint_call:
        input:
            gvcfs = get_all_tgs_gvcfs
        output:
            vcf = f"{TGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz"
        params:
            preset = _get_glnexus_preset(),
            samples = lambda wildcards: get_tgs_joint_samples()
        threads:
            config.get("tgs_joint_calling", {}).get("glnexus", {}).get("threads", 16)
        resources:
            mem_mb = lambda wildcards, attempt: int(
                config.get("tgs_joint_calling", {}).get("glnexus", {}).get("mem", "64G").rstrip("GgMmKk")
            ) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        container:
            config["containers"]["glnexus"]
        log:
            f"{LOGS_DIR}/tgs_joint_calling/glnexus_joint_call.log"
        benchmark:
            f"{LOGS_DIR}/benchmarks/tgs_joint_calling/glnexus_joint_call.benchmark.txt"
        shell:
            """
            SCRATCH_DIR="{TMP_DIR}/glnexus_scratch_$$"
            REHEADER_DIR="{TMP_DIR}/glnexus_reheader_$$"
            mkdir -p $(dirname {output.vcf}) "${{REHEADER_DIR}}"

            # Clair3 GVCFs use "SAMPLE" as sample name; reheader to unique names
            GVCF_ARGS=""
            for GVCF in {input.gvcfs}; do
                SAMPLE_NAME=$(basename $(dirname $(dirname "$GVCF")))
                REHEADERED="${{REHEADER_DIR}}/${{SAMPLE_NAME}}.gvcf.gz"
                echo "$SAMPLE_NAME" > "${{REHEADER_DIR}}/${{SAMPLE_NAME}}.txt"
                bcftools reheader -s "${{REHEADER_DIR}}/${{SAMPLE_NAME}}.txt" -o "$REHEADERED" "$GVCF"
                bcftools index -t "$REHEADERED"
                GVCF_ARGS="${{GVCF_ARGS}} ${{REHEADERED}}"
            done

            glnexus_cli \\
                --config {params.preset} \\
                --dir "${{SCRATCH_DIR}}" \\
                --threads {threads} \\
                --mem-gbytes $(({resources.mem_mb} / 1024)) \\
                ${{GVCF_ARGS}} \\
                2> {log} \\
            | bcftools view -O z -o {output.vcf}

            EXIT_CODE=$?
            rm -rf "${{SCRATCH_DIR}}" "${{REHEADER_DIR}}"
            exit $EXIT_CODE
            """


# =============================================================================
# GATK backend: GenomicsDBImport + GenotypeGVCFs + GatherVcfs (legacy, slower)
# =============================================================================

if TGS_JOINT_CALLING_TOOL == "gatk":

    # TGS-01: GATK GenomicsDBImport - Import TGS GVCFs into genomics database per chromosome
    rule tgs_genomics_db_import:
        input:
            gvcfs = get_all_tgs_gvcfs_for_chrom,
            ref = config["resolved_paths"]["reference_genome"],
            fai = config["resolved_paths"]["reference_genome"] + ".fai",
            dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict"
        output:
            db_dir = directory(f"{TGS_VARIANTS_COHORT_DIR}/genomicsdb/{{chromosome}}"),
            sample_map = f"{TGS_VARIANTS_COHORT_DIR}/sample_maps/{{chromosome}}.sample_map"
        params:
            java_opts = config["tgs_joint_calling"]["genomicsdbimport"]["java_opts"],
            batch_size = config["tgs_joint_calling"]["genomicsdbimport"]["batch_size"],
            reader_threads = config["tgs_joint_calling"]["genomicsdbimport"]["reader_threads"],
            tool = VARIANT_CALLING_TOOL,
            samples = lambda wildcards: get_tgs_joint_samples(),
            script = workflow_path("scripts/tgs_scripts/tgs_generate_sample_map.py")
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt: int(config["tgs_joint_calling"]["genomicsdbimport"]["mem"].rstrip("GgMmKk")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["java_module"],
            config["tool_modules"]["gatk_module"]
        conda:
            ENV_GATK4
        container:
            config["containers"]["gatk4"]
        log:
            f"{LOGS_DIR}/tgs_joint_calling/genomicsdb_import.{{chromosome}}.log"
        benchmark:
            f"{LOGS_DIR}/benchmarks/tgs_joint_calling/genomicsdb_import.{{chromosome}}.benchmark.txt"
        shell:
            """
            mkdir -p $(dirname {output.sample_map})

            # Generate sample map for TGS samples
            python3 {params.script} \
                --chromosome {wildcards.chromosome} \
                --output {output.sample_map} \
                --final-dir {TGS_VARIANTS_SINGLE_DIR} \
                --tool {params.tool} \
                --samples "{params.samples}"

            JOB_TMP_DIR="{TMP_DIR}/tgs_genomicsdb_{wildcards.chromosome}_$$"
            mkdir -p ${{JOB_TMP_DIR}}
            mkdir -p $(dirname {output.db_dir})

            gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" GenomicsDBImport \\
                -R {input.ref} \\
                --sample-name-map {output.sample_map} \\
                -L {wildcards.chromosome} \\
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


    # TGS-02: GATK GenotypeGVCFs - Perform joint genotyping on TGS genomics database
    rule tgs_genotype_gvcfs:
        input:
            db_dir = f"{TGS_VARIANTS_COHORT_DIR}/genomicsdb/{{chromosome}}",
            ref = config["resolved_paths"]["reference_genome"],
            fai = config["resolved_paths"]["reference_genome"] + ".fai",
            dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict"
        output:
            vcf = f"{TGS_VARIANTS_COHORT_DIR}/genotyped/{{chromosome}}.genotyped.vcf.gz"
        params:
            java_opts = config["tgs_joint_calling"]["genotypegvcfs"]["java_opts"]
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt: int(config["tgs_joint_calling"]["genotypegvcfs"]["mem"].rstrip("GgMmKk")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["java_module"],
            config["tool_modules"]["gatk_module"]
        conda:
            ENV_GATK4
        container:
            config["containers"]["gatk4"]
        log:
            f"{LOGS_DIR}/tgs_joint_calling/genotype_gvcfs.{{chromosome}}.log"
        benchmark:
            f"{LOGS_DIR}/benchmarks/tgs_joint_calling/genotype_gvcfs.{{chromosome}}.benchmark.txt"
        shell:
            """
            JOB_TMP_DIR="{TMP_DIR}/tgs_genotype_{wildcards.chromosome}_$$"
            mkdir -p ${{JOB_TMP_DIR}} $(dirname {output.vcf})

            gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" GenotypeGVCFs \\
                --reference {input.ref} \\
                --intervals {wildcards.chromosome} \\
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


    # TGS-03: GATK GatherVcfs - Gather per-chromosome VCFs into single cohort VCF
    rule tgs_gather_cohort_vcfs:
        input:
            vcfs = expand(f"{TGS_VARIANTS_COHORT_DIR}/genotyped/{{chromosome}}.genotyped.vcf.gz",
                         chromosome=config["chromosomes"])
        output:
            vcf = f"{TGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz"
        params:
            java_opts = config["tgs_joint_calling"]["gathervcfs"]["java_opts"],
            input_args = lambda wildcards, input: " ".join([f"--INPUT {vcf}" for vcf in input.vcfs])
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt: int(config["tgs_joint_calling"]["gathervcfs"]["mem"].rstrip("GgMmKk")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["java_module"],
            config["tool_modules"]["gatk_module"]
        conda:
            ENV_GATK4
        container:
            config["containers"]["gatk4"]
        log:
            f"{LOGS_DIR}/tgs_joint_calling/gather_vcfs.log"
        benchmark:
            f"{LOGS_DIR}/benchmarks/tgs_joint_calling/gather_vcfs.benchmark.txt"
        shell:
            """
            JOB_TMP_DIR="{TMP_DIR}/tgs_gather_vcfs_$$"
            mkdir -p ${{JOB_TMP_DIR}} $(dirname {output.vcf})

            gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" GatherVcfs \\
                {params.input_args} \\
                --OUTPUT {output.vcf} \\
                > {log} 2>&1

            EXIT_CODE=$?
            rm -rf ${{JOB_TMP_DIR}}
            exit $EXIT_CODE
            """


# TGS-03b: Index gathered cohort VCF with bcftools
rule tgs_index_cohort_vcf:
    input:
        vcf = f"{TGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz"
    output:
        idx = f"{TGS_VARIANTS_COHORT_DIR}/cohort/cohort.genotyped.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 2048,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BCFTOOLS
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/tgs_joint_calling/index_cohort_vcf.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_joint_calling/index_cohort_vcf.benchmark.txt"
    shell:
        """
        bcftools index --tbi --force {input.vcf} > {log} 2>&1
        """
