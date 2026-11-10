# TGS Variant Calling with DeepVariant for PacBio HiFi data
# This module performs SNP and Indel calling on PacBio HiFi alignments using DeepVariant
# Uses Singularity container for DeepVariant instead of conda

# Utility helpers
def _safe_int(value, default):
    """Return integer value if valid, otherwise fallback to default."""
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.isdigit():
            return int(stripped)
    return default


# Platform configuration
PLATFORM = config["tgs_config"]["platform"]
PLATFORM_SUFFIX = {
    "ONT": "ont",
    "PACBIO_HIFI": "pb",
    "PACBIO_CLR": "pb"
}.get(PLATFORM, "pb")

# Get directory configurations
TGS_PREPROCESSING_DIR = config["output"]["tgs"]["preprocessing"]
TGS_VARIANTS_SINGLE_DIR = config["output"]["tgs"]["variants"]["single_sample"]
TGS_VARIANTS_COHORT_DIR = config["output"]["tgs"]["variants"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# DeepVariant configuration (now under tgs_variant_calling)
DEEPVARIANT_CFG = config["tgs_variant_calling"]["deepvariant"]
DEEPVARIANT_THREADS = DEEPVARIANT_CFG["threads"]
DEEPVARIANT_MEM = DEEPVARIANT_CFG["mem"]
DEEPVARIANT_OPTIONS = DEEPVARIANT_CFG.get("options", "")
DEEPVARIANT_NUM_SHARDS = _safe_int(
    DEEPVARIANT_CFG.get("num_shards", DEEPVARIANT_THREADS),
    DEEPVARIANT_THREADS,
)
DEEPVARIANT_NUM_SHARDS = max(1, min(DEEPVARIANT_NUM_SHARDS, DEEPVARIANT_THREADS))

DEEPVARIANT_TF_INTRA_THREADS = max(1, _safe_int(DEEPVARIANT_CFG.get("tf_intra_threads", 1), 1))
DEEPVARIANT_TF_INTER_THREADS = max(1, _safe_int(DEEPVARIANT_CFG.get("tf_inter_threads", 1), 1))
DEEPVARIANT_OMP_THREADS = max(1, _safe_int(DEEPVARIANT_CFG.get("omp_threads", 1), 1))
DEEPVARIANT_MKL_THREADS = max(1, _safe_int(DEEPVARIANT_CFG.get("mkl_threads", DEEPVARIANT_OMP_THREADS), DEEPVARIANT_OMP_THREADS))
DEEPVARIANT_OPENBLAS_THREADS = max(1, _safe_int(DEEPVARIANT_CFG.get("openblas_threads", DEEPVARIANT_OMP_THREADS), DEEPVARIANT_OMP_THREADS))

# Container configuration
DEEPVARIANT_VERSION = DEEPVARIANT_CFG.get("version", "1.9.0")
# Handle empty string for container - use default if empty or missing
_container_cfg = DEEPVARIANT_CFG.get("container", "")
DEEPVARIANT_CONTAINER = _container_cfg if _container_cfg else f"docker://google/deepvariant:{DEEPVARIANT_VERSION}"

# DeepVariant model types: PACBIO for HiFi, ONT_R104 for ONT
# Automatically select model based on platform
if PLATFORM == "PACBIO_HIFI":
    DEFAULT_MODEL_TYPE = "PACBIO"
elif PLATFORM == "ONT":
    DEFAULT_MODEL_TYPE = "ONT_R104"
else:
    DEFAULT_MODEL_TYPE = "PACBIO"

# Handle empty string or missing model_type - use default based on platform
_model_type_cfg = DEEPVARIANT_CFG.get("model_type", "")
DEEPVARIANT_MODEL_TYPE = _model_type_cfg if _model_type_cfg else DEFAULT_MODEL_TYPE


# Rule: Run DeepVariant variant calling on aligned BAM files using Singularity container
rule deepvariant_variant_calling:
    input:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        bai = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam.bai",
        ref = config["reference"]["genome"],
        ref_fai = config["reference"]["genome"] + ".fai"
    output:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.vcf.gz",
        vcf_tbi = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.vcf.gz.tbi",
        gvcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.gvcf.gz",
        gvcf_tbi = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.gvcf.gz.tbi",
        visual_report = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.visual_report.html"
    params:
        model_type = DEEPVARIANT_MODEL_TYPE,
        options = DEEPVARIANT_OPTIONS,
        output_dir = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant",
        tmp_dir = f"{TMP_DIR}/{{sample}}_deepvariant",
        container = DEEPVARIANT_CONTAINER,
        num_shards = DEEPVARIANT_NUM_SHARDS,
        tf_intra = DEEPVARIANT_TF_INTRA_THREADS,
        tf_inter = DEEPVARIANT_TF_INTER_THREADS,
        omp_threads = DEEPVARIANT_OMP_THREADS,
        mkl_threads = DEEPVARIANT_MKL_THREADS,
        openblas_threads = DEEPVARIANT_OPENBLAS_THREADS
    threads: DEEPVARIANT_THREADS
    resources:
        mem_mb = lambda wildcards, attempt: _retry_mem_mb(
            deepvariant_memory_mb(DEEPVARIANT_MEM, DEEPVARIANT_NUM_SHARDS), attempt
        ),
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    container:
        DEEPVARIANT_CONTAINER
    log:
        f"{LOGS_DIR}/deepvariant/{{sample}}.deepvariant.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/deepvariant/{{sample}}.deepvariant.benchmark.txt"
    shell:
        """
        set -euo pipefail
        # Create necessary directories
        mkdir -p {params.output_dir}
        mkdir -p {params.tmp_dir}

        # Validate model type
        if [ -z "{params.model_type}" ]; then
            echo "Error: model_type is empty. Please check configuration." | tee {log}
            exit 1
        fi
        
        # Get absolute paths for container binding
        INPUT_BAM=$(readlink -f {input.bam})
        INPUT_BAI=$(readlink -f {input.bai})
        INPUT_REF=$(readlink -f {input.ref})
        OUTPUT_DIR=$(readlink -f {params.output_dir})
        TMP_DIR_ABS=$(readlink -f {params.tmp_dir})
        
        # Run DeepVariant using official run_deepvariant script
        echo "Running DeepVariant inside container..." | tee {log}
        echo "Model type: {params.model_type}" | tee -a {log}
        echo "Threads: {threads}" | tee -a {log}
        echo "Num shards: {params.num_shards}" | tee -a {log}
        echo "TF intra/inter threads: {params.tf_intra}/{params.tf_inter}" | tee -a {log}
        echo "OMP/MKL/OPENBLAS threads: {params.omp_threads}/{params.mkl_threads}/{params.openblas_threads}" | tee -a {log}

        export OMP_NUM_THREADS={params.omp_threads}
        export MKL_NUM_THREADS={params.mkl_threads}
        export OPENBLAS_NUM_THREADS={params.openblas_threads}
        export TF_NUM_INTRAOP_THREADS={params.tf_intra}
        export TF_NUM_INTEROP_THREADS={params.tf_inter}
        
        /opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model_type} \
            --ref="$INPUT_REF" \
            --reads="$INPUT_BAM" \
            --output_vcf="$OUTPUT_DIR/{wildcards.sample}.deepvariant.vcf.gz" \
            --output_gvcf="$OUTPUT_DIR/{wildcards.sample}.deepvariant.gvcf.gz" \
            --num_shards={params.num_shards} \
            --intermediate_results_dir="$TMP_DIR_ABS" \
            {params.options} \
            2>&1 | tee -a {log}
        
        # Index outputs (if not already done by run_deepvariant)
        echo "Ensuring output files are indexed..." | tee -a {log}
        if [ ! -f {output.vcf_tbi} ]; then
            tabix -p vcf {output.vcf}
        fi
        if [ ! -f {output.gvcf_tbi} ]; then
            tabix -p vcf {output.gvcf}
        fi
        
        # Generate visual report
        echo "<html><body><h1>DeepVariant Run Completed</h1>" > {output.visual_report}
        echo "<p>Sample: {wildcards.sample}</p>" >> {output.visual_report}
        echo "<p>Model Type: {params.model_type}</p>" >> {output.visual_report}
        echo "<p>Threads: {threads}</p>" >> {output.visual_report}
        echo "<p>Container: {params.container}</p>" >> {output.visual_report}
        echo "<p>Completed at: $(date)</p>" >> {output.visual_report}
        echo "</body></html>" >> {output.visual_report}
        
        # Clean up temporary directory
        echo "Cleaning up temporary files..." | tee -a {log}
        rm -rf {params.tmp_dir}
        
        echo "DeepVariant completed successfully!" | tee -a {log}
        """


FIX_DEEPVARIANT_GVCF = workflow_path("scripts/tgs_scripts/fix_deepvariant_gvcf.py")

# Rule: Prepare DeepVariant GVCF for GATK consumption by enforcing <NON_REF>
rule deepvariant_prepare_gatk_gvcf:
    input:
        gvcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.gvcf.gz"
    output:
        gvcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.gatk.gvcf.gz",
        gvcf_tbi = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.gatk.gvcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        "../envs/basic.yaml"
    container:
        config["containers"]["cyvcf2"]
    log:
        f"{LOGS_DIR}/deepvariant/{{sample}}.gvcf_nonref_fix.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/deepvariant/{{sample}}.gvcf_nonref_fix.benchmark.txt"
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.gvcf})
        mkdir -p $(dirname {log})

        TMP_TXT=$(mktemp)

        python3 {FIX_DEEPVARIANT_GVCF} \
            {input.gvcf} \
            "$TMP_TXT" \
            {log}

        bgzip -f -c "$TMP_TXT" > {output.gvcf}
        tabix -f -p vcf {output.gvcf}

        rm -f "$TMP_TXT"
        """


# Rule: Normalize and decompose DeepVariant VCF output
rule deepvariant_normalize_vcf:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.vcf.gz",
        ref = config["reference"]["genome"]
    output:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.norm.vcf.gz",
        vcf_tbi = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.norm.vcf.gz.tbi"
    threads: 2
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        "../envs/basic.yaml"
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/deepvariant/{{sample}}.normalize.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/deepvariant/{{sample}}.normalize.benchmark.txt"
    shell:
        """        
        # Normalize and split multi-allelic variants
        bcftools norm -m -both -f {input.ref} {input.vcf} -Oz -o {output.vcf} 2> {log}
        
        # Index the normalized VCF
        bcftools index -t {output.vcf} 2>> {log}
        """


# Rule: Generate statistics for DeepVariant variant calls
rule deepvariant_vcf_stats:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.norm.vcf.gz"
    output:
        stats = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/deepvariant/{{sample}}.deepvariant.stats.txt"
    threads: 1
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        "../envs/basic.yaml"
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/deepvariant/{{sample}}.stats.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/deepvariant/{{sample}}.stats.benchmark.txt"
    shell:
        """
        bcftools stats {input.vcf} > {output.stats} 2> {log}
        """


# Helper function to get TGS DeepVariant outputs
def get_tgs_deepvariant_outputs():
    """Get TGS DeepVariant outputs based on configuration"""
    if not config["workflows"].get("enable_tgs_variant_calling", False):
        return []
    
    # Only return DeepVariant outputs if tool resolves to "deepvariant"
    # Supports "Default" which auto-selects based on platform
    selected_tool = get_tgs_variant_calling_tool()
    if selected_tool != "deepvariant":
        return []
    
    outputs = []
    ont_samples = get_ont_samples()
    
    if not ont_samples:
        return []
    
    for sample in ont_samples:
        outputs.extend([
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.vcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.vcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.gvcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.gvcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.gatk.gvcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.gatk.gvcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.visual_report.html",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.norm.vcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.norm.vcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.stats.txt",
        ])
    
    return outputs
