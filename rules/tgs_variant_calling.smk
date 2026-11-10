# TGS Variant Calling with Clair3 for ONT/PacBio HiFi data
# This module performs SNP and Indel calling on long-read alignments using Clair3

CONDA_CLAIR3_ENV = "../envs/long_read/ont_clair3.yaml"

# Platform configuration
PLATFORM = config["tgs_config"]["platform"]
PLATFORM_SUFFIX = {
    "ONT": "ont",
    "PACBIO_HIFI": "pb",
    "PACBIO_CLR": "pb"
}.get(PLATFORM, "ont")

# Get directory configurations
TGS_PREPROCESSING_DIR = config["output"]["tgs"]["preprocessing"]
TGS_VARIANTS_SINGLE_DIR = config["output"]["tgs"]["variants"]["single_sample"]
TGS_VARIANTS_COHORT_DIR = config["output"]["tgs"]["variants"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# Clair3 configuration
CLAIR3_CFG = config["tgs_variant_calling"]["clair3"]
CLAIR3_THREADS = CLAIR3_CFG["threads"]
CLAIR3_MEM = CLAIR3_CFG["mem"]
CLAIR3_MIN_COVERAGE = CLAIR3_CFG.get("min_coverage", 2)
CLAIR3_MIN_MQ = CLAIR3_CFG.get("min_mq", 5)
CLAIR3_CHUNK_SIZE = CLAIR3_CFG.get("chunk_size", 5000000)
CLAIR3_OPTIONS = CLAIR3_CFG.get("options", "")

# Optional contig restriction (best practice: main chromosomes)
# If enabled, will pass --ctgName to run_clair3.sh.
CLAIR3_RESTRICT_TO_MAIN_CHROMS = bool(CLAIR3_CFG.get("restrict_to_main_chroms", True))
# If user sets ctg_name explicitly, use it as-is; otherwise auto-generate.
CLAIR3_CTG_NAME = str(CLAIR3_CFG.get("ctg_name", "")).strip()

# Determine model configuration
_user_model_path = CLAIR3_CFG.get("model_path", "").strip()
_user_model_name = CLAIR3_CFG.get("model_name", "").strip()

# Default model names for each platform (bioconda included models)
_DEFAULT_MODELS = {
    "ONT": "r941_prom_sup_g5014",
    "PACBIO_HIFI": "hifi_revio",
}

# Check if using container mode
_use_container = config.get("use_container", False)

if _user_model_path:
    # User specified full model path
    CLAIR3_MODEL_ARG = _user_model_path
    CLAIR3_MODEL_NAME = ""
elif _user_model_name:
    # User specified model name
    CLAIR3_MODEL_ARG = ""
    CLAIR3_MODEL_NAME = _user_model_name
else:
    # Use default model for platform
    CLAIR3_MODEL_ARG = ""
    CLAIR3_MODEL_NAME = _DEFAULT_MODELS.get(PLATFORM, "r941_prom_sup_g5014")

# Platform-specific options for Clair3
CLAIR3_PLATFORM_FLAG = {
    "ONT": "--platform=ont",
    "PACBIO_HIFI": "--platform=hifi",
}.get(PLATFORM, "--platform=ont")


# Rule: Run Clair3 variant calling on aligned BAM files
rule clair3_variant_calling:
    input:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        bai = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam.bai",
        ref = config["reference"]["genome"],
        ref_fai = config["reference"]["genome"] + ".fai"
    output:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.vcf.gz",
        vcf_tbi = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.vcf.gz.tbi",
        gvcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.gvcf.gz",
        gvcf_tbi = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.gvcf.gz.tbi",
        done = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/.clair3_done"
    params:
        model_path_arg = CLAIR3_MODEL_ARG,
        model_name = CLAIR3_MODEL_NAME,
        platform = CLAIR3_PLATFORM_FLAG,
        min_coverage = CLAIR3_MIN_COVERAGE,
        min_mq = CLAIR3_MIN_MQ,
        chunk_size = CLAIR3_CHUNK_SIZE,
        options = CLAIR3_OPTIONS,
        restrict_to_main_chroms = CLAIR3_RESTRICT_TO_MAIN_CHROMS,
        ctg_name = CLAIR3_CTG_NAME,
        output_dir = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3",
        tmp_dir = f"{TMP_DIR}/{{sample}}_clair3",
        logs_dir = LOGS_DIR
    threads: CLAIR3_THREADS
    resources:
        mem_mb = lambda wildcards, attempt: int(CLAIR3_MEM.rstrip("g")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["clair3_module"]
    conda:
        CONDA_CLAIR3_ENV
    container:
        config["containers"]["clair3"]
    log:
        f"{LOGS_DIR}/clair3/{{sample}}.clair3.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/clair3/{{sample}}.clair3.benchmark.txt"
    shell:
        """
        # Create necessary directories
        mkdir -p {params.output_dir}
        mkdir -p {params.tmp_dir}
        
        # Determine model path
        if [ -n "{params.model_path_arg}" ]; then
            MODEL_PATH="{params.model_path_arg}"
        else
            # Check common locations for the model
            if [ -d "/opt/models/{params.model_name}" ]; then
                MODEL_PATH="/opt/models/{params.model_name}"
            elif [ -d "/usr/local/bin/models/{params.model_name}" ]; then
                MODEL_PATH="/usr/local/bin/models/{params.model_name}"
            elif [ -n "$CONDA_PREFIX" ] && [ -d "$CONDA_PREFIX/bin/models/{params.model_name}" ]; then
                MODEL_PATH="$CONDA_PREFIX/bin/models/{params.model_name}"
            else
                echo "Error: Model {params.model_name} not found in /opt/models, /usr/local/bin/models or $CONDA_PREFIX/bin/models"
                exit 1
            fi
        fi

        # Determine python path
        if [ -n "${{CONDA_PREFIX:-}}" ]; then
            PYTHON_EXEC="$CONDA_PREFIX/bin/python"
        else
            PYTHON_EXEC=$(which python)
        fi

        # Build --ctgName argument (main chromosomes) if requested
        # Auto-detect chr-prefix based on reference .fai
        CTGNAME_ARG=""
        if [ "{params.restrict_to_main_chroms}" = "True" ]; then
            if [ -n "{params.ctg_name}" ]; then
                CTGNAME_ARG="--ctg_name={params.ctg_name}"
            else
                if [ ! -f "{input.ref_fai}" ]; then
                    echo "Error: reference fai not found: {input.ref_fai}" >&2
                    exit 1
                fi

                if awk 'NR==1{{print $1}}' "{input.ref_fai}" | grep -q '^chr'; then
                    CTGNAME_ARG="--ctg_name=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
                else
                    CTGNAME_ARG="--ctg_name=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
                fi
            fi
        fi
        
        # Run Clair3
        run_clair3.sh \\
            --bam_fn={input.bam} \\
            --ref_fn={input.ref} \\
            --threads={threads} \\
            {params.platform} \\
            --model_path=$MODEL_PATH \\
            --output={params.tmp_dir} \\
            --min_coverage={params.min_coverage} \\
            --min_mq={params.min_mq} \\
            --chunk_size={params.chunk_size} \\
            --python=$PYTHON_EXEC \\
            --gvcf \\
            $CTGNAME_ARG \\
            {params.options} 2>&1 | tee {log}

        # Check for incomplete Clair3 output (insufficient reads scenario)
        MISSING=""
        for f in merge_output.vcf.gz.tbi merge_output.gvcf.gz merge_output.gvcf.gz.tbi; do
            [ ! -f "{params.tmp_dir}/$f" ] && MISSING="$MISSING $f"
        done

        if [ -n "$MISSING" ]; then
            mkdir -p {params.logs_dir}
            printf "[%s] %s | Clair3 incomplete output (insufficient reads). Missing:%s\n" \
                "$(date '+%Y-%m-%d %H:%M:%S')" "{wildcards.sample}" "$MISSING" \
                >> {params.logs_dir}/error_clair3.err
            exit 1
        fi

        # Copy output files from tmp to final directory
        cp {params.tmp_dir}/merge_output.vcf.gz {output.vcf}
        cp {params.tmp_dir}/merge_output.vcf.gz.tbi {output.vcf_tbi}
        cp {params.tmp_dir}/merge_output.gvcf.gz {output.gvcf}
        cp {params.tmp_dir}/merge_output.gvcf.gz.tbi {output.gvcf_tbi}
        
        # Create completion flag
        touch {output.done}
        
        # Clean up temporary directory (optional, can be commented out for debugging)
        # rm -rf {params.tmp_dir}
        """


# Rule: Normalize and decompose Clair3 VCF output
rule clair3_normalize_vcf:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.vcf.gz",
        ref = config["reference"]["genome"]
    output:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.norm.vcf.gz",
        vcf_tbi = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.norm.vcf.gz.tbi"
    threads: 2
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        CONDA_CLAIR3_ENV
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/clair3/{{sample}}.normalize.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/clair3/{{sample}}.normalize.benchmark.txt"
    shell:
        """        
        # Normalize and split multi-allelic variants
        bcftools norm -m -both -f {input.ref} {input.vcf} -Oz -o {output.vcf} 2> {log}
        
        # Index the normalized VCF
        bcftools index -t {output.vcf} 2>> {log}
        """


# Rule: Generate statistics for Clair3 variant calls
rule clair3_vcf_stats:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.norm.vcf.gz"
    output:
        stats = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/clair3/{{sample}}.clair3.stats.txt"
    threads: 1
    resources:
        mem_mb = 2048,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        CONDA_CLAIR3_ENV
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/clair3/{{sample}}.stats.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/clair3/{{sample}}.stats.benchmark.txt"
    shell:
        """
        bcftools stats {input.vcf} > {output.stats} 2> {log}
        """


# Helper function to get TGS variant calling outputs
def get_tgs_variant_calling_outputs():
    """Get TGS variant calling outputs based on configuration"""
    if not config["workflows"].get("enable_tgs_variant_calling", False):
        return []
    
    # Only return Clair3 outputs if tool resolves to "clair3"
    # Supports "Default" which auto-selects based on platform
    selected_tool = get_tgs_variant_calling_tool()
    if selected_tool != "clair3":
        return []
    
    outputs = []
    ont_samples = get_ont_samples()
    
    if not ont_samples:
        return []
    
    for sample in ont_samples:
        outputs.extend([
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.vcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.vcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.gvcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.gvcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.norm.vcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.norm.vcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.stats.txt",
        ])
    
    return outputs
