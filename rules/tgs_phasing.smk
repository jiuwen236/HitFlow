# rules/tgs_phasing.smk
# TGS cohort-level phasing using SHAPEIT5 phase_common

ENV_BCFTOOLS = "../envs/basic.yaml"

TGS_VARIANTS_SINGLE_DIR = config["output"]["tgs"]["variants"]["single_sample"]
TGS_VARIANTS_COHORT_DIR = config["output"]["tgs"]["variants"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

TGS_PHASING_CFG = config.get("tgs_joint_calling", {}).get("phasing", {})
PHASING_THREADS = TGS_PHASING_CFG.get("threads", config["phasing_paras"]["shapeit"]["threads"])
PBWT_MAC = TGS_PHASING_CFG.get("pbwt_mac", config["phasing_paras"]["shapeit"].get("pbwt_mac", 1))
USE_MAP = TGS_PHASING_CFG.get("use_map", config.get("phasing", {}).get("use_map", True))
USE_REFERENCE = TGS_PHASING_CFG.get("use_reference", config.get("phasing", {}).get("use_reference", False))
GENETIC_MAP_TEMPLATE = TGS_PHASING_CFG.get("shapeit", {}).get(
    "genetic_map",
    config.get("phasing", {}).get("shapeit", {}).get("genetic_map", "data/phasing_resources/chr{chrom}.b38.gmap.gz"),
)
REFERENCE_TEMPLATE = TGS_PHASING_CFG.get("shapeit", {}).get(
    "reference",
    config.get("phasing", {}).get("shapeit", {}).get("reference", ""),
)
CONTAINERS = config.get("containers", {})
DEFAULT_SHAPEIT5_CONTAINER_KEY = "shapeit5"


def _resolve_container_image(conf_value, default_key):
    if conf_value:
        return CONTAINERS.get(conf_value, conf_value)
    if default_key in CONTAINERS:
        return CONTAINERS[default_key]
    raise ValueError(
        f"Missing container image for SHAPEIT5 phasing. Please define containers['{default_key}'] or set tgs_joint_calling.phasing.container."
    )


PHASING_CONTAINER = _resolve_container_image(TGS_PHASING_CFG.get("container"), DEFAULT_SHAPEIT5_CONTAINER_KEY)


def _resolve_genetic_map(wildcards):
    chrom_token = wildcards.chromosome.replace("chr", "")
    map_path = GENETIC_MAP_TEMPLATE.format(chrom=chrom_token)
    return resolve_input_file_path(map_path, config["input_data"]["base_dir"])


def _resolve_reference(wildcards):
    if not REFERENCE_TEMPLATE:
        return ""
    chrom_token = wildcards.chromosome.replace("chr", "")
    ref_path = REFERENCE_TEMPLATE.format(chrom=chrom_token)
    return resolve_input_file_path(ref_path, config["input_data"]["base_dir"])


rule tgs_phase_variants:
    """Phase cohort-level TGS variants per chromosome with SHAPEIT5."""
    input:
        vcf=f"{TGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz",
        idx=f"{TGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.vcf.gz.tbi"
    output:
        bcf=f"{TGS_VARIANTS_COHORT_DIR}/phased/{{chromosome}}.phased.bcf"
    params:
        threads=PHASING_THREADS,
        pbwt_mac=PBWT_MAC,
        use_map=USE_MAP,
        genetic_map=lambda wildcards: _resolve_genetic_map(wildcards),
        use_reference=USE_REFERENCE,
        reference=lambda wildcards: _resolve_reference(wildcards)
    threads: lambda wildcards: PHASING_THREADS
    resources:
        mem_mb=16384,
        partition=config["cluster"]["partition"],
        time=config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BCFTOOLS
    container:
        PHASING_CONTAINER
    log:
        f"{LOGS_DIR}/tgs_joint_calling/phase_variants.{{chromosome}}.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_joint_calling/phase_variants.{{chromosome}}.benchmark.txt"
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

        $PHASE_CMD \
            --input {input.vcf} \
            --region {wildcards.chromosome} \
            --output {output.bcf} \
            --thread {params.threads} \
            $MAP_ARG \
            $PBWT_ARG \
            $REF_ARG \
            > {log} 2>&1
        """


rule tgs_convert_phased_bcf:
    """Convert phased BCF to compressed VCF and create index."""
    input:
        bcf=f"{TGS_VARIANTS_COHORT_DIR}/phased/{{chromosome}}.phased.bcf"
    output:
        vcf=f"{TGS_VARIANTS_COHORT_DIR}/phased/{{chromosome}}.phased.vcf.gz",
        idx=f"{TGS_VARIANTS_COHORT_DIR}/phased/{{chromosome}}.phased.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb=4096,
        partition=config["cluster"]["partition"],
        time=config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BCFTOOLS
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/tgs_joint_calling/convert_phased.{{chromosome}}.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_joint_calling/convert_phased.{{chromosome}}.benchmark.txt"
    shell:
        """
        bcftools view -O z -o {output.vcf} {input.bcf} > {log} 2>&1
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """
