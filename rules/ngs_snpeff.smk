# rules/ngs_snpeff.smk
# snpEff/SnpSift annotation for single-sample variants

ENV_SNPEFF = "../envs/ngs/snpeff.yaml"

# Container configurations
CONTAINER_SNPEFF = config.get("containers", {}).get("snpeff", "")
# The snpEff image has neither SnpSift nor bcftools, so those steps use their own images.
CONTAINER_SNPSIFT = config.get("containers", {}).get("snpsift", "") or CONTAINER_SNPEFF
CONTAINER_BCFTOOLS = config.get("containers", {}).get("bcftools", "")

import os
import shlex
from pathlib import Path

# Get directory configurations
NGS_VARIANTS_SINGLE_DIR = config["output"]["ngs"]["variants"]["single_sample"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# Get snpEff configuration
SNPEFF_CFG = config.get("ngs_single_sample_variants", {}).get("annotation", {}).get("snpeff", {})


def _sanitize_snpeff_options(raw_opts: str) -> str:
    """Remove conflicting stats/export flags while preserving extra options."""
    tokens = shlex.split(raw_opts) if raw_opts else []
    sanitized = []
    i = 0

    while i < len(tokens):
        token = tokens[i]
        if token in {"-stats", "-htmlStats", "-csvStats"}:
            # Skip companion argument when present (and not another flag).
            if i + 1 < len(tokens) and not tokens[i + 1].startswith("-"):
                i += 2
            else:
                i += 1
            continue

        sanitized.append(token)
        i += 1

    return " ".join(sanitized)


SNPEFF_OPTIONS = _sanitize_snpeff_options(SNPEFF_CFG.get("options", ""))


def _resolve_snpeff_data_dir(raw_dir: str, genome: str) -> str:
    """Normalize snpEff data directory to point at the genome folder root."""
    if not raw_dir:
        return ""

    candidate = Path(os.path.expanduser(os.path.expandvars(raw_dir)))
    candidate_resolved = candidate.resolve(strict=False)

    direct_genome = candidate_resolved / genome
    if direct_genome.exists():
        return str(candidate_resolved)

    nested_root = candidate_resolved / "data"
    if (nested_root / genome).exists():
        return str(nested_root)

    return str(candidate_resolved)


SNPEFF_GENOME = SNPEFF_CFG.get("genome", "GRCh38.99")
SNPEFF_DATA_DIR = _resolve_snpeff_data_dir(SNPEFF_CFG.get("data_dir", ""), SNPEFF_GENOME)


# S01: snpEff annotation
rule snpeff_annotate:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.vcf.gz.tbi"
    output:
        vcf = temp(f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.vcf"),
        html = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff_summary.html",
        csv = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff_genes.txt"
    params:
        genome = SNPEFF_GENOME,
        options = SNPEFF_OPTIONS,
        java_opts = SNPEFF_CFG.get("java_opts", "-Xmx8G"),
        data_dir = SNPEFF_DATA_DIR
    threads: lambda wildcards: SNPEFF_CFG.get("threads", 4)
    resources:
        mem_mb = lambda wildcards, attempt: int(SNPEFF_CFG.get("mem", "10G").rstrip("GgMmKk")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"].get("snpeff_module", "")
    conda:
        ENV_SNPEFF
    container:
        CONTAINER_SNPEFF
    benchmark:
        f"{LOGS_DIR}/benchmarks/snpeff/{{sample_id}}.snpeff.tsv"
    log:
        f"{LOGS_DIR}/snpeff/{{sample_id}}.snpeff.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})

        # Build data directory argument
        DATA_DIR_ARG=""
        if [ -n "{params.data_dir}" ]; then
            DATA_DIR_ARG="-dataDir {params.data_dir}"
        fi

        # Run snpEff (the snpEff image has no bcftools; compression is a separate rule)
        snpEff {params.options} \
            -stats {output.html} \
            -csvStats {output.csv} \
            {params.java_opts} \
            $DATA_DIR_ARG \
            -nodownload \
            {params.genome} \
            {input.vcf} > {output.vcf} 2> {log}
        """


# S01b: Compress and index the snpEff VCF
rule snpeff_compress:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.vcf"
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        "../envs/basic.yaml"
    container:
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/snpeff/{{sample_id}}.snpeff_compress.tsv"
    log:
        f"{LOGS_DIR}/snpeff/{{sample_id}}.snpeff_compress.log"
    shell:
        """
        bcftools view -O z -o {output.vcf} {input.vcf} > {log} 2>&1
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """


# S02: Extract snpEff annotations to TSV table
rule extract_snpeff_table:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.vcf.gz.tbi"
    output:
        tsv = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.tsv"
    params:
        # Extract key snpEff annotation fields
        extract_fields = "CHROM POS ID REF ALT QUAL FILTER ANN[*].ALLELE ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].HGVS_C ANN[*].HGVS_P"
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"].get("snpeff_module", "")
    conda:
        ENV_SNPEFF
    container:
        CONTAINER_SNPSIFT
    benchmark:
        f"{LOGS_DIR}/benchmarks/snpeff/{{sample_id}}.extract_table.tsv"
    log:
        f"{LOGS_DIR}/snpeff/{{sample_id}}.extract_table.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})

        # Use SnpSift extractFields to extract annotations
        SnpSift extractFields -s "," -e "." \
            {input.vcf} \
            {params.extract_fields} \
            > {output.tsv} 2> {log}

        # If extraction fails or produces empty output, create minimal output
        if [ ! -s {output.tsv} ]; then
            echo "Warning: No variants with snpEff annotations found" >> {log}
            echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tANN[*].ALLELE\\tANN[*].EFFECT\\tANN[*].IMPACT\\tANN[*].GENE\\tANN[*].GENEID\\tANN[*].FEATURE\\tANN[*].FEATUREID\\tANN[*].BIOTYPE\\tANN[*].HGVS_C\\tANN[*].HGVS_P" > {output.tsv}
        fi
        """


# S03: Filter variants by snpEff impact (optional)
rule filter_by_impact:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.vcf.gz.tbi"
    output:
        vcf = temp(f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.high_moderate.vcf")
    params:
        # Filter for HIGH and MODERATE impact variants
        impact_filter = "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')",
        java_opts = SNPEFF_CFG.get("java_opts", "-Xmx4G")
    threads: 1
    resources:
        mem_mb = 5120,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"].get("snpeff_module", "")
    conda:
        ENV_SNPEFF
    container:
        CONTAINER_SNPSIFT
    benchmark:
        f"{LOGS_DIR}/benchmarks/snpeff/{{sample_id}}.filter_impact.tsv"
    log:
        f"{LOGS_DIR}/snpeff/{{sample_id}}.filter_impact.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})

        # Use SnpSift filter to select high/moderate impact variants
        SnpSift filter {params.java_opts} \
            "{params.impact_filter}" \
            {input.vcf} > {output.vcf} 2> {log}
        """


# S03b: Compress and index the impact-filtered VCF
rule snpeff_compress_impact_vcf:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.high_moderate.vcf"
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.high_moderate.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/annotated/{{sample_id}}.snpeff.high_moderate.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        "../envs/basic.yaml"
    container:
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/snpeff/{{sample_id}}.compress_impact.tsv"
    log:
        f"{LOGS_DIR}/snpeff/{{sample_id}}.compress_impact.log"
    shell:
        """
        bcftools view -O z -o {output.vcf} {input.vcf} > {log} 2>&1
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """
