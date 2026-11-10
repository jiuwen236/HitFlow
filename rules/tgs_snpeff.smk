# rules/tgs_snpeff.smk
# snpEff annotation for TGS (Third-Generation Sequencing) single-sample variants
# Supports Clair3 and DeepVariant outputs

ENV_SNPEFF = "../envs/ngs/snpeff.yaml"

import os
import shlex
from pathlib import Path

# Get directory configurations
TGS_VARIANTS_SINGLE_DIR = config["output"]["tgs"]["variants"]["single_sample"]
TGS_VARIANTS_COHORT_DIR = config["output"]["tgs"]["variants"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# Get TGS annotation configuration
TGS_ANNO_CFG = config.get("tgs_variant_annotation", {})
SNPEFF_CFG = TGS_ANNO_CFG.get("snpeff", {})

# Determine which variant calling tool is being used
# Supports "Default" which auto-selects based on platform (ONT->Clair3, PACBIO_HIFI->DeepVariant)
VARIANT_CALLING_TOOL = get_tgs_variant_calling_tool()


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


# Helper function to get input VCF path based on variant calling tool
def get_tgs_vcf_for_annotation(wildcards):
    """Get the normalized VCF path for TGS variant annotation based on the tool used."""
    sample = wildcards.sample
    tool = VARIANT_CALLING_TOOL.lower()
    
    if tool == "clair3":
        # Use Clair3 normalized output
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.norm.vcf.gz"
    elif tool == "deepvariant":
        # Use DeepVariant output (already normalized through GATK processing)
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.vcf.gz"
    else:
        raise ValueError(f"Unsupported TGS variant calling tool: {tool}")


def get_tgs_vcf_index_for_annotation(wildcards):
    """Get the VCF index path for TGS variant annotation."""
    sample = wildcards.sample
    tool = VARIANT_CALLING_TOOL.lower()
    
    if tool == "clair3":
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.norm.vcf.gz.tbi"
    elif tool == "deepvariant":
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.vcf.gz.tbi"
    else:
        raise ValueError(f"Unsupported TGS variant calling tool: {tool}")


def get_tgs_annotation_output_dir(wildcards):
    """Get the output directory for TGS annotations based on the tool used."""
    sample = wildcards.sample
    tool = VARIANT_CALLING_TOOL.lower()
    
    if tool == "clair3":
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/annotated"
    elif tool == "deepvariant":
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/annotated"
    else:
        raise ValueError(f"Unsupported TGS variant calling tool: {tool}")


# TGS-ANNO-01: snpEff annotation for TGS variants
rule tgs_snpeff_annotate:
    input:
        vcf = get_tgs_vcf_for_annotation,
        idx = get_tgs_vcf_index_for_annotation
    output:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.vcf",
        html = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff_summary.html",
        csv = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff_genes.txt"
    params:
        genome = SNPEFF_GENOME,
        options = SNPEFF_OPTIONS,
        java_opts = SNPEFF_CFG.get("java_opts", "-Xmx8G"),
        data_dir = SNPEFF_DATA_DIR,
        tool_name = VARIANT_CALLING_TOOL.lower()
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
        config["containers"]["snpeff"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.snpeff.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.snpeff.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        # Build data directory argument
        DATA_DIR_ARG=""
        if [ -n "{params.data_dir}" ]; then
            DATA_DIR_ARG="-dataDir {params.data_dir}"
        fi
        
        echo "Annotating TGS variants from {params.tool_name} with snpEff..." > {log}
        echo "Input VCF: {input.vcf}" >> {log}
        echo "Output VCF: {output.vcf}" >> {log}
        
        # Run snpEff
        snpEff {params.options} \\
            -stats {output.html} \\
            -csvStats {output.csv} \\
            {params.java_opts} \\
            $DATA_DIR_ARG \\
            -nodownload \\
            {params.genome} \\
            {input.vcf} > {output.vcf} 2>> {log}
        """


# TGS-ANNO-01b: Compress and index snpEff annotated VCF
rule tgs_snpeff_compress:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.vcf"
    output:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.vcf.gz",
        idx = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.vcf.gz.tbi"
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
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.snpeff_compress.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.snpeff_compress.log"
    shell:
        """
        bcftools view -O z -o {output.vcf} {input.vcf} > {log} 2>&1
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """


# TGS-ANNO-02: Extract snpEff annotations to TSV table
rule tgs_extract_snpeff_table:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.vcf.gz",
        idx = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.vcf.gz.tbi"
    output:
        tsv = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.tsv"
    params:
        # Extract key snpEff annotation fields
        extract_fields = "CHROM POS ID REF ALT QUAL FILTER ANN[*].ALLELE ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].HGVS_C ANN[*].HGVS_P",
        tool_name = VARIANT_CALLING_TOOL.lower()
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
        config["containers"]["snpsift"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.extract_table.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.extract_table.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        
        echo "Extracting snpEff annotations from TGS {params.tool_name} results to TSV..." > {log}
        
        # Use SnpSift extractFields to extract annotations
        SnpSift extractFields -s "," -e "." \\
            {input.vcf} \\
            {params.extract_fields} \\
            > {output.tsv} 2>> {log}
        
        # If extraction fails or produces empty output, create minimal output
        if [ ! -s {output.tsv} ]; then
            echo "Warning: No variants with snpEff annotations found" >> {log}
            echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tANN[*].ALLELE\\tANN[*].EFFECT\\tANN[*].IMPACT\\tANN[*].GENE\\tANN[*].GENEID\\tANN[*].FEATURE\\tANN[*].FEATUREID\\tANN[*].BIOTYPE\\tANN[*].HGVS_C\\tANN[*].HGVS_P" > {output.tsv}
        fi
        """


# TGS-ANNO-03: Filter variants by snpEff impact (HIGH and MODERATE) - Run SnpSift filter
rule tgs_filter_by_impact:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.vcf.gz",
        idx = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.vcf.gz.tbi"
    output:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.high_moderate.vcf"
    params:
        # Filter for HIGH and MODERATE impact variants
        impact_filter = "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')",
        java_opts = SNPEFF_CFG.get("java_opts", "-Xmx4G"),
        tool_name = VARIANT_CALLING_TOOL.lower()
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
        config["containers"]["snpsift"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.filter_impact.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.filter_impact.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        echo "Filtering TGS {params.tool_name} variants for HIGH and MODERATE impact..." > {log}
        
        # Use SnpSift filter to select high/moderate impact variants
        SnpSift filter {params.java_opts} \\
            "{params.impact_filter}" \\
            {input.vcf} \\
            > {output.vcf} 2>> {log}
        """


# TGS-ANNO-03b: Compress and index filtered impact VCF
rule tgs_compress_impact_vcf:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.high_moderate.vcf"
    output:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.high_moderate.vcf.gz",
        idx = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.snpeff.high_moderate.vcf.gz.tbi"
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
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.compress_impact.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.compress_impact.log"
    shell:
        """
        bcftools view -O z -o {output.vcf} {input.vcf} > {log} 2>&1
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """


# Helper function to get TGS annotation outputs
def get_tgs_annotation_outputs():
    """Get TGS annotation outputs based on configuration."""
    if not config.get("workflows", {}).get("enable_tgs_variant_annotation", False):
        return []
    
    # Check if TGS variant calling is enabled
    if not config.get("workflows", {}).get("enable_tgs_variant_calling", False):
        return []
    
    outputs = []
    ont_samples = get_ont_samples()
    
    if not ont_samples:
        return []
    
    tool = VARIANT_CALLING_TOOL.lower()
    
    for sample in ont_samples:
        outputs.extend([
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.snpeff.vcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.snpeff.vcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.snpeff_summary.html",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.snpeff_genes.txt",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.snpeff.tsv",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.snpeff.high_moderate.vcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.snpeff.high_moderate.vcf.gz.tbi",
        ])
    
    return outputs
