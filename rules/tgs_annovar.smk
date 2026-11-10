# rules/tgs_annovar.smk
# ANNOVAR annotation for TGS (Third-Generation Sequencing) single-sample variants
# Supports Clair3 and DeepVariant outputs

ENV_ANNOVAR = "../envs/ngs/annovar.yaml"

import os
from pathlib import Path

# Get directory configurations
TGS_VARIANTS_SINGLE_DIR = config["output"]["tgs"]["variants"]["single_sample"]
TGS_VARIANTS_COHORT_DIR = config["output"]["tgs"]["variants"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# Get TGS annotation configuration
TGS_ANNO_CFG = config.get("tgs_variant_annotation", {})
ANNOVAR_CFG = TGS_ANNO_CFG.get("annovar", {})

# Determine which variant calling tool is being used
# Supports "Default" which auto-selects based on platform (ONT->Clair3, PACBIO_HIFI->DeepVariant)
VARIANT_CALLING_TOOL = get_tgs_variant_calling_tool()

# ANNOVAR configuration
ANNOVAR_BUILDVER = ANNOVAR_CFG.get("buildver", "hg38")
ANNOVAR_DB_DIR = ANNOVAR_CFG.get("database_dir", "ref_data/annovar/humandb")
ANNOVAR_PROTOCOL = ANNOVAR_CFG.get("protocol", "refGene,avsnp150,clinvar_20221231,gnomad41_genome")
ANNOVAR_OPERATION = ANNOVAR_CFG.get("operation", "g,f,f,f")
ANNOVAR_OPTIONS = ANNOVAR_CFG.get("options", "-nastring . -polish")
ANNOVAR_PATH = ANNOVAR_CFG.get("annovar_path") or workflow_path("bin/annovar")


# Helper function to get input VCF path based on variant calling tool
def get_tgs_vcf_for_annovar(wildcards):
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


def get_tgs_vcf_index_for_annovar(wildcards):
    """Get the VCF index path for TGS variant annotation."""
    sample = wildcards.sample
    tool = VARIANT_CALLING_TOOL.lower()
    
    if tool == "clair3":
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.norm.vcf.gz.tbi"
    elif tool == "deepvariant":
        return f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.vcf.gz.tbi"
    else:
        raise ValueError(f"Unsupported TGS variant calling tool: {tool}")


# TGS-ANNOVAR-01: ANNOVAR annotation for TGS variants
rule tgs_annovar_annotate:
    input:
        vcf = get_tgs_vcf_for_annovar,
        idx = get_tgs_vcf_index_for_annovar
    output:
        multianno_txt = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.txt",
        multianno_vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.vcf"
    params:
        buildver = ANNOVAR_BUILDVER,
        database_dir = ANNOVAR_DB_DIR,
        protocol = ANNOVAR_PROTOCOL,
        operation = ANNOVAR_OPERATION,
        options = ANNOVAR_OPTIONS,
        annovar_path = ANNOVAR_PATH,
        output_prefix = lambda wildcards: f"{TGS_VARIANTS_SINGLE_DIR}/{wildcards.sample}/{VARIANT_CALLING_TOOL.lower()}/annotated/{wildcards.sample}.annovar",
        tool_name = VARIANT_CALLING_TOOL.lower()
    threads: lambda wildcards: ANNOVAR_CFG.get("threads", 4)
    resources:
        mem_mb = lambda wildcards, attempt: int(ANNOVAR_CFG.get("mem", "10G").rstrip("GgMmKk")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"].get("perl_module", ""),
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_ANNOVAR
    container:
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.annovar.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.annovar.log"
    shell:
        """
        mkdir -p $(dirname {output.multianno_txt})
        
        echo "Annotating TGS variants from {params.tool_name} with ANNOVAR..." > {log}
        echo "Input VCF: {input.vcf}" >> {log}
        echo "Output prefix: {params.output_prefix}" >> {log}
        echo "Build version: {params.buildver}" >> {log}
        echo "Database directory: {params.database_dir}" >> {log}
        echo "Protocol: {params.protocol}" >> {log}
        echo "Operation: {params.operation}" >> {log}
        
        # Check if ANNOVAR is installed
        if [ ! -f "{params.annovar_path}/table_annovar.pl" ]; then
            echo "ERROR: ANNOVAR not found at {params.annovar_path}/table_annovar.pl" >> {log}
            echo "Please install ANNOVAR to {params.annovar_path}/" >> {log}
            exit 1
        fi
        
        # Check if database directory exists
        if [ ! -d "{params.database_dir}" ]; then
            echo "ERROR: ANNOVAR database directory not found: {params.database_dir}" >> {log}
            echo "Please download ANNOVAR databases using annotate_variation.pl -downdb" >> {log}
            exit 1
        fi
        
        # Run table_annovar.pl
        perl {params.annovar_path}/table_annovar.pl \\
            {input.vcf} \\
            {params.database_dir} \\
            -buildver {params.buildver} \\
            -out {params.output_prefix} \\
            -remove \\
            -protocol {params.protocol} \\
            -operation {params.operation} \\
            -vcfinput \\
            -thread {threads} \\
            {params.options} \\
            >> {log} 2>&1
        
        EXIT_CODE=$?
        
        if [ $EXIT_CODE -ne 0 ]; then
            echo "ERROR: ANNOVAR annotation failed with exit code $EXIT_CODE" >> {log}
            exit $EXIT_CODE
        fi
        
        echo "ANNOVAR annotation completed successfully" >> {log}
        """


# TGS-ANNOVAR-02: Compress and index ANNOVAR VCF output
rule tgs_annovar_compress_vcf:
    input:
        vcf = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.vcf"
    output:
        vcf_gz = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.vcf.gz",
        idx = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 2048,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_ANNOVAR
    container:
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.annovar_compress.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.annovar_compress.log"
    shell:
        """
        echo "Compressing and indexing ANNOVAR VCF output..." > {log}
        
        # Compress VCF
        bcftools view -O z -o {output.vcf_gz} {input.vcf} >> {log} 2>&1
        EXIT_CODE=$?
        
        # Index VCF
        if [ $EXIT_CODE -eq 0 ]; then
            bcftools index --tbi --force {output.vcf_gz} >> {log} 2>&1
            EXIT_CODE=$?
        fi
        
        exit $EXIT_CODE
        """


# TGS-ANNOVAR-03: Extract key annotations to simplified TSV table
rule tgs_annovar_extract_table:
    input:
        txt = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.txt"
    output:
        tsv = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.summary.tsv"
    params:
        tool_name = VARIANT_CALLING_TOOL.lower()
    threads: 1
    resources:
        mem_mb = 2048,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        ENV_ANNOVAR
    container:
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.annovar_extract.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.annovar_extract.log"
    shell:
        """
        echo "Extracting key ANNOVAR annotations to summary TSV..." > {log}
        
        # Extract key columns (adjust column numbers based on your protocol)
        # Common columns: Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, 
        #                 ExonicFunc.refGene, AAChange.refGene, avsnp150, etc.
        
        if [ -s {input.txt} ]; then
            # Create simplified summary (extract first 20 columns as example)
            head -1 {input.txt} > {output.tsv}
            tail -n +2 {input.txt} >> {output.tsv}
            echo "Extracted $(tail -n +2 {output.tsv} | wc -l) variants" >> {log}
        else
            echo "Warning: Input file is empty" >> {log}
            touch {output.tsv}
        fi
        """


# TGS-ANNOVAR-04: Filter for high-impact variants (exonic, splicing, etc.)
rule tgs_annovar_filter_high_impact:
    input:
        vcf_gz = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.vcf.gz",
        idx = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.vcf.gz.tbi",
        txt = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.{ANNOVAR_BUILDVER}_multianno.txt"
    output:
        high_impact_txt = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.high_impact.txt",
        high_impact_vcf_gz = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.high_impact.vcf.gz",
        high_impact_idx = f"{TGS_VARIANTS_SINGLE_DIR}/{{sample}}/{VARIANT_CALLING_TOOL.lower()}/annotated/{{sample}}.annovar.high_impact.vcf.gz.tbi"
    params:
        tool_name = VARIANT_CALLING_TOOL.lower(),
        temp_positions = lambda wildcards: f"{TMP_DIR}/{wildcards.sample}_annovar_high_impact_positions.txt"
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_ANNOVAR
    container:
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/tgs_annotation/{{sample}}.annovar_filter.tsv"
    log:
        f"{LOGS_DIR}/tgs_annotation/{{sample}}.annovar_filter.log"
    shell:
        """
        echo "Filtering for high-impact variants..." > {log}
        
        # Filter for exonic, splicing, and other high-impact variants from TXT file
        # Func.refGene column typically contains: exonic, splicing, ncRNA_exonic, etc.
        
        head -1 {input.txt} > {output.high_impact_txt}
        
        # Filter variants (exonic, splicing, exonic;splicing)
        awk -F'\\t' 'NR>1 && ($6 ~ /exonic/ || $6 ~ /splicing/)' {input.txt} >> {output.high_impact_txt}
        
        FILTERED_COUNT=$(tail -n +2 {output.high_impact_txt} | wc -l)
        echo "Filtered $FILTERED_COUNT high-impact variants" >> {log}
        
        # Extract positions for VCF filtering
        tail -n +2 {output.high_impact_txt} | awk -F'\\t' '{{print $1":"$2"-"$3}}' > {params.temp_positions}
        
        # Filter VCF if there are high-impact variants
        if [ -s {params.temp_positions} ]; then
            bcftools view -R {params.temp_positions} -O z -o {output.high_impact_vcf_gz} {input.vcf_gz} >> {log} 2>&1
            bcftools index --tbi --force {output.high_impact_vcf_gz} >> {log} 2>&1
        else
            echo "No high-impact variants found, creating empty VCF" >> {log}
            # Create empty VCF with header only
            bcftools view -h {input.vcf_gz} | bcftools view -O z -o {output.high_impact_vcf_gz} >> {log} 2>&1
            bcftools index --tbi --force {output.high_impact_vcf_gz} >> {log} 2>&1
        fi
        
        rm -f {params.temp_positions}
        """


# Helper function to get TGS ANNOVAR annotation outputs
def get_tgs_annovar_outputs():
    """Get TGS ANNOVAR annotation outputs based on configuration."""
    if not config.get("workflows", {}).get("enable_tgs_variant_annotation", False):
        return []
    
    # Check if TGS variant calling is enabled
    if not config.get("workflows", {}).get("enable_tgs_variant_calling", False):
        return []
    
    # Check if ANNOVAR is selected as annotation tool
    annotation_tool = config.get("tgs_variant_annotation", {}).get("tool", "snpeff").strip().lower()
    if annotation_tool != "annovar":
        return []
    
    outputs = []
    ont_samples = get_ont_samples()
    
    if not ont_samples:
        return []
    
    tool = VARIANT_CALLING_TOOL.lower()
    buildver = ANNOVAR_BUILDVER
    
    for sample in ont_samples:
        outputs.extend([
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.annovar.{buildver}_multianno.txt",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.annovar.{buildver}_multianno.vcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.annovar.{buildver}_multianno.vcf.gz.tbi",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.annovar.summary.tsv",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.annovar.high_impact.txt",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.annovar.high_impact.vcf.gz",
            f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{tool}/annotated/{sample}.annovar.high_impact.vcf.gz.tbi",
        ])
    
    return outputs
