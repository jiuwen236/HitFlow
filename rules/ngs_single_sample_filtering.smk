# rules/ngs_single_sample_filtering.smk
# NGS single-sample variant filtering with hard filters (GATK Best Practices)

ENV_GATK4 = "../envs/ngs/gatk4.yaml"
ENV_BCFTOOLS = "../envs/basic.yaml"

# Container configurations
CONTAINER_GATK4 = config.get("containers", {}).get("gatk4", "")
CONTAINER_BCFTOOLS = config.get("containers", {}).get("bcftools", "")

import os

# Get directory configurations
NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
NGS_VARIANTS_SINGLE_DIR = config["output"]["ngs"]["variants"]["single_sample"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# Get reference paths
REF_GENOME = config["resolved_paths"]["reference_genome"]
REF_FAI = REF_GENOME + ".fai"
REF_DICT = os.path.splitext(REF_GENOME)[0] + ".dict"

# Helper function to get all chromosomes for a sample
def get_sample_gvcfs(wildcards):
    """Get all per-chromosome GVCF files for a sample"""
    return [f"{NGS_VARIANTS_SINGLE_DIR}/{wildcards.sample_id}/gvcf/{wildcards.sample_id}.{chrom}.hc.g.vcf.bgz"
            for chrom in config["chromosomes"]]


# S01a: Gather per-chromosome GVCFs into single GVCF (unindexed)
rule gather_sample_gvcfs_gatk:
    input:
        gvcfs = get_sample_gvcfs,
        ref = REF_GENOME,
        fai = REF_FAI,
        dict = REF_DICT
    output:
        gvcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/gvcf/{{sample_id}}.g.unindexed.vcf.gz"
    params:
        java_opts = "-Xmx8G -XX:ParallelGCThreads=1"
    threads: 1
    resources:
        mem_mb = 10240,
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.gather_gvcfs_gatk.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.gather_gvcfs_gatk.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/gather_gvcfs_{wildcards.sample_id}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.gvcf})
        
        INPUT_SYMLINK_DIR="${{JOB_TMP_DIR}}/gvcf_symlinks"
        mkdir -p "${{INPUT_SYMLINK_DIR}}"

        INPUT_ARGS=()
        for SRC in {input.gvcfs}; do
            BASENAME=$(basename "${{SRC}}" .bgz)
            LINK_PATH="${{INPUT_SYMLINK_DIR}}/${{BASENAME}}.gz"
            ln -sf "$(realpath "${{SRC}}")" "${{LINK_PATH}}"
            INPUT_ARGS+=(--INPUT "${{LINK_PATH}}")
        done

        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" GatherVcfs \
            "${{INPUT_ARGS[@]}}" \
            --OUTPUT {output.gvcf} \
            > {log} 2>&1
        
        rm -rf ${{JOB_TMP_DIR}}
        """

# S01b: Index gathered GVCF with bcftools
rule gather_sample_gvcfs:
    input:
        gvcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/gvcf/{{sample_id}}.g.unindexed.vcf.gz"
    output:
        gvcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/gvcf/{{sample_id}}.g.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/gvcf/{{sample_id}}.g.vcf.gz.tbi"
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.gather_gvcfs.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.gather_gvcfs.log"
    shell:
        """
        mkdir -p $(dirname {output.gvcf})
        ln -f {input.gvcf} {output.gvcf}

        bcftools index --tbi --force {output.gvcf} > {log} 2>&1
        """


# S02: GenotypeGVCFs - Call variants from single-sample GVCF
rule genotype_single_sample_gvcf:
    input:
        gvcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/gvcf/{{sample_id}}.g.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/gvcf/{{sample_id}}.g.vcf.gz.tbi",
        ref = REF_GENOME,
        fai = REF_FAI,
        dict = REF_DICT
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.vcf.gz"
    params:
        java_opts = config.get("ngs_single_sample_variants", {}).get("genotypegvcfs", {}).get("java_opts", "-Xmx8G -XX:ParallelGCThreads=1"),
        wes_bed = config.get("resolved_paths", {}).get("wes", {}).get("targets_bed", ""),
        wes_pad = config.get("resolved_paths", {}).get("wes", {}).get("interval_padding", 0)
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(config.get("ngs_single_sample_variants", {}).get("genotypegvcfs", {}).get("mem", "10G").rstrip("GgMmKk")) * 1024,
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.genotype_gvcf.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.genotype_gvcf.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/genotype_gvcf_{wildcards.sample_id}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.vcf})
        
        INTERVALS=""
        if [ -n "{params.wes_bed}" ]; then
            INTERVALS="-L {params.wes_bed} --interval-padding {params.wes_pad}"
        fi
        
        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" GenotypeGVCFs \\
            --reference {input.ref} \\
            --variant {input.gvcf} \\
            --output {output.vcf} \\
            $INTERVALS \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """


rule index_genotype_single_sample_gvcf:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.vcf.gz"
    output:
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.vcf.gz.tbi"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.genotype_gvcf_index.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.genotype_gvcf_index.log"
    shell:
        """
        mkdir -p $(dirname {output.idx})
        bcftools index --tbi --force {input.vcf} > {log} 2>&1
        """


# S03: Select SNPs
rule select_snps:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.vcf.gz.tbi",
        ref = REF_GENOME,
        fai = REF_FAI,
        dict = REF_DICT
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.snps.vcf.gz"
    params:
        java_opts = "-Xmx4G -XX:ParallelGCThreads=1"
    threads: 1
    resources:
        mem_mb = 5120,
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.select_snps.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.select_snps.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/select_snps_{wildcards.sample_id}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        
        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" SelectVariants \\
            --reference {input.ref} \\
            --variant {input.vcf} \\
            --select-type-to-include SNP \\
            --output {output.vcf} \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """


rule index_select_snps:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.snps.vcf.gz"
    output:
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.snps.vcf.gz.tbi"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.select_snps_index.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.select_snps_index.log"
    shell:
        """
        mkdir -p $(dirname {output.idx})
        bcftools index --tbi --force {input.vcf} > {log} 2>&1
        """


# S04: Select INDELs
rule select_indels:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.vcf.gz.tbi",
        ref = REF_GENOME,
        fai = REF_FAI,
        dict = REF_DICT
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.indels.vcf.gz"
    params:
        java_opts = "-Xmx4G -XX:ParallelGCThreads=1"
    threads: 1
    resources:
        mem_mb = 5120,
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.select_indels.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.select_indels.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/select_indels_{wildcards.sample_id}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        
        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" SelectVariants \\
            --reference {input.ref} \\
            --variant {input.vcf} \\
            --select-type-to-include INDEL \\
            --output {output.vcf} \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """


rule index_select_indels:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.indels.vcf.gz"
    output:
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.indels.vcf.gz.tbi"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.select_indels_index.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.select_indels_index.log"
    shell:
        """
        mkdir -p $(dirname {output.idx})
        bcftools index --tbi --force {input.vcf} > {log} 2>&1
        """


# S05: Hard filter SNPs (GATK Best Practices thresholds)
rule filter_snps:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.snps.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.snps.vcf.gz.tbi",
        ref = REF_GENOME,
        fai = REF_FAI,
        dict = REF_DICT
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.snps.vcf.gz"
    params:
        java_opts = "-Xmx4G -XX:ParallelGCThreads=1",
        filters = lambda wildcards: " ".join([f'--filter-expression "{f}" --filter-name "SNP_{i}"' 
                  for i, f in enumerate(config.get("ngs_single_sample_variants", {}).get("filtering", {}).get("snp_filters", [
                      "QD < 2.0", "FS > 60.0", "MQ < 40.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0"
                  ]), 1)])
    threads: 1
    resources:
        mem_mb = 5120,
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.filter_snps.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.filter_snps.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/filter_snps_{wildcards.sample_id}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.vcf})
        
        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" VariantFiltration \\
            --reference {input.ref} \\
            --variant {input.vcf} \\
            --output {output.vcf} \\
            {params.filters} \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """


rule index_filter_snps:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.snps.vcf.gz"
    output:
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.snps.vcf.gz.tbi"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.filter_snps_index.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.filter_snps_index.log"
    shell:
        """
        mkdir -p $(dirname {output.idx})
        bcftools index --tbi --force {input.vcf} > {log} 2>&1
        """


# S06: Hard filter INDELs (GATK Best Practices thresholds)
rule filter_indels:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.indels.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/raw/{{sample_id}}.raw.indels.vcf.gz.tbi",
        ref = REF_GENOME,
        fai = REF_FAI,
        dict = REF_DICT
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.indels.vcf.gz"
    params:
        java_opts = "-Xmx4G -XX:ParallelGCThreads=1",
        filters = lambda wildcards: " ".join([f'--filter-expression "{f}" --filter-name "INDEL_{i}"' 
                  for i, f in enumerate(config.get("ngs_single_sample_variants", {}).get("filtering", {}).get("indel_filters", [
                      "QD < 2.0", "FS > 200.0", "ReadPosRankSum < -20.0"
                  ]), 1)])
    threads: 1
    resources:
        mem_mb = 5120,
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.filter_indels.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.filter_indels.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/filter_indels_{wildcards.sample_id}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.vcf})
        
        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" VariantFiltration \\
            --reference {input.ref} \\
            --variant {input.vcf} \\
            --output {output.vcf} \\
            {params.filters} \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """


rule index_filter_indels:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.indels.vcf.gz"
    output:
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.indels.vcf.gz.tbi"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.filter_indels_index.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.filter_indels_index.log"
    shell:
        """
        mkdir -p $(dirname {output.idx})
        bcftools index --tbi --force {input.vcf} > {log} 2>&1
        """


# S07: Merge filtered SNPs and INDELs
rule merge_filtered_variants:
    input:
        snps = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.snps.vcf.gz",
        snps_idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.snps.vcf.gz.tbi",
        indels = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.indels.vcf.gz",
        indels_idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.indels.vcf.gz.tbi",
        ref = REF_GENOME,
        fai = REF_FAI,
        dict = REF_DICT
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.vcf.gz"
    params:
        java_opts = "-Xmx4G -XX:ParallelGCThreads=1"
    threads: 1
    resources:
        mem_mb = 5120,
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.merge_filtered.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.merge_filtered.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/merge_filtered_{wildcards.sample_id}_$$"
        mkdir -p ${{JOB_TMP_DIR}}
        
        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" MergeVcfs \\
            --INPUT {input.snps} \\
            --INPUT {input.indels} \\
            --OUTPUT {output.vcf} \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        rm -rf ${{JOB_TMP_DIR}}
        exit $EXIT_CODE
        """


rule index_merge_filtered_variants:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.vcf.gz"
    output:
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.vcf.gz.tbi"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.merge_filtered_index.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.merge_filtered_index.log"
    shell:
        """
        mkdir -p $(dirname {output.idx})
        bcftools index --tbi --force {input.vcf} > {log} 2>&1
        """


# S08: Normalize VCF with bcftools norm
rule normalize_filtered_vcf:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/filtered/{{sample_id}}.filtered.vcf.gz.tbi",
        ref = REF_GENOME,
        fai = REF_FAI
    output:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.vcf.gz.tbi"
    params:
        options = config.get("ngs_single_sample_variants", {}).get("normalize", {}).get("options", "-m -both")
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
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.normalize.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.normalize.log"
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


# S09: Summary statistics of the filtered, normalized VCF for MultiQC
rule ngs_filtered_vcf_stats:
    input:
        vcf = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.vcf.gz",
        idx = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.vcf.gz.tbi"
    output:
        stats = f"{NGS_VARIANTS_SINGLE_DIR}/{{sample_id}}/variant_calling/normalized/{{sample_id}}.filtered.normalized.stats.txt"
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
        CONTAINER_BCFTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/single_sample_filtering/{{sample_id}}.stats.tsv"
    log:
        f"{LOGS_DIR}/single_sample_filtering/{{sample_id}}.stats.log"
    shell:
        """
        bcftools stats {input.vcf} > {output.stats} 2> {log}
        """
