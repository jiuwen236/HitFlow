# rules/bqsr.smk
# Snakemake rules for GATK Base Quality Score Recalibration (BQSR)
ENV_SAMTOOLS = "../envs/basic.yaml"
ENV_PICARD = "../envs/ngs/picard.yaml"
ENV_GATK4 = "../envs/ngs/gatk4.yaml"

CONTAINER_SAMTOOLS = config.get("containers", {}).get("samtools", "")
CONTAINER_PICARD = config.get("containers", {}).get("picard", "")
CONTAINER_GATK4 = config.get("containers", {}).get("gatk4", "")

# Get directory configurations from main Snakefile
NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
LOGS_DIR = config["output"]["logs"]
KNOWN_SITES_VCF = config["resolved_paths"].get("known_sites_vcf", "")

# Rule order: prefer symlink_user_bam over picard_merge_recal_bams when user provides BAM
ruleorder: symlink_user_bam > picard_merge_recal_bams

# Check for and create reference index and dictionary files if they don't exist
rule reference_fai:
    input:
        ref = config["resolved_paths"]["reference_genome"]
    output:
        fai = config["resolved_paths"]["reference_genome"] + ".fai"
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["samtools_module"]
    conda:
        ENV_SAMTOOLS
    container:
        CONTAINER_SAMTOOLS
    benchmark:
        f"{LOGS_DIR}/benchmarks/reference/create_reference_fai.txt"
    log:
        f"{LOGS_DIR}/reference/create_reference_fai.log"
    shell:
        """
        # Check if FAI already exists (to prevent overwriting)
        if [ ! -f {output.fai} ]; then
            echo "Creating reference index file {output.fai}" > {log}
            samtools faidx {input.ref} 2>> {log}
        else
            echo "Reference index file {output.fai} already exists" > {log}
        fi
        """

rule reference_dict:
    input:
        ref = config["resolved_paths"]["reference_genome"]
    output:
        dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict"
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["picard_module"],
        config["tool_modules"]["java_module"]
    conda:
        ENV_PICARD
    container:
        CONTAINER_PICARD
    benchmark:
        f"{LOGS_DIR}/benchmarks/reference/create_reference_dict.txt"
    log:
        f"{LOGS_DIR}/reference/create_reference_dict.log"
    shell:
        """
        # Check if DICT already exists (to prevent overwriting)
        if [ ! -f {output.dict} ]; then
            echo "Creating reference dictionary file {output.dict}" > {log}
            picard CreateSequenceDictionary \
                R={input.ref} \
                O={output.dict} \
                2>> {log}
        else
            echo "Reference dictionary file {output.dict} already exists" > {log}
        fi
        """

# BWA index for alignment
rule bwa_index:
    input:
        ref = config["resolved_paths"]["reference_genome"]
    output:
        amb = config["resolved_paths"]["reference_genome"] + ".amb",
        ann = config["resolved_paths"]["reference_genome"] + ".ann",
        bwt = config["resolved_paths"]["reference_genome"] + ".bwt",
        pac = config["resolved_paths"]["reference_genome"] + ".pac",
        sa = config["resolved_paths"]["reference_genome"] + ".sa"
    resources:
        mem_mb = 32768,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    threads: 4
    envmodules:
        config["tool_modules"]["bwa_module"]
    conda:
        "../envs/ngs/bwa.yaml"
    container:
        config.get("containers", {}).get("bwa_mem", "")
    benchmark:
        f"{LOGS_DIR}/benchmarks/reference/create_bwa_index.txt"
    log:
        f"{LOGS_DIR}/reference/create_bwa_index.log"
    shell:
        """
        # Check if BWA index already exists (all 5 files must exist)
        if [ ! -f {output.bwt} ] || [ ! -f {output.amb} ] || [ ! -f {output.ann} ] || [ ! -f {output.pac} ] || [ ! -f {output.sa} ]; then
            echo "Creating BWA index for {input.ref}" > {log}
            bwa index {input.ref} 2>> {log}
        else
            echo "BWA index files already exist" > {log}
        fi
        """

# Make the reference files as dependencies for BaseRecalibrator and ApplyBQSR
# S11 GATK BaseRecalibrator: Generate recalibration table
rule gatk_base_recalibrator:
    input:
        bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bam",
        bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bai",
        ref = config["resolved_paths"]["reference_genome"],
        fai = config["resolved_paths"]["reference_genome"] + ".fai",
        dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict",
        known_sites = lambda wildcards: KNOWN_SITES_VCF if KNOWN_SITES_VCF else []
    output:
        recal_report = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bqsr/{{sample_id}}.dup.sorted.recalibration_report.grp"
    params:
        java_opts = config["gatk_paras"]["java_opts_baserecalibrator"],
        known_sites_param = lambda wildcards, input: ("--known-sites " + input.known_sites) if input.known_sites else ""
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(config["gatk_paras"]["mem_baserecalibrator_job"].rstrip("GgMmKk")) * 1024,
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
        f"{LOGS_DIR}/benchmarks/gatk_bqsr/{{sample_id}}.baserecalibrator.txt"
    log:
        f"{LOGS_DIR}/gatk_bqsr/{{sample_id}}.baserecalibrator.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/gatk_br_{wildcards.sample_id}_{jobid}"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.recal_report})

        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" BaseRecalibrator \\
            -R {input.ref} \\
            -I {input.bam} \\
            {params.known_sites_param} \\
            --use-original-qualities \\
            -O {output.recal_report} \\
            > {log} 2>&1

        EXIT_CODE=$?
        exit $EXIT_CODE
        """

# S11 GATK ApplyBQSR: Apply recalibration per chromosome
rule gatk_apply_bqsr_per_chromosome:
    input:
        bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bam",
        bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bai",
        recal_report = rules.gatk_base_recalibrator.output.recal_report,
        ref = config["resolved_paths"]["reference_genome"],
        fai = config["resolved_paths"]["reference_genome"] + ".fai",
        dict = os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict"
    output:
        recal_bam_chr = f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.dup.sorted.{{chromosome}}.recal.bam",
        recal_bai_chr = f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.dup.sorted.{{chromosome}}.recal.bai"
    params:
        java_opts = config["gatk_paras"]["java_opts_applybqsr"],
        static_quals_params = " ".join([f"--static-quantized-quals {q}" for q in config["gatk_paras"]["static_quals"]]),
        extra_options = config["gatk_paras"]["applybqsr_options"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(config["gatk_paras"]["mem_applybqsr_job"].rstrip("GgMmKk")) * 1024,
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
        f"{LOGS_DIR}/benchmarks/gatk_bqsr/{{sample_id}}.{{chromosome}}.applybqsr.txt"
    log:
        f"{LOGS_DIR}/gatk_bqsr/{{sample_id}}.{{chromosome}}.applybqsr.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/gatk_abqsr_{wildcards.sample_id}_{wildcards.chromosome}_{jobid}"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.recal_bam_chr})

        gatk --java-options "{params.java_opts} -Djava.io.tmpdir=${{JOB_TMP_DIR}}" ApplyBQSR \\
            -R {input.ref} \\
            -I {input.bam} \\
            --intervals {wildcards.chromosome} \\
            -O {output.recal_bam_chr} \\
            --bqsr-recal-file {input.recal_report} \\
            {params.static_quals_params} \\
            {params.extra_options} \\
            > {log} 2>&1
        
        EXIT_CODE=$?
        exit $EXIT_CODE
        """

# S11 Picard MergeSamFiles: Merge per-chromosome recalibrated BAMs
rule picard_merge_recal_bams:
    input:
        recal_bams_chr = lambda wildcards: expand(
            f"{TMP_DIR}/{wildcards.sample_id}/alignment/{wildcards.sample_id}.dup.sorted.{{chrom}}.recal.bam", # Note: {{chrom}} here is for expand, not a rule wildcard
            sample_id=wildcards.sample_id,
            chrom=config["chromosomes"] # Use the global list of chromosomes for expansion
        ),
    output:
        merged_recal_bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.recal.bam",
        merged_recal_bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.recal.bai"
    params:
        java_opts = config["picard"]["java_opts_merge_recal_bams"],
        max_records = config["picard"]["max_records"],
        input_list_cmd = lambda wildcards, input: " ".join([f"INPUT={bam}" for bam in input.recal_bams_chr])
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(config["picard"]["mem_merge_recal_bams_job"].rstrip("GgMmKk")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["java_module"],
        config["tool_modules"]["picard_module"]
    conda:
        ENV_PICARD
    container:
        CONTAINER_PICARD
    benchmark:
        f"{LOGS_DIR}/benchmarks/picard_merge/{{sample_id}}.merge_recal_bams.txt"
    log:
        f"{LOGS_DIR}/picard_merge/{{sample_id}}.merge_recal_bams.log"
    shell:
        """
        JOB_TMP_DIR="{TMP_DIR}/picard_merge_recal_{wildcards.sample_id}_{jobid}"
        mkdir -p ${{JOB_TMP_DIR}}
        mkdir -p $(dirname {output.merged_recal_bam})

        picard MergeSamFiles \\
            VALIDATION_STRINGENCY=SILENT \\
            SORT_ORDER=coordinate \\
            CREATE_INDEX=true \\
            TMP_DIR=${{JOB_TMP_DIR}} \\
            {params.input_list_cmd} \\
            OUTPUT={output.merged_recal_bam} \\
            MAX_RECORDS_IN_RAM={params.max_records} \\
            > {log} 2>&1
        """


# Rule to create symlinks for user-provided BAM files (skip alignment and BQSR)
# This rule will be preferred over picard_merge_recal_bams when user provides a BAM file
rule symlink_user_bam:
    wildcard_constraints:
        sample_id = SYMLINK_BAM_WILDCARD_REGEX
    input:
        bam = lambda wildcards: get_sample_bam_path(wildcards.sample_id),
    output:
        bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.recal.bam",
        bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.recal.bai"
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/symlink_user_bam/{{sample_id}}.txt"
    log:
        f"{LOGS_DIR}/symlink_user_bam/{{sample_id}}.log"
    run:
        import os
        
        user_bam = input.bam
        
        # Create output directory
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        
        # Create symlink for BAM
        if os.path.islink(output.bam) or os.path.exists(output.bam):
            os.remove(output.bam)
        os.symlink(os.path.abspath(user_bam), output.bam)
        
        # Find and create symlink for BAI
        user_bai = user_bam + ".bai"
        if not os.path.exists(user_bai):
            user_bai = user_bam.replace(".bam", ".bai")
        
        if os.path.exists(user_bai):
            if os.path.islink(output.bai) or os.path.exists(output.bai):
                os.remove(output.bai)
            os.symlink(os.path.abspath(user_bai), output.bai)
        else:
            # If BAI doesn't exist, create index using samtools
            shell(f"samtools index {output.bam} {output.bai}")
        
        with open(log[0], 'w') as f:
            f.write(f"Created symlink: {output.bam} -> {user_bam}\n")
            f.write(f"BAI handled: {output.bai}\n")