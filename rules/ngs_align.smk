# WGS Pipeline Snakemake workflow rules
# Extracted from main Snakefile
ENV_FASTQC = "../envs/ngs/fastqc.yaml"
ENV_PERL = "../envs/basic.yaml"
ENV_BWA = "../envs/ngs/bwa.yaml"
ENV_PICARD = "../envs/ngs/picard.yaml"

CONTAINER_FASTQC = config.get("containers", {}).get("fastqc", "")
CONTAINER_PERL = config.get("containers", {}).get("perl_json", "")
CONTAINER_BWA = config.get("containers", {}).get("bwa_mem", "")
CONTAINER_PICARD = config.get("containers", {}).get("picard", "")
CONTAINER_SAMBLASTER = config.get("containers", {}).get("samblaster", "")
CONTAINER_SAMBAMBA = config.get("containers", {}).get("sambamba", "")

# Get directory configurations from main Snakefile
NGS_PREPROCESSING_DIR = config["output"]["ngs"]["preprocessing"]
NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
TMP_DIR = config["output"]["tmp"]
LOGS_DIR = config["output"]["logs"]

# Pipe outputs group the existing container-specific rules without materializing SAM.
NGS_STREAM_ALIGNMENT = config["bwa"].get("stream_alignment", False)
if not isinstance(NGS_STREAM_ALIGNMENT, bool):
    raise ValueError("bwa.stream_alignment must be a YAML boolean")

# Mixed FASTQ/BAM samples use both alignment and symlink producers of the same BAM.
# Retain the legacy DAG rather than grouping those overlapping outputs into a pipe.
if NGS_STREAM_ALIGNMENT and any(
    sample_has_bam(sample_id) and sample_has_fastq(sample_id)
    for sample_id in get_all_sample_ids()
):
    print("[NGS] Mixed FASTQ/BAM sample detected; using disk-backed alignment.")
    NGS_STREAM_ALIGNMENT = False

# Rule order: prefer symlink_user_bam_presort over bwa_mem_step3_sort when user provides BAM
ruleorder: symlink_user_bam_presort > bwa_mem_step3_sort

# S04 FastQC - Direct processing without readset subdirectory
rule fastqc_qc:
    input:
        r1 = lambda wildcards: get_fastq_path(wildcards, "1"),
        r2 = lambda wildcards: get_fastq_path(wildcards, "2")
    output:
        html1 = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc/{{sample_id}}_1_fastqc.html",
        zip1 = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc/{{sample_id}}_1_fastqc.zip",
        html2 = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc/{{sample_id}}_2_fastqc.html",
        zip2 = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc/{{sample_id}}_2_fastqc.zip"
    params:
        outdir = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc",
        opts = config["fastqc"]["options"],
        # Convert list to space-separated string
        r1_input = lambda wildcards, input: " ".join(input.r1) if isinstance(input.r1, list) else input.r1,
        r2_input = lambda wildcards, input: " ".join(input.r2) if isinstance(input.r2, list) else input.r2
    threads: config["fastqc"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["fastqc"]["mem"].rstrip("g")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["fastqc_module"],
        config["tool_modules"]["java_module"]
    conda:
        ENV_FASTQC
    container:
        CONTAINER_FASTQC
    benchmark:
        f"{LOGS_DIR}/benchmarks/fastqc/{{sample_id}}.txt"
    log:
        f"{LOGS_DIR}/fastqc/{{sample_id}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        # FastQC can process multiple files at once
        fastqc -t {threads} {params.opts} {params.r1_input} {params.r2_input} -o {params.outdir} 2> {log}
        
        # Rename FastQC output files to match expected sample_id naming convention
        # FastQC names output based on input filename, we need to rename to sample_id format
        cd {params.outdir}
        
        # Find and rename R1 files (first fastqc output for read 1)
        for f in *_fastqc.html *_fastqc.zip; do
            if [ -f "$f" ]; then
                # Determine if it's R1 or R2 based on filename pattern
                if echo "$f" | grep -qE '_R1[_.]|_1[_.]|\\.R1\\.|_read1|_1_fastqc'; then
                    ext="${{f##*_fastqc}}"
                    if [ ! -f "{wildcards.sample_id}_1_fastqc$ext" ]; then
                        mv "$f" "{wildcards.sample_id}_1_fastqc$ext"
                    fi
                elif echo "$f" | grep -qE '_R2[_.]|_2[_.]|\\.R2\\.|_read2|_2_fastqc'; then
                    ext="${{f##*_fastqc}}"
                    if [ ! -f "{wildcards.sample_id}_2_fastqc$ext" ]; then
                        mv "$f" "{wildcards.sample_id}_2_fastqc$ext"
                    fi
                fi
            fi
        done
        """

# S05 FastQC Filter Step 1 - Simplified without readset
rule fastqc_filter_step1:
    input:
        zip1 = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc/{{sample_id}}_1_fastqc.zip",
        zip2 = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc/{{sample_id}}_2_fastqc.zip"
    output:
        txt = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc.txt"
    params:
        outdir = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}",
        ngs_scripts = config["ngs_scripts"],
        conf = config["configure_path"]
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["perl_module"]
    conda:
        ENV_PERL
    container:
        CONTAINER_PERL
    benchmark:
        f"{LOGS_DIR}/benchmarks/fastqc_filter_step1/{{sample_id}}.txt"
    log:
        f"{LOGS_DIR}/fastqc_filter_step1/{{sample_id}}.log"
    shell:
        """
        perl {params.ngs_scripts}/fastqc_step1.pl \
          -fastq1 {input.zip1} \
          -fastq2 {input.zip2} \
          -conf {params.conf} \
          -dirout {params.outdir} > {log} 2>&1
          
        sleep 1
        """

# S05 FastQC Filter Step 2 - Simplified input
rule fastqc_filter_step2:
    input:
        txt = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc.txt"
    output:
        report = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc/{{sample_id}}.fastqc.report",
        json = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc/{{sample_id}}.fastqc.report.json",
    params:
        ngs_scripts = config["ngs_scripts"],
        flag_success = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc.success",
        flag_failed = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}/fastqc.failed",
        indir = f"{NGS_PREPROCESSING_DIR}/{{sample_id}}",
        conf = config["configure_path"],
        qc_results_dir = config["output"]["ngs"]["qc_reports"],
        failed_samples_dir = config["output"]["ngs"]["qc_reports"] + "/failed_samples",
        fail_log = config["output"]["ngs"]["qc_reports"] + "/failed_samples.txt"
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["perl_module"]
    conda:
        ENV_PERL
    container:
        CONTAINER_PERL
    benchmark:
        f"{LOGS_DIR}/benchmarks/fastqc_filter_step2/{{sample_id}}.txt"
    log:
        f"{LOGS_DIR}/fastqc_filter_step2/{{sample_id}}.log"
    shell:
        """
        perl {params.ngs_scripts}/fastqc_step2.pl \
          -indir {params.indir} \
          -conf {params.conf} \
          -outfile {output.report} \
          -filename {wildcards.sample_id} > {log} 2>&1
          
        sleep 1

        mkdir -p {params.qc_results_dir}
        mkdir -p {params.failed_samples_dir}

        copy_dir() {{
            local src="$1"
            local dest="$2"
            if command -v rsync >/dev/null 2>&1; then
                rsync -a "$src" "$dest"
            else
                echo "[fastqc_filter_step2] rsync not found, fallback to cp -a." >&2
                cp -a "$src" "$dest"
            fi
        }}
        
        if [ $? -eq 0 ]; then
            echo "success" > {log}
            copy_dir {params.indir} {params.qc_results_dir}
            touch {params.flag_success}
        else
            echo "{wildcards.sample_id} fastqc '$(date)'" >> {params.fail_log}
            copy_dir {params.indir} {params.failed_samples_dir}
            touch {params.flag_failed}
            exit 1
        fi
        """

# S06-1 BWA-MEM alignment - Step 1: Alignment
rule bwa_mem_step1:
    input:
        r1 = lambda wildcards: get_input_fastq_for_alignment(wildcards, 1),
        r2 = lambda wildcards: get_input_fastq_for_alignment(wildcards, 2),
        ref = config["resolved_paths"]["reference_genome"],
        # BWA index files as dependencies
        amb = config["resolved_paths"]["reference_genome"] + ".amb",
        ann = config["resolved_paths"]["reference_genome"] + ".ann",
        bwt = config["resolved_paths"]["reference_genome"] + ".bwt",
        pac = config["resolved_paths"]["reference_genome"] + ".pac",
        sa = config["resolved_paths"]["reference_genome"] + ".sa"
    output:
        sam = (
            pipe(f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.raw.sam")
            if NGS_STREAM_ALIGNMENT
            else f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.raw.sam"
        )
    params:
        rg = lambda wildcards: f"@RG\\tID:{wildcards.sample_id}\\tSM:{wildcards.sample_id}\\tLB:wgs\\tPL:Illumina",
        opts = config["bwa"]["options"]
    threads: config["bwa"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["bwa"]["mem"].rstrip("g")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bwa_module"]
    conda:
        ENV_BWA
    container:
        CONTAINER_BWA
    benchmark:
        f"{LOGS_DIR}/benchmarks/bwa_mem/{{sample_id}}_step1.txt"
    log:
        f"{LOGS_DIR}/bwa_mem/{{sample_id}}_step1.log"
    shell:
        """
        bwa mem -t {threads} {params.opts} \
          -R '{params.rg}' \
          {input.ref} \
          {input.r1} \
          {input.r2} > {output.sam} 2> {log}
        """

# S06-2 BWA-MEM alignment - Step 2: Samblaster
rule bwa_mem_step2_samblaster:
    input:
        sam = f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.raw.sam"
    output:
        split = f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.split.sam",
        discordant = f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.discordant.sam",
        sam = (
            pipe(f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.marked.sam")
            if NGS_STREAM_ALIGNMENT
            else f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.marked.sam"
        )
    # Samblaster is single-threaded; do not reserve four concurrent pipe cores.
    threads: 1 if NGS_STREAM_ALIGNMENT else 4
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["samblaster_module"]
    conda:
        ENV_BWA
    container:
        CONTAINER_SAMBLASTER
    benchmark:
        f"{LOGS_DIR}/benchmarks/bwa_mem/{{sample_id}}_step2.txt"
    log:
        f"{LOGS_DIR}/bwa_mem/{{sample_id}}_step2.log"
    shell:
        """
        samblaster \
          -s {output.split} \
          -d {output.discordant} \
          -i {input.sam} \
          -o {output.sam} 2> {log}
        """

# S06-3 BWA-MEM alignment - Step 3: Sambamba Sort & Index
rule bwa_mem_step3_sort:
    input:
        sam = f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.marked.sam"
    output:
        bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bam",
        bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bai",
        stream_complete = (
            touch(temp(f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.streamed.complete"))
            if NGS_STREAM_ALIGNMENT else []
        )
    params:
        sort_mem = config["sambamba"]["sort_mem"]
    threads: config["bwa"]["threads"]
    resources:
        mem_mb = lambda wildcards, threads, attempt: _retry_mem_mb(
            sambamba_sort_memory_mb(config["sambamba"]["sort_mem"], threads, config["bwa"]["mem"]), attempt
        ),
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["sambamba_module"]
    conda:
        ENV_BWA
    container:
        CONTAINER_SAMBAMBA
    benchmark:
        f"{LOGS_DIR}/benchmarks/bwa_mem/{{sample_id}}_step3.txt"
    log:
        f"{LOGS_DIR}/bwa_mem/{{sample_id}}_step3.log"
    shell:
        """
        mkdir -p {NGS_ALIGNMENT_DIR}/{wildcards.sample_id}/bam
        TMP_DIR_SORT=$(mktemp -d -p .)
        trap 'rm -rf "$TMP_DIR_SORT"' EXIT

        {{
          sambamba view -S -f bam -t {threads} {input.sam} | \
          sambamba sort -m {params.sort_mem} -t {threads} \
            --tmpdir "$TMP_DIR_SORT" \
            /dev/stdin \
            -o {output.bam}
          sambamba index -t {threads} {output.bam} {output.bai}
        }} > {log} 2>&1
        """

# S07 Convert SAM to BAM and create indices - No merging needed since no readsets
rule convert_sam_to_bam:
    input:
        split_sam = f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.split.sam",
        discordant_sam = f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.discordant.sam",
        # Keep the marked-SAM pipe consumed even when a user BAM has priority.
        stream_complete = (
            f"{TMP_DIR}/{{sample_id}}/alignment/{{sample_id}}.streamed.complete"
            if NGS_STREAM_ALIGNMENT else []
        )
    output:
        split_bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.split.bam",
        split_bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.split.bai",
        discordant_bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.discordant.bam",
        discordant_bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.discordant.bai"
    params:
        mem = config["picard"]["mem"]
    threads: 1
    resources:
        mem_mb = 20480,
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
        f"{LOGS_DIR}/benchmarks/convert_sam_to_bam/{{sample_id}}.txt"
    log:
        f"{LOGS_DIR}/convert_sam_to_bam/{{sample_id}}.log"
    shell:
        """
        # Convert split SAM to BAM
        picard -Xmx{params.mem} SortSam \
          INPUT={input.split_sam} \
          OUTPUT={output.split_bam} \
          SORT_ORDER=coordinate \
          CREATE_INDEX=true \
          VALIDATION_STRINGENCY=SILENT >> {log} 2>&1
          
        # Convert discordant SAM to BAM  
        picard -Xmx{params.mem} SortSam \
          INPUT={input.discordant_sam} \
          OUTPUT={output.discordant_bam} \
          SORT_ORDER=coordinate \
          CREATE_INDEX=true \
          VALIDATION_STRINGENCY=SILENT >> {log} 2>&1
        """


# Rule to create symlinks for user-provided BAM files (skip alignment)
# This rule provides .dup.sorted.bam for QC steps (qualimap) when user provides BAM
rule symlink_user_bam_presort:
    wildcard_constraints:
        sample_id = SYMLINK_BAM_WILDCARD_REGEX
    input:
        bam = lambda wildcards: get_sample_bam_path(wildcards.sample_id),
    output:
        bam = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bam",
        bai = f"{NGS_ALIGNMENT_DIR}/{{sample_id}}/bam/{{sample_id}}.dup.sorted.bai"
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/symlink_user_bam_presort/{{sample_id}}.txt"
    log:
        f"{LOGS_DIR}/symlink_user_bam_presort/{{sample_id}}.log"
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
