# ONT SV Detection Workflow Part 2: Force calling and population analysis

CONDA_CUTESV_ENV = "../envs/long_read/ont_cutesv.yaml"
CONDA_BCFTOOLS_ENV = "../envs/basic.yaml"
CONDA_SURVIVOR_ENV = "../envs/long_read/ont_survivor.yaml"

# Platform configuration (keep consistent with tgs_sv.smk)
PLATFORM = config["tgs_config"]["platform"]
PLATFORM_SUFFIX = {
    "ONT": "ont",
    "PACBIO_HIFI": "pb",
    "PACBIO_CLR": "pb"
}.get(PLATFORM, "ont")

# SV detection tool switches (keep consistent with tgs_sv.smk)
_SV_CFG = config["tgs_config"].get("sv_detection", {})
ENABLE_MERGE_PER_SAMPLE = _SV_CFG.get("merge_per_sample", True)
ENABLE_FORCE_CALLING = _SV_CFG.get("force_calling", True)
ENABLE_POPULATION_MERGE = _SV_CFG.get("population_merge", True)

# Survivor config (population merge)
_SURV_CFG = config["tgs_config"].get("survivor", {})
SURV_POP_SUPPORT = _SURV_CFG.get("population_support", 2)
SURV_MAX_DISTANCE = _SURV_CFG.get("max_distance", 1000)
SURV_MIN_LENGTH = _SURV_CFG.get("min_length", 50)
PRE_FORCE_SUPPORT = 1  # Minimal support for pre-force population merge (union of sites)

SV_FILTER_ENABLED = config["tgs_config"].get("sv_filter", {}).get("enabled", True)

# Get directory configurations  
TGS_SV_SINGLE_DIR = config["output"]["tgs"]["sv"]["single_sample"]
TGS_SV_COHORT_DIR = config["output"]["tgs"]["sv"]["cohort"]
TGS_PREPROCESSING_DIR = config["output"]["tgs"]["preprocessing"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# S10-pre: Pre-force population merge (merge all per-sample merged SVs before force calling)
# This collects the union of variant sites across samples (support >=1) to provide a comprehensive
# candidate list for force calling. The final population merge (S11) still occurs after force calling.
rule merge_population_variants_pre_force:
    input:
        merged_vcfs = expand(f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.vcf", sample=get_ont_samples())
    output:
        merge_list = f"{TGS_SV_SINGLE_DIR}/population/pre_force/population_pre_force_merge_files.txt",
        merged_vcf = f"{TGS_SV_SINGLE_DIR}/population/pre_force/population.pre_force.merged.vcf"
    params:
        max_distance = SURV_MAX_DISTANCE,
        min_support = PRE_FORCE_SUPPORT,
        min_length = SURV_MIN_LENGTH
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["survivor_module"]
    conda:
        CONDA_SURVIVOR_ENV
    container:
        config["containers"]["survivor"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/survivor_population/population_pre_force_merge.txt"
    log:
        f"{LOGS_DIR}/survivor_population/population_pre_force_merge.log"
    shell:
        """
        mkdir -p $(dirname {output.merge_list})
        : > {output.merge_list}
        for vcf in {input.merged_vcfs}; do
            echo "$vcf" >> {output.merge_list}
        done

        # Pre-force population merge (support=1 to retain all discovered variants)
        TMP_MERGED={output.merged_vcf}.tmp
        SURVIVOR merge \
            {output.merge_list} \
            {params.max_distance} {params.min_support} 1 1 0 {params.min_length} \
            $TMP_MERGED 2> {log}
        
        # Filter out non-standard chromosomes (HLA, ALT contigs with special characters)
        # Keep only standard chromosomes: chr1-22, chrX, chrY, chrM and autosomes without 'chr' prefix
        echo "Filtering non-standard chromosomes from pre-force merged VCF..." >> {log}
        grep "^#" $TMP_MERGED > {output.merged_vcf}
        grep -v "^#" $TMP_MERGED | \
            awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/ {{print}}' \
            >> {output.merged_vcf}
        
        # Log filtering results
        BEFORE=$(grep -v "^#" $TMP_MERGED | wc -l)
        AFTER=$(grep -v "^#" {output.merged_vcf} | wc -l)
        echo "Variants before filtering: $BEFORE" >> {log}
        echo "Variants after filtering: $AFTER" >> {log}
        echo "Removed $(($BEFORE - $AFTER)) variants from non-standard chromosomes" >> {log}
        
        rm -f $TMP_MERGED
        """

# S10-pre2: Get chromosomes with read coverage (for filtering VCF before cuteFC)
rule cutesv_force_get_chroms:
    input:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        bai = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam.bai",
        input_vcf = f"{TGS_SV_SINGLE_DIR}/population/pre_force/population.pre_force.merged.vcf"
    output:
        chroms_file = f"{TMP_DIR}/{{sample}}_cutesv_force/chroms_with_coverage.txt",
        filtered_vcf = f"{TMP_DIR}/{{sample}}_cutesv_force/filtered_input.vcf"
    threads: 2
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["samtools_module"]
    conda:
        CONDA_CUTESV_ENV
    container:
        config["containers"]["samtools"]
    log:
        f"{LOGS_DIR}/cutesv_force/{{sample}}.get_chroms.log"
    shell:
        """
        mkdir -p $(dirname {output.chroms_file})
        
        echo "Getting chromosomes with read coverage..." > {log}
        
        # Get list of chromosomes with read coverage from BAM
        samtools idxstats {input.bam} | \
            awk '$3 > 0 {{print $1}}' > {output.chroms_file}
        
        # Filter input VCF to only include these chromosomes
        CHROMS=$(cat {output.chroms_file})
        if [ -n "$CHROMS" ]; then
            # Keep header lines
            grep "^#" {input.input_vcf} > {output.filtered_vcf}
            # Filter data records using awk (compatible with BusyBox)
            # Read chromosome list into awk and filter VCF records
            awk -v chroms_file="{output.chroms_file}" '
                BEGIN {{
                    while ((getline chrom < chroms_file) > 0) {{
                        chroms[chrom] = 1
                    }}
                    close(chroms_file)
                }}
                !/^#/ && ($1 in chroms) {{print}}
            ' {input.input_vcf} >> {output.filtered_vcf}
        else
            echo "Warning: No chromosomes with coverage found!" >> {log}
            cp {input.input_vcf} {output.filtered_vcf}
        fi
        
        TOTAL_VARS=$(grep -v "^#" {input.input_vcf} | wc -l || echo 0)
        FILTERED_VARS=$(grep -v "^#" {output.filtered_vcf} | wc -l || echo 0)
        echo "Total variants in input: $TOTAL_VARS" >> {log}
        echo "Variants after filtering: $FILTERED_VARS" >> {log}
        """

# S10a: cuteSV force calling with population-level pre-merged VCF using cuteFC
# NOTE: Using cuteFC for force calling since cuteSV has disabled this module
rule cutesv_force_calling:
    input:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        ref = config["reference"]["genome"],
        filtered_vcf = f"{TMP_DIR}/{{sample}}_cutesv_force/filtered_input.vcf"
    output:
        raw_vcf = f"{TMP_DIR}/{{sample}}_cutesv_force/force_called.raw.vcf"
    params:
        tmpdir = f"{TMP_DIR}/{{sample}}_cutesv_force",
        min_support = config["cutesv"]["min_support"],
        opts = config["cutesv"]["options"]
    threads: config["cutesv"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: _retry_mem_mb(
            max(memory_mb(config["cutesv"]["mem"]), 32 * 1024), attempt
        ),
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["cutesv_module"]
    conda:
        CONDA_CUTESV_ENV
    container:
        config["containers"]["cutefc"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/cutesv_force/{{sample}}.txt"
    log:
        f"{LOGS_DIR}/cutesv_force/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {output.raw_vcf})
        
        echo "Running cuteFC force calling..." > {log}

        # Clean up old temp files from previous runs to avoid conflicts
        rm -rf {params.tmpdir}/signatures {params.tmpdir}/*.pickle 2>> {log} || true
        echo "Cleaned up old temp files" >> {log}

        # Use cuteFC for force calling since cuteSV force calling is disabled
        cuteFC \
            {input.bam} \
            {input.ref} \
            {output.raw_vcf} \
            {params.tmpdir} \
            -Ivcf {input.filtered_vcf} \
            -t {threads} \
            -s {params.min_support} \
            2>> {log}

        echo "cuteFC force calling completed!" >> {log}
        """

# S10b: Post-process cuteFC output with bcftools
# Separate rule to use bcftools container for reheader, bgzip, and tabix
rule cutefc_post_process:
    input:
        raw_vcf = f"{TMP_DIR}/{{sample}}_cutesv_force/force_called.raw.vcf"
    output:
        vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/force_called/{{sample}}.force_called.vcf",
        vcf_gz = f"{TGS_SV_SINGLE_DIR}/{{sample}}/force_called/{{sample}}.force_called.vcf.gz",
        tbi = f"{TGS_SV_SINGLE_DIR}/{{sample}}/force_called/{{sample}}.force_called.vcf.gz.tbi"
    params:
        tmpdir = f"{TMP_DIR}/{{sample}}_cutesv_force"
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = "4:00:00"
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        CONDA_BCFTOOLS_ENV
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/cutesv_force/{{sample}}.post_process.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        echo "Post-processing cuteFC output..." > {log}

        SAMPLE_LIST={params.tmpdir}/force_called.samples.txt
        echo {wildcards.sample} > $SAMPLE_LIST

        # Reheader with correct sample name
        bcftools reheader --samples $SAMPLE_LIST {input.raw_vcf} | \\
        # Sort the VCF file by genomic position
        bcftools sort -Ov -o {output.vcf} 2>> {log}
        
        # Compress and index
        bgzip -f -c {output.vcf} > {output.vcf_gz} 2>> {log}
        tabix -f -p vcf {output.vcf_gz} 2>> {log}

        rm -f $SAMPLE_LIST
        echo "Post-processing completed successfully!" >> {log}
        """

# S11: Merge population-level variants using SURVIVOR
rule merge_population_variants:
    input:
        force_called_vcfs = expand(f"{TGS_SV_SINGLE_DIR}/{{sample}}/force_called/{{sample}}.force_called.vcf", 
                                 sample=get_ont_samples())
    output:
        merge_list = f"{TGS_SV_COHORT_DIR}/population_merge_files.txt",
        merged_vcf = f"{TGS_SV_COHORT_DIR}/population.merged.vcf"
    params:
        max_distance = SURV_MAX_DISTANCE,
        population_support = SURV_POP_SUPPORT,
        min_length = SURV_MIN_LENGTH
    resources:
        mem_mb = 16384,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["survivor_module"]
    conda:
        CONDA_SURVIVOR_ENV
    container:
        config["containers"]["survivor"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/survivor_population/population_merge.txt"
    log:
        f"{LOGS_DIR}/survivor_population/population_merge.log"
    shell:
        """
        mkdir -p $(dirname {output.merge_list})
        
        # Create merge file list
        for vcf in {input.force_called_vcfs}; do
            echo "$vcf" >> {output.merge_list}
        done
        
        # Merge population variants
        SURVIVOR merge \
            {output.merge_list} \
            {params.max_distance} {params.population_support} 1 1 0 {params.min_length} \
            {output.merged_vcf} 2> {log}
        """

FINAL_SV_VCF = f"{TGS_SV_SINGLE_DIR}/{{sample}}/final/{{sample}}.sv.final.vcf.gz"
FINAL_SV_TBI = f"{TGS_SV_SINGLE_DIR}/{{sample}}/final/{{sample}}.sv.final.vcf.gz.tbi"

if ENABLE_FORCE_CALLING:
    rule finalize_sv_delivery:
        input:
            force_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/force_called/{{sample}}.force_called.vcf.gz",
            force_tbi = f"{TGS_SV_SINGLE_DIR}/{{sample}}/force_called/{{sample}}.force_called.vcf.gz.tbi"
        output:
            final_vcf = FINAL_SV_VCF,
            final_tbi = FINAL_SV_TBI
        resources:
            mem_mb = 1024,
            partition = config["cluster"]["partition"],
            time = "0:30:00"
        shell:
            """
            mkdir -p $(dirname {output.final_vcf})
            ln -sf $(readlink -f {input.force_vcf}) {output.final_vcf}
            ln -sf $(readlink -f {input.force_tbi}) {output.final_tbi}
            """
else:
    if SV_FILTER_ENABLED:
        rule finalize_sv_delivery:
            input:
                filtered_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.filtered.vcf.gz",
                filtered_tbi = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.filtered.vcf.gz.tbi"
            output:
                final_vcf = FINAL_SV_VCF,
                final_tbi = FINAL_SV_TBI
            resources:
                mem_mb = 1024,
                partition = config["cluster"]["partition"],
                time = "0:30:00"
            shell:
                """
                mkdir -p $(dirname {output.final_vcf})
                ln -sf $(readlink -f {input.filtered_vcf}) {output.final_vcf}
                ln -sf $(readlink -f {input.filtered_tbi}) {output.final_tbi}
                """
    else:
        rule finalize_sv_delivery:
            input:
                merged_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.vcf"
            output:
                final_vcf = FINAL_SV_VCF,
                final_tbi = FINAL_SV_TBI
            resources:
                mem_mb = 2048,
                partition = config["cluster"]["partition"],
                time = "1:00:00"
            envmodules:
                config["tool_modules"]["bcftools_module"]
            conda:
                CONDA_CUTESV_ENV
            container:
                config["containers"]["bcftools"]
            log:
                f"{LOGS_DIR}/survivor_filter/{{sample}}.finalize.log"
            shell:
                """
                mkdir -p $(dirname {output.final_vcf})
                bcftools view -Oz {input.merged_vcf} -o {output.final_vcf} 2> {log}
                tabix -f -p vcf {output.final_vcf} 2>> {log}
                """

# Use the same target definition as rule all, including caller-only mode.
rule all_sv_detection:
    input:
        lambda w: get_sv_detection_outputs()
