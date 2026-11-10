# ONT/TGS Long-read Structural Variant Detection Workflow

CONDA_ALIGN_ENV = "../envs/long_read/ont_align.yaml"
CONDA_QC_ENV = "../envs/long_read/ont_qc.yaml"
CONDA_SNIFFLES_ENV = "../envs/long_read/ont_sniffles.yaml"
CONDA_CUTESV_ENV = "../envs/long_read/ont_cutesv.yaml"
CONDA_SVIM_ENV = "../envs/long_read/ont_svim.yaml"
CONDA_NANOVAR_ENV = "../envs/long_read/ont_nanovar.yaml"
CONDA_PBSV_ENV = "../envs/long_read/pb_pbsv.yaml"
CONDA_SURVIVOR_ENV = "../envs/long_read/ont_survivor.yaml"
CONDA_TRUVARI_ENV = "../envs/long_read/ont_truvari.yaml"

# Rule priority: prefer pre-aligned BAM if available, convert unaligned BAM to FASTQ, otherwise use FASTQ directly
ruleorder: use_prealigned_bam > minimap2_from_bam > minimap2_alignment

# Platform configuration (ONT, PACBIO_HIFI, PACBIO_CLR)
PLATFORM = config["tgs_config"]["platform"]

# SV detection tool switches
_SV_CFG = config["tgs_config"].get("sv_detection", {})
ENABLE_SNIFFLES = _SV_CFG.get("sniffles", {}).get("enabled", True)
ENABLE_CUTESV = _SV_CFG.get("cutesv", {}).get("enabled", True)
ENABLE_SVIM = _SV_CFG.get("svim", {}).get("enabled", True)
ENABLE_NANOVAR = _SV_CFG.get("nanovar", {}).get("enabled", True)
ENABLE_PBSV = _SV_CFG.get("pbsv", {}).get("enabled", True)  # PacBio HiFi only
ENABLE_MERGE_PER_SAMPLE = _SV_CFG.get("merge_per_sample", True)
ENABLE_FORCE_CALLING = _SV_CFG.get("force_calling", True)
ENABLE_POPULATION_MERGE = _SV_CFG.get("population_merge", True)
# Per-sample multi-caller merge method: "truvari" (default, truvari collapse) or "survivor".
# Both produce {sample}.merged.vcf with INFO/NUM_CALLERS for filter_sv_variants.
TGS_SV_MERGE_METHOD = str(_SV_CFG.get("merge_method", "truvari")).lower()
if TGS_SV_MERGE_METHOD not in ("survivor", "truvari"):
    raise ValueError(f"tgs_config.sv_detection.merge_method must be 'survivor' or 'truvari', got '{TGS_SV_MERGE_METHOD}'")

# Platform-specific tool options
# Sniffles: use --non-germline for HiFi CCS reads (Sniffles2 auto-detects, but explicit is safer)
SNIFFLES_PLATFORM_OPTS = {
    "ONT": "",
    "PACBIO_HIFI": "",  # Sniffles2 auto-detects HiFi; no special flag needed
    "PACBIO_CLR": ""
}.get(PLATFORM, "")

# cuteSV: different parameters for different platforms
CUTESV_PLATFORM_OPTS = {
    "ONT": "--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3",
    "PACBIO_HIFI": "--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5",
    "PACBIO_CLR": "--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5"
}.get(PLATFORM, "")

# SVIM: different presets for different platforms
SVIM_PLATFORM_OPTS = {
    "ONT": "",
    "PACBIO_HIFI": "--types DEL,INS,INV,DUP:TANDEM,DUP:INT,BND --min_mapq 20",
    "PACBIO_CLR": ""
}.get(PLATFORM, "")
# File suffix for BAM naming
PLATFORM_SUFFIX = {
    "ONT": "ont",
    "PACBIO_HIFI": "pb",
    "PACBIO_CLR": "pb"
}.get(PLATFORM, "ont")
# Read group PL tag value
RG_PL = {
    "ONT": "ONT",
    "PACBIO_HIFI": "PACBIO",
    "PACBIO_CLR": "PACBIO"
}.get(PLATFORM, "ONT")
# Minimap2 preset by platform
MM2_PRESET = {
    "ONT": "-x map-ont",
    "PACBIO_HIFI": "-x map-hifi",
    "PACBIO_CLR": "-x map-pb"
}.get(PLATFORM, "-x map-ont")

# Survivor config (from tgs_config)
_SURV_CFG = config["tgs_config"].get("survivor", {})
SURV_MAX_DISTANCE = _SURV_CFG.get("max_distance", 1000)
SURV_MIN_SUPPORT = _SURV_CFG.get("min_support", 2)
SURV_MIN_LENGTH = _SURV_CFG.get("min_length", 50)
SURV_POP_SUPPORT = _SURV_CFG.get("population_support", 2)

_SV_FILTER_CFG = config["tgs_config"].get("sv_filter", {})
SV_FILTER_ENABLED = _SV_FILTER_CFG.get("enabled", True)
SV_FILTER_MIN_LENGTH = _SV_FILTER_CFG.get("min_length", max(50, SURV_MIN_LENGTH))

_COLLAPSE_CFG = config["tgs_config"].get("collapse", {})
COLLAPSE_PCTSEQ = _COLLAPSE_CFG.get("pctseq", 0)
COLLAPSE_PCTSIZE = _COLLAPSE_CFG.get("pctsize", 0.7)
COLLAPSE_REFDIST = _COLLAPSE_CFG.get("refdist", 500)
COLLAPSE_SIZEMIN = _COLLAPSE_CFG.get("sizemin", 30)
COLLAPSE_SIZEMAX = _COLLAPSE_CFG.get("sizemax", 50000)
COLLAPSE_DUP_TO_INS = _COLLAPSE_CFG.get("dup_to_ins", True)
COLLAPSE_KEEP = _COLLAPSE_CFG.get("keep", "first")
COLLAPSE_PASSONLY = _COLLAPSE_CFG.get("passonly", True)
_COLLAPSE_PRIORITY = _COLLAPSE_CFG.get("priority", ["sniffles", "cutesv", "pbsv", "svim"])
COLLAPSE_PRIORITY = ",".join(_COLLAPSE_PRIORITY) if isinstance(_COLLAPSE_PRIORITY, (list, tuple)) else str(_COLLAPSE_PRIORITY or "")
# Shared post-merge support threshold for either backend. The legacy
# tgs_config.collapse.min_callers key remains a fallback for existing configurations.
_MERGE_CFG = config["tgs_config"].get("merge", {})
TGS_MERGE_MIN_CALLERS = _MERGE_CFG.get(
    "min_callers", _COLLAPSE_CFG.get("min_callers", 2)
)

# Get directory configurations
TGS_SV_SINGLE_DIR = config["output"]["tgs"]["sv"]["single_sample"]
TGS_SV_COHORT_DIR = config["output"]["tgs"]["sv"]["cohort"]
TGS_PREPROCESSING_DIR = config["output"]["tgs"]["preprocessing"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]

# Helper: strip any existing -x preset flags from minimap2 options to avoid conflicts
_def_mm2_opts = config["minimap2"]["options"]
# naive removal of known presets to prevent duplicates
for _preset in ["-x map-ont", "-x map-hifi", "-x map-pb"]:
    _def_mm2_opts = _def_mm2_opts.replace(_preset, "").strip()


def _streamed_align_mem_mb(threads, attempt=1):
    """Allow for sorting, minimap2, and the measured 64.1 GiB ONT peak."""
    base_mb = max(
        memory_mb(config["minimap2"]["mem"]),
        samtools_sort_memory_mb(
            config["minimap2"].get("sort_mem", "3G"), threads, 24 * 1024
        ),
    )
    return _retry_mem_mb(base_mb, attempt)

# S01a: Use pre-aligned BAM file (for PACBIO_HIFI or other pre-aligned data)
# NOTE: Uses samtools container for BAM processing
rule use_prealigned_bam:
    input:
        bam = lambda wildcards: get_ont_bam_path(wildcards),
        ref = config["reference"]["genome"]
    output:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        bai = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam.bai"
    threads: 4
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["samtools_module"]
    conda:
        CONDA_ALIGN_ENV
    container:
        config["containers"]["samtools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/minimap2/{{sample}}.use_bam.txt"
    log:
        f"{LOGS_DIR}/minimap2/{{sample}}.use_bam.log"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        # Check if input BAM has sequence dictionary (@SQ lines)
        if samtools view -H {input.bam} | grep -q '^@SQ'; then
            echo "Input BAM has sequence dictionary" > {log}
            # Check if input BAM is sorted
            if samtools view -H {input.bam} | grep -q '@HD.*SO:coordinate'; then
                echo "Input BAM is sorted, creating symlink..." >> {log}
                ln -sf $(readlink -f {input.bam}) {output.bam}
            else
                echo "Input BAM is not sorted, sorting..." >> {log}
                samtools sort -@ {threads} -m 2G -o {output.bam} {input.bam} 2>> {log}
            fi
        else
            echo "Input BAM lacks sequence dictionary, adding header from reference and sorting..." >> {log}
            # Add sequence dictionary from reference and sort
            # Step 1: Extract original header (excluding @SQ lines if any)
            # Step 2: Add @SQ lines from reference
            # Step 3: Sort and save
            {{
                samtools view -H {input.bam} | grep -v '^@SQ' || true
                samtools dict {input.ref} | grep '^@SQ'
            }} > {TMP_DIR}/{{sample}}_new_header.sam
            
            samtools reheader {TMP_DIR}/{{sample}}_new_header.sam {input.bam} | \
            samtools sort -@ {threads} -m 2G -o {output.bam} - 2>> {log}
            
            rm -f {TMP_DIR}/{{sample}}_new_header.sam
        fi
        
        # Index the BAM file
        samtools index -@ {threads} {output.bam} 2>> {log}
        echo "Pre-aligned BAM ready at {output.bam}" >> {log}
        """

# S01b: Convert unaligned BAM to FASTQ (Step 1: samtools extraction)
rule minimap2_from_bam_extract:
    input:
        bam = lambda wildcards: get_ont_unaligned_bam_path(wildcards)
    output:
        fastq = f"{TMP_DIR}/{{sample}}.{PLATFORM_SUFFIX}.from_bam.fastq.gz"
    threads: 4
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["samtools_module"]
    conda:
        CONDA_ALIGN_ENV
    container:
        config["containers"]["samtools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/minimap2/{{sample}}.from_bam_extract.txt"
    log:
        f"{LOGS_DIR}/minimap2/{{sample}}.from_bam_extract.log"
    shell:
        """
        mkdir -p {TMP_DIR}
        echo "Extracting FASTQ from unaligned BAM..." > {log}
        
        samtools fastq -@ {threads} {input.bam} 2>> {log} | gzip > {output.fastq}
        
        echo "FASTQ extraction completed" >> {log}
        """

# S01b: Align extracted FASTQ and sort in one streamed step (minimap2 | samtools sort).
# Streaming avoids materializing a huge uncompressed SAM; uses a combined container.
rule minimap2_from_bam:
    input:
        fastq = f"{TMP_DIR}/{{sample}}.{PLATFORM_SUFFIX}.from_bam.fastq.gz",
        ref = config["reference"]["genome"]
    output:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        bai = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam.bai"
    params:
        read_group = f"@RG\\tID:{{sample}}\\tSM:{{sample}}\\tPL:{RG_PL}",
        opts = _def_mm2_opts,
        mm2_preset = MM2_PRESET,
        tmp_prefix = f"{TMP_DIR}/{{sample}}_sort",
        sort_mem = config["minimap2"].get("sort_mem", "3G")
    threads: config["minimap2"]["threads"]
    resources:
        mem_mb = lambda wildcards, threads, attempt: _streamed_align_mem_mb(threads, attempt),
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["minimap2_module"]
    conda:
        CONDA_ALIGN_ENV
    container:
        config["containers"]["minimap2_samtools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/minimap2/{{sample}}.from_bam.txt"
    log:
        f"{LOGS_DIR}/minimap2/{{sample}}.from_bam.log"
    shell:
        """
        mkdir -p $(dirname {output.bam}) {TMP_DIR}
        echo "Aligning extracted FASTQ (streamed minimap2 | samtools sort)..." > {log}

        minimap2 {input.ref} {input.fastq} \
            -t {threads} {params.mm2_preset} {params.opts} \
            -a \
            -R '{params.read_group}' 2>> {log} \
        | samtools sort -m {params.sort_mem} -@ {threads} \
            -T {params.tmp_prefix} \
            -o {output.bam} - 2>> {log}

        samtools index -@ {threads} {output.bam} {output.bai} 2>> {log}

        echo "Alignment + sort + index completed" >> {log}
        """

# S01c: minimap2 alignment + sort in one streamed step (minimap2 | samtools sort).
# Supports single or multiple FASTQ files (concatenated). Streaming avoids
# materializing a huge uncompressed SAM; uses a combined container.
rule minimap2_alignment:
    input:
        fastq = lambda wildcards: get_ont_fastq_path(wildcards),
        ref = config["reference"]["genome"]
    output:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        bai = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam.bai"
    params:
        read_group = f"@RG\\tID:{{sample}}\\tSM:{{sample}}\\tPL:{RG_PL}",
        opts = _def_mm2_opts,
        mm2_preset = MM2_PRESET,
        tmp_prefix = f"{TMP_DIR}/{{sample}}_sort",
        sort_mem = config["minimap2"].get("sort_mem", "3G"),
        # Convert input.fastq to space-separated string (handles both single file and list)
        fastq_input = lambda wildcards, input: " ".join(input.fastq) if isinstance(input.fastq, list) else input.fastq
    threads: config["minimap2"]["threads"]
    resources:
        mem_mb = lambda wildcards, threads, attempt: _streamed_align_mem_mb(threads, attempt),
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["minimap2_module"]
    conda:
        CONDA_ALIGN_ENV
    container:
        config["containers"]["minimap2_samtools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/minimap2/{{sample}}.align.txt"
    log:
        f"{LOGS_DIR}/minimap2/{{sample}}.align.log"
    shell:
        """
        mkdir -p $(dirname {output.bam}) {TMP_DIR}

        # Log input information
        echo "Sample: {wildcards.sample}" > {log}
        echo "Input FASTQ files: {params.fastq_input}" >> {log}

        # Streamed minimap2 | samtools sort (no intermediate SAM on disk)
        minimap2 {input.ref} {params.fastq_input} \
            -t {threads} {params.mm2_preset} {params.opts} \
            -a \
            -R '{params.read_group}' 2>> {log} \
        | samtools sort -m {params.sort_mem} -@ {threads} \
            -T {params.tmp_prefix} \
            -o {output.bam} - 2>> {log}

        samtools index -@ {threads} {output.bam} {output.bai} 2>> {log}

        echo "Alignment + sort + index completed successfully" >> {log}
        """

# S02: NanoPlot QC analysis
rule nanoplot_qc:
    input:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam"
    output:
        stats = f"{TGS_PREPROCESSING_DIR}/{{sample}}/qc/nanoplot/{{sample}}NanoPlot-report.html",
        dir_flag = f"{TGS_PREPROCESSING_DIR}/{{sample}}/qc/nanoplot/.nanoplot_done"
    params:
        outdir = f"{TGS_PREPROCESSING_DIR}/{{sample}}/qc/nanoplot",
        prefix = "{sample}",
        opts = config["nanoplot"]["options"]
    threads: config["nanoplot"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["nanoplot"]["mem"].rstrip("g")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["nanoplot_module"]
    conda:
        CONDA_QC_ENV
    container:
        config["containers"]["nanoplot"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/nanoplot/{{sample}}.txt"
    log:
        f"{LOGS_DIR}/nanoplot/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        NanoPlot {params.opts} \
            --tsv_stats \
            --info_in_report \
            -f png \
            --plots dot \
            --legacy \
            -t {threads} \
            -o {params.outdir} \
            -p {params.prefix} \
            --bam {input.bam} 2> {log}
        touch {output.dir_flag}
        """

# S04: Sniffles SV detection
rule sniffles_sv:
    input:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam"
    output:
        vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/sniffles/{{sample}}.sniffles.vcf"
    params:
        min_support = config["sniffles"]["min_support"],
        min_length = config["sniffles"]["min_length"],
        opts = config["sniffles"]["options"],
        platform_opts = SNIFFLES_PLATFORM_OPTS
    threads: config["sniffles"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["sniffles"]["mem"].rstrip("g")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["sniffles_module"]
    conda:
        CONDA_SNIFFLES_ENV
    container:
        config["containers"]["sniffles"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/sniffles/{{sample}}.txt"
    log:
        f"{LOGS_DIR}/sniffles/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        sniffles \
            --input {input.bam} \
            --vcf {output.vcf} \
            --sample-id {wildcards.sample} \
            --threads {threads} \
            --minsupport {params.min_support} \
            --minsvlen {params.min_length} \
            {params.platform_opts} {params.opts} 2> {log}
        """

# S05: cuteSV SV detection
rule cutesv_sv:
    input:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        ref = config["reference"]["genome"]
    output:
        vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/cutesv/{{sample}}.cutesv.raw.vcf"
    params:
        tmpdir = f"{TMP_DIR}/{{sample}}_cutesv",
        min_support = config["cutesv"]["min_support"],
        opts = config["cutesv"]["options"],
        platform_opts = CUTESV_PLATFORM_OPTS
    threads: config["cutesv"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["cutesv"]["mem"].rstrip("g")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["cutesv_module"]
    conda:
        CONDA_CUTESV_ENV
    container:
        config["containers"]["cutesv"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/cutesv/{{sample}}.txt"
    log:
        f"{LOGS_DIR}/cutesv/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf}) {params.tmpdir}
        
        cuteSV \
            {input.bam} \
            {input.ref} \
            {output.vcf} \
            {params.tmpdir} \
            {params.platform_opts} {params.opts} \
            --threads {threads} \
            -s {params.min_support} 2> {log}
        
        rm -rf {params.tmpdir}
        """

# S05b: Reheader cuteSV output to set correct sample name (uses bcftools container)
rule cutesv_reheader:
    input:
        vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/cutesv/{{sample}}.cutesv.raw.vcf"
    output:
        vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/cutesv/{{sample}}.cutesv.vcf"
    threads: 1
    resources:
        mem_mb = 2048,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        CONDA_CUTESV_ENV
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/cutesv/{{sample}}.reheader.log"
    shell:
        """
        # cuteSV does not support --sample parameter, reheader to set correct sample name
        TMPDIR=$(dirname {output.vcf})/tmp_reheader_$$
        mkdir -p $TMPDIR
        echo "{wildcards.sample}" > $TMPDIR/sample_name.txt
        
        bcftools reheader --samples $TMPDIR/sample_name.txt {input.vcf} > {output.vcf} 2> {log}
        
        rm -rf $TMPDIR
        """

# S06: SVIM SV detection
rule svim_sv:
    input:
        bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
        ref = config["reference"]["genome"]
    output:
        vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/svim/variants.vcf"
    params:
        outdir = f"{TGS_SV_SINGLE_DIR}/{{sample}}/svim",
        min_sv_size = config["svim"]["min_sv_size"],
        minimum_depth = config["svim"]["minimum_depth"],
        opts = config["svim"]["options"],
        platform_opts = SVIM_PLATFORM_OPTS
    threads: config["svim"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: int(config["svim"]["mem"].rstrip("g")) * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["svim_module"]
    conda:
        CONDA_SVIM_ENV
    container:
        config["containers"]["svim"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/svim/{{sample}}.txt"
    log:
        f"{LOGS_DIR}/svim/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        svim alignment \
            {params.outdir} \
            {input.bam} \
            {input.ref} \
            --sample {wildcards.sample} \
            {params.platform_opts} {params.opts} \
            --min_sv_size {params.min_sv_size} \
            --minimum_depth {params.minimum_depth} 2> {log}
        """

# S07: NanoVar SV detection (ONT only)
if PLATFORM == "ONT":
    rule nanovar_sv:
        input:
            bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
            ref = config["reference"]["genome"]
        output:
            vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/nanovar/{{sample}}.nanovar.pass.vcf"
        params:
            outdir = f"{TGS_SV_SINGLE_DIR}/{{sample}}/nanovar",
            opts = config["nanovar"]["options"]
        threads: config["nanovar"]["threads"]
        resources:
            mem_mb = lambda wildcards, attempt: int(config["nanovar"]["mem"].rstrip("g")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["nanovar_module"]
        conda:
            CONDA_NANOVAR_ENV
        container:
            config["containers"]["nanovar"]
        log:
            f"{LOGS_DIR}/nanovar/{{sample}}.log"
        shell:
            """
            mkdir -p {params.outdir}
            nanovar {params.opts} \
                -t {threads} \
                {input.bam} \
                {input.ref} \
                {params.outdir} 2> {log}
            # NanoVar generates output based on BAM filename, so we need to rename it
            if [ -f "{params.outdir}/{wildcards.sample}.{PLATFORM_SUFFIX}.sorted.nanovar.pass.vcf" ]; then
                mv "{params.outdir}/{wildcards.sample}.{PLATFORM_SUFFIX}.sorted.nanovar.pass.vcf" "{output.vcf}"
            fi
            """

# S08: PBSV SV detection (PacBio HiFi only)
# pbsv is a two-step process: discover (BAM -> svsig.gz) and call (svsig.gz -> VCF)
if PLATFORM == "PACBIO_HIFI" and ENABLE_PBSV:
    rule pbsv_discover:
        input:
            bam = f"{TGS_PREPROCESSING_DIR}/{{sample}}/alignment/{{sample}}.{PLATFORM_SUFFIX}.sorted.bam",
            ref = config["reference"]["genome"]
        output:
            svsig = f"{TGS_SV_SINGLE_DIR}/{{sample}}/pbsv/{{sample}}.svsig.gz"
        params:
            opts = config["pbsv"]["discover_options"]
        threads: config["pbsv"]["threads"]
        resources:
            mem_mb = lambda wildcards, attempt: int(config["pbsv"]["mem"].rstrip("g")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["pbsv_module"]
        conda:
            CONDA_PBSV_ENV
        container:
            config["containers"]["pbsv"]
        log:
            f"{LOGS_DIR}/pbsv/{{sample}}_discover.log"
        shell:
            """
            mkdir -p $(dirname {output.svsig})
            pbsv discover \
                {params.opts} \
                {input.bam} \
                {output.svsig} 2> {log}
            """

    rule pbsv_call:
        input:
            svsig = f"{TGS_SV_SINGLE_DIR}/{{sample}}/pbsv/{{sample}}.svsig.gz",
            ref = config["reference"]["genome"]
        output:
            vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/pbsv/{{sample}}.pbsv.raw.vcf"
        params:
            opts = config["pbsv"]["call_options"],
            min_sv_length = config["pbsv"]["min_sv_length"]
        threads: config["pbsv"]["threads"]
        resources:
            mem_mb = lambda wildcards, attempt: int(config["pbsv"]["mem"].rstrip("g")) * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["pbsv_module"]
        conda:
            CONDA_PBSV_ENV
        container:
            config["containers"]["pbsv"]
        log:
            f"{LOGS_DIR}/pbsv/{{sample}}_call.log"
        shell:
            """
            pbsv call \
                {params.opts} \
                -j {threads} \
                --min-sv-length {params.min_sv_length} \
                {input.ref} \
                {input.svsig} \
                {output.vcf} 2> {log}
            """

    rule pbsv_reheader:
        input:
            vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/pbsv/{{sample}}.pbsv.raw.vcf"
        output:
            vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/pbsv/{{sample}}.pbsv.vcf"
        threads: 1
        resources:
            mem_mb = 2048,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["bcftools_module"]
        conda:
            CONDA_PBSV_ENV
        container:
            config["containers"]["bcftools"]
        log:
            f"{LOGS_DIR}/pbsv/{{sample}}_reheader.log"
        shell:
            """
            # pbsv may not support --sample-name parameter, reheader to set correct sample name
            TMPDIR=$(dirname {output.vcf})/tmp_reheader_$$
            mkdir -p $TMPDIR
            echo "{wildcards.sample}" > $TMPDIR/sample_name.txt
            
            bcftools reheader --samples $TMPDIR/sample_name.txt {input.vcf} > {output.vcf} 2> {log}
            
            rm -rf $TMPDIR
            """

# Helper to construct input VCF list for merging (conditionally include tools based on config)
def _merged_vcf_inputs(wildcards):
    vcfs = []
    if ENABLE_SNIFFLES:
        vcfs.append(f"{TGS_SV_SINGLE_DIR}/{wildcards.sample}/sniffles/{wildcards.sample}.sniffles.vcf")
    if ENABLE_CUTESV:
        vcfs.append(f"{TGS_SV_SINGLE_DIR}/{wildcards.sample}/cutesv/{wildcards.sample}.cutesv.vcf")
    if ENABLE_SVIM:
        vcfs.append(f"{TGS_SV_SINGLE_DIR}/{wildcards.sample}/svim/variants.vcf")
    if PLATFORM == "ONT" and ENABLE_NANOVAR:
        vcfs.append(f"{TGS_SV_SINGLE_DIR}/{wildcards.sample}/nanovar/{wildcards.sample}.nanovar.pass.vcf")
    if PLATFORM == "PACBIO_HIFI" and ENABLE_PBSV:
        vcfs.append(f"{TGS_SV_SINGLE_DIR}/{wildcards.sample}/pbsv/{wildcards.sample}.pbsv.vcf")
    if _SV_CFG.get("hitsv", {}).get("enabled", False):
        vcfs.append(f"{TGS_SV_SINGLE_DIR}/{wildcards.sample}/hitsv/{wildcards.sample}.hitsv.vcf.gz")
    return vcfs

# S09a: Prepare SV VCFs for merging - unify sample names using bcftools
rule prepare_sv_vcfs_for_merge:
    input:
        vcfs = _merged_vcf_inputs
    output:
        merge_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_files.txt",
        caller_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_callers.txt",
        flag = touch(f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/.reheader_done")
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        CONDA_CUTESV_ENV  # Contains bcftools
    container:
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/survivor_prepare/{{sample}}.txt"
    log:
        f"{LOGS_DIR}/survivor_prepare/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {output.merge_list})
        : > {output.merge_list}
        : > {output.caller_list}
        
        # Ensure all input VCFs have the same sample name before merging
        TMPDIR=$(dirname {output.merge_list})/tmp_reheader
        mkdir -p $TMPDIR
        echo "{wildcards.sample}" > $TMPDIR/sample_name.txt
        
        for vcf in {input.vcfs}; do
            # Record the caller name (parent directory of the original VCF, e.g.
            # .../sniffles/HG002.sniffles.vcf -> "sniffles"). This survives reheadering,
            # unlike the file basename (svim's file is just "variants.vcf"), and lets
            # sv_tag_callers stamp INFO/CALLER in the same order as merge_files.txt.
            caller=$(basename $(dirname "$vcf"))
            echo "$caller" >> {output.caller_list}
            # Check if VCF already has the correct sample name
            sample_in_vcf=$(bcftools query -l $vcf 2>/dev/null | head -n1 || echo "UNKNOWN")
            if [ "$sample_in_vcf" != "{wildcards.sample}" ]; then
                # Reheader to unify sample name
                reheadered_vcf=$TMPDIR/$(basename $vcf)
                echo "Reheadering $vcf: $sample_in_vcf -> {wildcards.sample}" >> {log}
                bcftools reheader --samples $TMPDIR/sample_name.txt $vcf > $reheadered_vcf 2>> {log}
                echo "$reheadered_vcf" >> {output.merge_list}
            else
                echo "Sample name already correct in $vcf" >> {log}
                echo "$vcf" >> {output.merge_list}
            fi
        done
        
        echo "Sample name unification completed" >> {log}
        """

# S09b: Merge per-caller SV variants using truvari collapse
# SURVIVOR merge previously distorted variant representations (it keeps one representative
# per cluster, often a poorly-resolved record), yielding a merged call set whose recall was
# LOWER than any single caller. truvari collapse clusters by sequence/size similarity, keeps
# a faithful representative and annotates INFO/NumCollapsed.
# NOTE: The collapse workflow is split into three container-homogeneous rules because the
# biocontainers used with `--sdm apptainer` are single-tool images: the truvari image ships
# ONLY truvari (no bcftools/bgzip/tabix), and the bcftools image ships no truvari. A single
# rule that shells out to both fails in container mode with "bcftools: command not found".
# Splitting keeps each rule inside one image (bcftools OR truvari), matching how the rest of
# the pipeline pairs a conda env with a matching container.

# S09b-1a: Tag each per-caller record with INFO/CALLER (cyvcf2 image: python3, no bcftools).
# truvari collapse only records a folded-record COUNT (NumCollapsed), which over-counts
# support when one caller emits several near-identical records (observed max 22 with only
# 4 callers). Stamping INFO/CALLER before collapse lets sv_annotate_support recover the true
# distinct-caller count via the CollapseId/MatchId link. Kept in its own container-homogeneous
# rule because the bcftools image (used next) ships no python3.
rule sv_tag_callers:
    input:
        merge_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_files.txt",
        caller_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_callers.txt",
        flag = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/.reheader_done"
    output:
        tagged_dir = directory(f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/tagged"),
        tagged_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/tagged_files.txt"
    params:
        script=os.path.join(workflow.basedir, "scripts/sv_pipeline_scripts/tag_caller_source.py")
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        CONDA_CUTESV_ENV
    container:
        config["containers"]["cyvcf2"]
    log:
        f"{LOGS_DIR}/sv_collapse/{{sample}}.tag.log"
    shell:
        """
        set -o pipefail
        : > {log}
        rm -rf {output.tagged_dir}; mkdir -p {output.tagged_dir}
        : > {output.tagged_list}

        # merge_files.txt and merge_callers.txt are written in the same order by
        # prepare_sv_vcfs_for_merge, so paste pairs each VCF with its caller name.
        paste {input.merge_list} {input.caller_list} | while IFS="$(printf '\t')" read -r vcf caller; do
            [ -z "$vcf" ] && continue
            out={output.tagged_dir}/$caller.$(basename "$vcf" .vcf).tagged.vcf
            python3 {params.script:q} \
                --caller "$caller" --input "$vcf" --output "$out" 2>> {log}
            echo "$out" >> {output.tagged_list}
        done
        echo "Tagged per-caller VCFs: $(wc -l < {output.tagged_list})" >> {log}
        """

# S09b-1b: Normalise (QUAL fix), sort and concatenate the tagged VCFs (bcftools image).
rule sv_concat_callers:
    input:
        tagged_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/tagged_files.txt",
        tagged_dir = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/tagged"
    output:
        concat_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.concat.vcf.gz",
        concat_tbi = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.concat.vcf.gz.tbi"
    params:
        filter_opt = "-f PASS,." if COLLAPSE_PASSONLY else ""
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        CONDA_CUTESV_ENV
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/sv_collapse/{{sample}}.concat.log"
    shell:
        """
        set -o pipefail
        WORK=$(dirname {output.concat_vcf})/concat_tmp
        rm -rf $WORK; mkdir -p $WORK
        : > {log}

        # Set missing QUAL to 0 (avoids downstream truvari comparison errors), sort, bgzip
        # and index each tagged VCF. INFO/CALLER was already added by sv_tag_callers.
        : > $WORK/concat_list.txt
        while read -r vcf; do
            [ -z "$vcf" ] && continue
            base=$(basename "$vcf" .vcf)
            awk 'BEGIN{{OFS="\t"}} /^#/{{print;next}} {{if($6==".")$6=0; print}}' "$vcf" \
                | bcftools view {params.filter_opt} -Ou 2>> {log} \
                | bcftools sort -Oz -o $WORK/$base.norm.vcf.gz 2>> {log}
            tabix -f -p vcf $WORK/$base.norm.vcf.gz 2>> {log}
            echo "$WORK/$base.norm.vcf.gz" >> $WORK/concat_list.txt
        done < {input.tagged_list}

        # Concatenate all caller records into a single multi-record, single-sample VCF.
        bcftools concat -a -f $WORK/concat_list.txt -Oz -o {output.concat_vcf} 2>> {log}
        tabix -f -p vcf {output.concat_vcf} 2>> {log}

        echo "Concatenated caller records: $(zcat {output.concat_vcf} | grep -vc '^#' || true)" >> {log}
        rm -rf $WORK
        """

# S09b-2: Collapse redundant multi-caller calls (truvari image only).
rule sv_collapse_callers:
    input:
        concat_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.concat.vcf.gz",
        concat_tbi = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.concat.vcf.gz.tbi",
        # Reference declared as input (like all other rules here) so Snakemake bind-mounts it
        # into the container; truvari reads it to resolve symbolic ALTs. It is an existing
        # source file (no rule produces it), so it is never treated as a build target.
        ref = config["reference"]["genome"]
    output:
        collapsed_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.collapsed.vcf",
        removed_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.collapsed.removed.vcf"
    params:
        pctseq = COLLAPSE_PCTSEQ,
        pctsize = COLLAPSE_PCTSIZE,
        refdist = COLLAPSE_REFDIST,
        sizemin = COLLAPSE_SIZEMIN,
        sizemax = COLLAPSE_SIZEMAX,
        dup_opt = "--dup-to-ins" if COLLAPSE_DUP_TO_INS else "",
        keep = COLLAPSE_KEEP
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        CONDA_TRUVARI_ENV
    container:
        config["containers"].get("truvari", "docker://quay.io/biocontainers/truvari:5.4.0--pyhdfd78af_0")
    benchmark:
        f"{LOGS_DIR}/benchmarks/sv_collapse/{{sample}}.txt"
    log:
        f"{LOGS_DIR}/sv_collapse/{{sample}}.collapse.log"
    shell:
        """
        set -o pipefail
        : > {log}
        # truvari refuses to overwrite existing outputs; clear them for reruns.
        rm -f {output.collapsed_vcf} {output.removed_vcf}

        # Collapse redundant calls. NumCollapsed is the number of additional input records
        # folded into each representative; it does not count distinct callers.
        truvari collapse -i {input.concat_vcf} -o {output.collapsed_vcf} -c {output.removed_vcf} \
            -f {input.ref} -k {params.keep} {params.dup_opt} \
            --pctseq {params.pctseq} --pctsize {params.pctsize} --refdist {params.refdist} \
            --sizemin {params.sizemin} --sizemax {params.sizemax} 2>> {log}
        """

# S09b-2b: Recover true distinct-caller support (truvari image: python3, no bcftools).
# Joins the collapsed representatives (INFO/CollapseId) with the removed/folded records
# (INFO/MatchId) over their per-record INFO/CALLER tags, and writes INFO/CALLERS and
# INFO/NUM_CALLERS. This replaces NumCollapsed (a record count that over-counts support)
# as the basis for the multi-caller filter. Empirically NUM_CALLERS is bounded by the number
# of enabled callers (<=4 on HG002 ONT), whereas NumCollapsed reached 22.
rule sv_annotate_support:
    input:
        collapsed_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.collapsed.vcf",
        removed_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.collapsed.removed.vcf"
    output:
        support_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.collapsed.support.vcf"
    params:
        script=os.path.join(workflow.basedir, "scripts/sv_pipeline_scripts/annotate_caller_support.py"),
        priority_opt = f"--priority {COLLAPSE_PRIORITY}" if COLLAPSE_PRIORITY else ""
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        CONDA_TRUVARI_ENV
    container:
        config["containers"].get("truvari", "docker://quay.io/biocontainers/truvari:5.4.0--pyhdfd78af_0")
    log:
        f"{LOGS_DIR}/sv_collapse/{{sample}}.support.log"
    shell:
        """
        set -o pipefail
        : > {log}
        python3 {params.script:q} \
            --collapsed {input.collapsed_vcf} --removed {input.removed_vcf} \
            {params.priority_opt} --output {output.support_vcf} 2>> {log}
        echo "NUM_CALLERS distribution:" >> {log}
        grep -v "^#" {output.support_vcf} \
            | sed -n 's/.*NUM_CALLERS=\\([0-9]*\\).*/\\1/p' \
            | sort -n | uniq -c >> {log} 2>&1 || true
        """

# S09b-3: Re-sort support-annotated output and keep standard chromosomes only (bcftools).
rule merge_sv_variants:
    input:
        collapsed_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.collapsed.support.vcf"
    output:
        merged_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.vcf"
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        CONDA_CUTESV_ENV
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/sv_collapse/{{sample}}.merge.log"
    shell:
        """
        set -o pipefail
        : > {log}

        # Re-sort (collapse output is unsorted) and keep only standard chromosomes
        # (chr1-22, X, Y, M/MT, with or without the "chr" prefix).
        bcftools sort -Ov {input.collapsed_vcf} 2>> {log} | \
            awk 'BEGIN{{OFS="\t"}} /^#/{{print;next}} $1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/ {{print}}' \
            > {output.merged_vcf}

        BEFORE=$(grep -vc "^#" {input.collapsed_vcf} || true)
        AFTER=$(grep -vc "^#" {output.merged_vcf} || true)
        echo "Collapsed records: $BEFORE" >> {log}
        echo "Representative variants after chromosome filter: $AFTER" >> {log}
        """

# S09b-alt: SURVIVOR merge (tgs_config.sv_detection.merge_method: survivor)
# Produces the same {sample}.merged.vcf as the truvari chain. SURVIVOR runs with support=1;
# INFO/SUPP (number of supporting input VCFs, one per caller) is copied to INFO/NUM_CALLERS so
# filter_sv_variants applies the same min_callers threshold to both methods.
if TGS_SV_MERGE_METHOD == "survivor":
    ruleorder: sv_survivor_merge > merge_sv_variants

    rule sv_survivor_sanitize:
        """Decompress caller VCFs and normalize REF/SVLEN/ALT for SURVIVOR 1.0.7 (bcftools)."""
        input:
            merge_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_files.txt",
            caller_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_callers.txt",
            flag = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/.reheader_done"
        output:
            sanitized_dir = directory(f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/survivor_input"),
            sanitized_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/survivor_files.txt",
            sanitized_callers = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/survivor_callers.txt"
        resources:
            mem_mb = 4096,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["bcftools_module"]
        conda:
            CONDA_CUTESV_ENV
        container:
            config["containers"]["bcftools"]
        log:
            f"{LOGS_DIR}/survivor_merge/{{sample}}.sanitize.log"
        shell:
            """
            set -o pipefail
            : > {log}
            rm -rf {output.sanitized_dir}; mkdir -p {output.sanitized_dir}
            : > {output.sanitized_list}
            : > {output.sanitized_callers}
            i=0
            paste {input.merge_list} {input.caller_list} | \
            while IFS="$(printf '\t')" read -r vcf caller; do
                [ -z "$vcf" ] && continue
                if ! bcftools view -H "$vcf" 2>> {log} | grep . > /dev/null; then
                    echo "Warning: skipping empty caller VCF: $vcf" >> {log}
                    continue
                fi
                i=$((i + 1))
                out={output.sanitized_dir}/$i.$(basename "$vcf" .gz)
                out=${{out%.vcf}}.vcf
                # SURVIVOR 1.0.7 cannot read gzipped input, REF='.', multi-valued SVLEN or
                # <DUP:TANDEM>. Keep the first integer SVLEN (it sizes DEL/INS/DUP).
                bcftools view -Ov "$vcf" 2>> {log} | awk 'BEGIN{{OFS="\t"}}
                    /^#/ {{print; next}}
                    {{
                        if ($4 == ".") $4 = "N"
                        gsub(/<DUP:TANDEM>/, "<DUP>", $5)
                        info = $8
                        if (info != ".") {{
                            n = split(info, kv, ";")
                            out = ""
                            for (k = 1; k <= n; k++) {{
                                if (kv[k] == "") continue
                                if (kv[k] ~ /^SVLEN=/) {{
                                    v = substr(kv[k], 7)
                                    sub(/,.*/, "", v)
                                    if (v !~ /^-?[0-9]+$/) continue
                                    kv[k] = "SVLEN=" v
                                }}
                                out = (out == "" ? kv[k] : out ";" kv[k])
                            }}
                            info = (out == "" ? "." : out)
                        }}
                        $8 = info
                        print
                    }}' > "$out"
                echo "$out" >> {output.sanitized_list}
                echo "$caller" >> {output.sanitized_callers}
            done
            if [ ! -s {output.sanitized_list} ]; then
                echo "Error: no non-empty caller VCFs available for SURVIVOR merge" >> {log}
                exit 1
            fi
            echo "Sanitized caller VCFs: $(wc -l < {output.sanitized_list})" >> {log}
            """

    rule sv_survivor_merge:
        """SURVIVOR merge, then reduce to one sample column and add INFO/CALLERS, NUM_CALLERS."""
        input:
            sanitized_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/survivor_files.txt",
            sanitized_dir = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/survivor_input",
            caller_list = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/survivor_callers.txt"
        output:
            merged_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.vcf"
        params:
            max_distance = SURV_MAX_DISTANCE,
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
            f"{LOGS_DIR}/benchmarks/survivor_merge/{{sample}}.txt"
        log:
            f"{LOGS_DIR}/survivor_merge/{{sample}}.log"
        shell:
            """
            set -o pipefail
            : > {log}
            TMP_MERGED={output.merged_vcf}.survivor.tmp
            # support=1: caller support is enforced later by filter_sv_variants (min_callers).
            SURVIVOR merge {input.sanitized_list} {params.max_distance} 1 1 1 0 {params.min_length} \
                $TMP_MERGED >> {log} 2>&1

            # One sample column (prefer a non-reference GT), standard chromosomes only, and
            # CALLERS/NUM_CALLERS from SUPP_VEC/SUPP. merge_callers.txt has the same order as
            # merge_files.txt and therefore as SUPP_VEC.
            awk -v sample={wildcards.sample} 'BEGIN{{OFS="\t"}}
                FNR == NR {{ caller[++nc] = $1; next }}
                /^##/ {{ print; next }}
                /^#CHROM/ {{
                    print "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\\"Distinct SV callers supporting this variant (from SURVIVOR SUPP_VEC)\\">"
                    print "##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description=\\"Number of distinct SV callers supporting this variant (SURVIVOR SUPP)\\">"
                    print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample
                    next
                }}
                $1 !~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/ {{ next }}
                {{
                    supp = 0; vec = ""
                    n = split($8, kv, ";")
                    for (k = 1; k <= n; k++) {{
                        if (kv[k] ~ /^SUPP=/) supp = substr(kv[k], 6)
                        else if (kv[k] ~ /^SUPP_VEC=/) vec = substr(kv[k], 10)
                    }}
                    callers = ""
                    for (k = 1; k <= length(vec); k++)
                        if (substr(vec, k, 1) == "1") callers = (callers == "" ? caller[k] : callers "," caller[k])
                    if (callers == "") callers = "unknown"
                    $8 = $8 ";CALLERS=" callers ";NUM_CALLERS=" supp
                    col = 0
                    fallback = 0
                    for (k = 10; k <= NF; k++) {{
                        split($k, f, ":")
                        if (f[1] != "./." && f[1] != ".") {{
                            if (fallback == 0) fallback = k
                            if (f[1] != "0/0" && f[1] != "0|0") {{ col = k; break }}
                        }}
                    }}
                    if (col == 0) col = (fallback == 0 ? 10 : fallback)
                    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $col
                }}' {input.caller_list} $TMP_MERGED > {output.merged_vcf}

            echo "SURVIVOR records: $(grep -vc '^#' $TMP_MERGED || true)" >> {log}
            echo "Records after chromosome filter: $(grep -vc '^#' {output.merged_vcf} || true)" >> {log}
            rm -f $TMP_MERGED
            """

if SV_FILTER_ENABLED:
    rule filter_sv_variants:
        input:
            merged_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.vcf"
        output:
            filtered_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.filtered.vcf",
            filtered_vcf_gz = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.filtered.vcf.gz",
            filtered_tbi = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.filtered.vcf.gz.tbi"
        params:
            min_sv_length = max(SURV_MIN_LENGTH, SV_FILTER_MIN_LENGTH),
            # Minimum number of DISTINCT callers required. Both backends write
            # INFO/NUM_CALLERS before this shared filter.
            min_callers = max(1, TGS_MERGE_MIN_CALLERS)
        threads: 2
        resources:
            mem_mb = 4096,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["bcftools_module"]
        conda:
            CONDA_CUTESV_ENV
        container:
            config["containers"]["bcftools"]
        log:
            f"{LOGS_DIR}/survivor_filter/{{sample}}.log"
        shell:
            """
            mkdir -p $(dirname {output.filtered_vcf})

            # Filter the collapsed call set:
            # 1. Require INFO/NUM_CALLERS >= min_callers, i.e. support from at least this many
            #    DISTINCT callers (recovered by sv_annotate_support via CollapseId/MatchId over
            #    per-record INFO/CALLER). This replaces the old NumCollapsed test, which counted
            #    folded records and over-counted support (one caller could satisfy it alone).
            # 2. Drop translocations / breakends (TRA, BND): not size-resolvable and dominated
            #    by false positives in long-read single-sample calling.
            # 3. Drop variants shorter than min_sv_length (allow missing SVLEN=".").
            echo "Filtering collapsed SV VCF..." > {log}
            echo "Min SV length: {params.min_sv_length}; min distinct callers: {params.min_callers}" >> {log}

            FILTER_EXPR='INFO/NUM_CALLERS >= {params.min_callers} && INFO/SVTYPE!="TRA" && INFO/SVTYPE!="BND" && (INFO/SVLEN=="." || abs(INFO/SVLEN) >= {params.min_sv_length})'

            BEFORE=$(grep -v "^#" {input.merged_vcf} | wc -l)
            echo "Variants before filtering: $BEFORE" >> {log}

            bcftools view -i "$FILTER_EXPR" -Ov {input.merged_vcf} | \
                bcftools sort -Ov -o {output.filtered_vcf} 2>> {log}

            AFTER=$(grep -v "^#" {output.filtered_vcf} | wc -l)
            echo "Variants after filtering: $AFTER" >> {log}
            echo "Filtered out $(($BEFORE - $AFTER)) variants" >> {log}

            # Distinct-caller support distribution in filtered VCF
            echo "" >> {log}
            echo "NUM_CALLERS distribution in filtered VCF:" >> {log}
            bcftools query -f '%INFO/NUM_CALLERS\n' {output.filtered_vcf} | sort -n | uniq -c >> {log} 2>&1 || true

            bgzip -f -c {output.filtered_vcf} > {output.filtered_vcf_gz} 2>> {log}
            tabix -f -p vcf {output.filtered_vcf_gz} 2>> {log}
            """
