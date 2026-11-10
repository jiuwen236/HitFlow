# NGS Short-read Structural Variant Detection Workflow (post-joint-calling)

import os

# Per-rule conda environments
ENV_BASIC = "../envs/basic.yaml"
ENV_DOCTOR_MANTA = "../envs/basic.yaml"
ENV_SMOOVE = "../envs/ngs_sv/ngs_sv.smoove.yaml"
ENV_COUNT_LUMPY = "../envs/ngs_sv/ngs_sv.count_lumpy.yaml"
ENV_COUNT_MANTA = "../envs/basic.yaml"
ENV_CNVPYTOR = "../envs/ngs_sv/ngs_sv.cnvpytor.yaml"
ENV_SURVIVOR = "../envs/ngs_sv/ngs_sv.survivor.yaml"

NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
NGS_VARIANTS_SINGLE_DIR = config["output"]["ngs"]["variants"]["single_sample"]
NGS_SV_SINGLE_DIR = config["output"]["ngs"]["sv"]["single_sample"]
NGS_SV_COHORT_DIR = config["output"]["ngs"]["sv"]["cohort"]
LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]
NGS_VARIANTS_COHORT_DIR = config["output"]["ngs"]["variants"]["cohort"]
NGS_PREPROCESSING_DIR = config["output"]["ngs"]["preprocessing"]

SAMPLES_VCF_CSV = config["samples_vcf_csv"]
SCAN_DIRS = config.get("scan_ngs_dirs", [])

# Discover samples reused from joint calling
import csv, glob

def _load_vcf_samples(csv_path):
    samples = []
    if os.path.isfile(csv_path):
        with open(csv_path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                samples.append(row["sample"])
    return set(samples)

# Union of samples known in final BAM directories and VCF CSV

def get_ngs_samples():
    """Discover NGS SV samples following the same source as bwa_mem_pipe.
    Prefer the main workflow's FASTQ-derived sample list, without scanning final/.
    """
    # Prefer using the main Snakefile's ALL_SAMPLES via helper to match bwa_mem_pipe
    try:
        if 'get_all_sample_ids' in globals():
            return sorted(get_all_sample_ids())
    except NameError:
        pass

    # Fallback: accept samples explicitly listed in samples_vcf_csv (if provided)
    csv_samples = _load_vcf_samples(SAMPLES_VCF_CSV)
    return sorted(csv_samples)

REFERENCE_FA = config["reference"]["genome"]

def _ngs_sv_caller_inputs(wildcards):
    """Return caller inputs, gating VerifyBamID2 on its workflow switch."""
    inputs = {
        "bam": f"{NGS_ALIGNMENT_DIR}/{wildcards.sample}/bam/{wildcards.sample}.dup.sorted.recal.bam",
        "ref": REFERENCE_FA,
    }
    if config.get("verifybamid2", {}).get("enabled", True):
        inputs["verify"] = f"{NGS_ALIGNMENT_DIR}/{wildcards.sample}/qc/{wildcards.sample}.verifybamid.report"
    return inputs

# SURVIVOR merge parameters for NGS SV
_SURVIVOR_NGS_CFG = config.get("ngs_sv", {}).get("survivor", {})
SURV_MAX_DISTANCE = _SURVIVOR_NGS_CFG.get("max_distance", 1000)
SURV_MIN_SUPPORT = _SURVIVOR_NGS_CFG.get("min_support", 2)
SURV_MIN_LENGTH = _SURVIVOR_NGS_CFG.get("min_length", 30)
SURV_POP_SUPPORT = _SURVIVOR_NGS_CFG.get("population_support", 2)
SURV_ENABLE_POP_MERGE = _SURVIVOR_NGS_CFG.get("enable_population_merge", True)

# Per-sample multi-caller merge method: "survivor" (default) or "truvari" (truvari collapse).
# Both produce {sample}.ngs_sv.merged.vcf.gz; the cohort-level merge always uses SURVIVOR.
NGS_SV_MERGE_METHOD = str(config.get("ngs_sv", {}).get("merge_method", "survivor")).lower()
if NGS_SV_MERGE_METHOD not in ("survivor", "truvari"):
    raise ValueError(f"ngs_sv.merge_method must be 'survivor' or 'truvari', got '{NGS_SV_MERGE_METHOD}'")

# truvari collapse parameters for NGS SV. Short-read calls are mostly symbolic
# (no inserted sequence), so sequence similarity is ignored by default.
_COLLAPSE_NGS_CFG = config.get("ngs_sv", {}).get("collapse", {})
NGS_COLLAPSE_PCTSEQ = _COLLAPSE_NGS_CFG.get("pctseq", 0)
NGS_COLLAPSE_PCTSIZE = _COLLAPSE_NGS_CFG.get("pctsize", 0.7)
NGS_COLLAPSE_REFDIST = _COLLAPSE_NGS_CFG.get("refdist", 500)
# truvari passes variants outside sizemin/sizemax through uncollapsed (NUM_CALLERS=1), so
# the default disables the upper bound to keep large read-depth CNVs comparable.
NGS_COLLAPSE_SIZEMAX = _COLLAPSE_NGS_CFG.get("sizemax", -1)
NGS_COLLAPSE_KEEP = _COLLAPSE_NGS_CFG.get("keep", "common")
NGS_COLLAPSE_MIN_LENGTH = _COLLAPSE_NGS_CFG.get("min_length", SURV_MIN_LENGTH)
# Shared support threshold for either backend. Legacy backend-specific keys remain fallbacks.
_MERGE_NGS_CFG = config.get("ngs_sv", {}).get("merge", {})
NGS_MERGE_MIN_CALLERS = _MERGE_NGS_CFG.get(
    "min_callers", _COLLAPSE_NGS_CFG.get("min_callers", SURV_MIN_SUPPORT)
)

# Tool toggles (unused)
_MANTA_ENABLED = config.get("ngs_sv", {}).get("manta", {}).get("enabled", True)
_SMOOVE_ENABLED = config.get("ngs_sv", {}).get("smoove", {}).get("enabled", True)
_CNVPYTOR_ENABLED = config.get("ngs_sv", {}).get("cnvpytor", {}).get("enabled", False)
_CNVPYTOR_TO_VCF = str(config.get("ngs_sv", {}).get("cnvpytor", {}).get("to_vcf_script", "")).strip()
_MANTA_CALL_REGIONS = str(config.get("ngs_sv", {}).get("manta", {}).get("call_regions", "")).strip()

def _manta_inputs(wildcards):
    inputs = _ngs_sv_caller_inputs(wildcards)
    if _MANTA_CALL_REGIONS:
        inputs["call_regions"] = _MANTA_CALL_REGIONS
        inputs["call_regions_tbi"] = f"{_MANTA_CALL_REGIONS}.tbi"
    return inputs

# S01: Manta germline workflow per sample
rule manta_germline:
    input:
        unpack(_manta_inputs)
    output:
        vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/manta/{{sample}}.manta.vcf.gz",
        tbi = f"{NGS_SV_SINGLE_DIR}/{{sample}}/manta/{{sample}}.manta.vcf.gz.tbi"
    params:
        rundir = f"{NGS_SV_SINGLE_DIR}/{{sample}}/manta",
        call_regions = _MANTA_CALL_REGIONS,
        import_manta = workflow_path("scripts/pipeline/import_manta.sh")
    threads: 10
    resources:
        mem_mb = 5120,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["manta_module"],
        config["tool_modules"].get("python_module", "pipeline/python/2.7.12")
    conda:
        ENV_BASIC
    container:
        config["containers"]["manta"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/manta/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/ngs_sv/manta/{{sample}}.log"
    shell:
        """
        # configManta.py in bin/manta/
        source {params.import_manta}
        
        JOB_TMP_DIR="{TMP_DIR}/manta_{wildcards.sample}_$$"
        mkdir -p ${{JOB_TMP_DIR}} {params.rundir} $(dirname {output.vcf})
        CALL_REGIONS_OPT=""
        if [ -n "{params.call_regions}" ]; then
            CALL_REGIONS_OPT="--callRegions {params.call_regions}"
        fi
        rm -rf {params.rundir}/*
        configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir {params.rundir} $CALL_REGIONS_OPT >> {log} 2>&1
        {params.rundir}/runWorkflow.py -m local -j {threads} >> {log} 2>&1
        mv {params.rundir}/results/variants/diploidSV.vcf.gz {output.vcf}
        mv {params.rundir}/results/variants/diploidSV.vcf.gz.tbi {output.tbi}
        """

# S02: doctor_manta (vcf normalization from upstream pipeline script)
DOCTOR_MANTA_CANDIDATES = [
    os.path.join(config["ngs_scripts"], "doctor_manta.1.py"),
    os.path.join(config["sv_pipeline_scripts"], "doctor_manta.1.py"),
]
DOCTOR_MANTA = next((p for p in DOCTOR_MANTA_CANDIDATES if os.path.exists(p)), DOCTOR_MANTA_CANDIDATES[0])

rule doctor_manta:
    input:
        vcf = rules.manta_germline.output.vcf
    output:
        vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/manta/{{sample}}.manta.doctored.vcf.gz",
        tbi = f"{NGS_SV_SINGLE_DIR}/{{sample}}/manta/{{sample}}.manta.doctored.vcf.gz.tbi"
    params:
        sv_pipeline_scripts = config["sv_pipeline_scripts"]
    threads: 1
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["manta_module"],
        config["tool_modules"]["tabix_module"]
    conda:
        ENV_DOCTOR_MANTA
    container:
        config["containers"]["cyvcf2"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/doctor_manta/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/ngs_sv/doctor_manta/{{sample}}.log"
    shell:
        """
        set -euo pipefail
        export SV_PIPELINE_SCRIPTS="{params.sv_pipeline_scripts}"
        zcat {input.vcf} | python {DOCTOR_MANTA} -m 700 | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# S03: Smoove call (lumpy wrapper)
rule smoove_call:
    input:
        unpack(_ngs_sv_caller_inputs)
    output:
        vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/smoove/{{sample}}-smoove.genotyped.vcf.gz"
    params:
        exclude_regions = config.get("ngs_sv", {}).get("smoove", {}).get("exclude_regions", "")
    threads: 1
    resources:
        mem_mb = 10240,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["smoove_module"],
        config["tool_modules"]["lumpy_module"],
        config["tool_modules"]["gsort_module"],
        config["tool_modules"]["bcftools_module"],
        config["tool_modules"]["svtyper_module"]
    conda:
        ENV_SMOOVE
    container:
        config["containers"]["smoove"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/smoove/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/ngs_sv/smoove/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        EXCLUDE_OPT=""
        if [ -n "{params.exclude_regions}" ] && [ -f "{params.exclude_regions}" ]; then
            EXCLUDE_OPT="--exclude {params.exclude_regions}"
        fi
        smoove call \
          -outdir $(dirname {output.vcf}) \
          --name {wildcards.sample} \
          --fasta {input.ref} \
          $EXCLUDE_OPT \
          --noextrafilters --genotype {input.bam} >> {log} 2>&1
        """

# S04: Count lumpy/smoove SVs
COUNT_LUMPY_SH = os.path.join(config["sv_pipeline_scripts"], "count_lumpy.sh")
COUNT_MANTA_SH = os.path.join(config["sv_pipeline_scripts"], "count_manta.sh")

rule count_lumpy:
    input:
        vcf = rules.smoove_call.output.vcf
    output:
        txt = f"{NGS_SV_SINGLE_DIR}/{{sample}}/lumpy/count_lumpy/{{sample}}.lumpy.counts.1.txt"
    params:
        sv_pipeline_scripts = config["sv_pipeline_scripts"]
    threads: 1
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_COUNT_LUMPY
    container:
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/count_lumpy/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/ngs_sv/count_lumpy/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {output.txt})
        export SV_PIPELINE_SCRIPTS="{params.sv_pipeline_scripts}"
        bash {COUNT_LUMPY_SH} {input.vcf} {wildcards.sample} $(dirname {output.txt}) >> {log} 2>&1
        """

rule count_manta:
    input:
        vcf = rules.manta_germline.output.vcf
    output:
        txt = f"{NGS_SV_SINGLE_DIR}/{{sample}}/manta/{{sample}}.manta.counts.1.txt"
    params:
        sv_pipeline_scripts = config["sv_pipeline_scripts"]
    threads: 1
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_COUNT_MANTA
    container:
        config["containers"]["bcftools"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/count_manta/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/ngs_sv/count_manta/{{sample}}.log"
    shell:
        """
        export SV_PIPELINE_SCRIPTS="{params.sv_pipeline_scripts}"
        mkdir -p $(dirname {output.txt})
        bash {COUNT_MANTA_SH} {input.vcf} {wildcards.sample} $(dirname {output.txt}) >> {log} 2>&1
        """

# S05: CNVpytor (optional, controlled by config)
if config.get("ngs_sv", {}).get("cnvpytor", {}).get("enabled", False):

    rule cnvpytor_download_resources:
        """Download CNVpytor reference genome resource files (GC content and mask files for hg19/hg38)"""
        output:
            flag = f"{LOGS_DIR}/.cnvpytor_resources_downloaded"
        threads: 1
        resources:
            mem_mb = 2048,
            partition = config["cluster"]["partition"],
            time = "01:00:00"
        envmodules:
            config["tool_modules"]["cnvpytor_module"]
        conda:
            ENV_CNVPYTOR
        container:
            config["containers"]["cnvpytor"]
        log:
            f"{LOGS_DIR}/ngs_sv/cnvpytor_download.log"
        shell:
            """
            # Download reference genome resource files if missing
            cnvpytor -download >> {log} 2>&1
            touch {output.flag}
            """

    rule cnvpytor:
        input:
            bam = f"{NGS_ALIGNMENT_DIR}/{{sample}}/bam/{{sample}}.dup.sorted.recal.bam",
            resources_flag = rules.cnvpytor_download_resources.output.flag
        output:
            summary = f"{NGS_SV_SINGLE_DIR}/{{sample}}/cnvpytor/{{sample}}.summary.txt"
        params:
            outdir = f"{NGS_SV_SINGLE_DIR}/{{sample}}/cnvpytor",
            bin_size = config.get("ngs_sv", {}).get("cnvpytor", {}).get("bin_size", 100),
            # Only process main chromosomes to avoid performance issues with HLA/decoy sequences
            chromosomes = config.get("ngs_sv", {}).get("cnvpytor", {}).get("chromosomes", 
                "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY")
        threads: 1
        resources:
            mem_mb = 30 * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["cnvpytor_module"]
        conda:
            ENV_CNVPYTOR
        container:
            config["containers"]["cnvpytor"]
        benchmark:
            f"{LOGS_DIR}/benchmarks/ngs_sv/cnvpytor/{{sample}}.tsv"
        log:
            f"{LOGS_DIR}/ngs_sv/cnvpytor/{{sample}}.log"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir}
            # Build CNVpytor root, histogram, partition, and call
            # Use -chrom to restrict to main chromosomes, avoiding HLA and decoy sequences
            cnvpytor -root {params.outdir}/out.pytor -rd {input.bam} -chrom {params.chromosomes} >> {log} 2>&1
            cnvpytor -root {params.outdir}/out.pytor -his {params.bin_size} >> {log} 2>&1
            cnvpytor -root {params.outdir}/out.pytor -partition {params.bin_size} >> {log} 2>&1
            cnvpytor -root {params.outdir}/out.pytor -call {params.bin_size} > {output.summary} 2>> {log}
            """

    # Optional: convert CNVpytor outputs to VCF if converter provided, include in merge
    if str(config.get("ngs_sv", {}).get("cnvpytor", {}).get("to_vcf_script", "")).strip() != "":
        rule cnvpytor_to_vcf:
            input:
                summary = rules.cnvpytor.output.summary
            output:
                vcfgz = f"{NGS_SV_SINGLE_DIR}/{{sample}}/cnvpytor/{{sample}}.cnvpytor.vcf.gz",
                tbi = f"{NGS_SV_SINGLE_DIR}/{{sample}}/cnvpytor/{{sample}}.cnvpytor.vcf.gz.tbi"
            params:
                script = str(config.get("ngs_sv", {}).get("cnvpytor", {}).get("to_vcf_script")),
                sample = "{sample}",
                ref_base = str(config.get("ngs_sv", {}).get("cnvpytor", {}).get("ref_base", ".")),
                vcf_reference_header = str(config.get("ngs_sv", {}).get("cnvpytor", {}).get("vcf_reference_header", ""))
            threads: 1
            resources:
                mem_mb = 1024,
                partition = config["cluster"]["partition"],
                time = config["cluster"]["time"]
            envmodules:
                config["tool_modules"]["bcftools_module"],
                config["tool_modules"]["tabix_module"]
            conda:
                ENV_CNVPYTOR
            container:
                config["containers"]["cyvcf2"]
            benchmark:
                f"{LOGS_DIR}/benchmarks/ngs_sv/cnvpytor_to_vcf/{{sample}}.tsv"
            log:
                f"{LOGS_DIR}/ngs_sv/cnvpytor_to_vcf/{{sample}}.log"
            shell:
                """
                mkdir -p $(dirname {output.vcfgz})
                # Expect converter script to read CNVpytor summary and emit VCF to stdout
                {params.script} --sample {params.sample} --input {input.summary} \
                  --ref-base {params.ref_base} \
                  --reference "{params.vcf_reference_header}" \
                  | bgzip -c > {output.vcfgz}
                tabix -f -p vcf {output.vcfgz}
                """

# S06: SURVIVOR merges for NGS SVs
# Split into 3 rules to use correct containers:
# 1. ngs_sv_survivor_sanitize_vcfs: preprocess VCFs using bcftools container (has bgzip/tabix)
# 2. ngs_sv_survivor_merge_run: run SURVIVOR merge using survivor container
# 3. ngs_sv_survivor_merge_sample: convert to one sample, compress and index (bcftools)

def _ngs_sv_per_sample_vcfs(wildcards):
    vcfs = []
    if config.get("ngs_sv", {}).get("manta", {}).get("enabled", True):
        vcfs.append(f"{NGS_SV_SINGLE_DIR}/{wildcards.sample}/manta/{wildcards.sample}.manta.doctored.vcf.gz")
    if config.get("ngs_sv", {}).get("smoove", {}).get("enabled", True):
        vcfs.append(f"{NGS_SV_SINGLE_DIR}/{wildcards.sample}/smoove/{wildcards.sample}-smoove.genotyped.vcf.gz")
    if config.get("ngs_sv", {}).get("hitsv", {}).get("enabled", False):
        vcfs.append(f"{NGS_SV_SINGLE_DIR}/{wildcards.sample}/hitsv/{wildcards.sample}.hitsv.vcf.gz")
    # Optionally include CNVpytor-derived VCF if converter provided
    cnv_cfg = config.get("ngs_sv", {}).get("cnvpytor", {})
    to_vcf = str(cnv_cfg.get("to_vcf_script", "")).strip()
    if cnv_cfg.get("enabled", False) and to_vcf != "":
        vcfs.append(f"{NGS_SV_SINGLE_DIR}/{wildcards.sample}/cnvpytor/{wildcards.sample}.cnvpytor.vcf.gz")
    return vcfs

rule ngs_sv_survivor_sanitize_vcfs:
    """Preprocess VCFs for SURVIVOR: fix REF/SVLEN and normalize ALT types.
    Uses bcftools container which includes bgzip/tabix."""
    input:
        vcfs = _ngs_sv_per_sample_vcfs
    output:
        merge_list = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_files.txt",
        caller_list = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_callers.txt",
        done_flag = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/.sanitize_done"
    threads: 1
    resources:
        mem_mb = 4 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["tabix_module"],
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BASIC
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/ngs_sv/survivor_merge/{{sample}}.sanitize.log"
    shell:
        """
        mkdir -p $(dirname {output.merge_list})
        TMPDIR=$(dirname {output.merge_list})/tmp_sanitize_$$
        mkdir -p $TMPDIR
        : > {output.merge_list}
        : > {output.caller_list}
        
        # Filter out empty VCF files and fix incompatible patterns for SURVIVOR
        # 1) Replace REF '.' with 'N' (SURVIVOR can't handle '.')
        # 2) Reduce INFO/SVLEN to its first integer value (multi-valued/invalid SVLEN can crash
        #    old SURVIVOR). Do not drop it: SURVIVOR derives DEL/INS/DUP length from SVLEN and
        #    would otherwise discard all of them via min_length.
        # 3) Normalize ALT by mapping '<DUP:TANDEM>' to '<DUP>' to avoid parser bugs
        for vcf in {input.vcfs}; do
            # Count non-header lines (actual variant records)
            variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l)
            if [ "$variant_count" -gt 0 ]; then
                # Check if VCF has '.' in REF column (SURVIVOR incompatible)
                ref_dot_count=$(zcat "$vcf" | grep -v "^#" | awk '{{if ($4 == ".") print}}' | wc -l)
                
                # choose source vcf for next sanitization stage
                chosen_vcf="$vcf"
                if [ "$ref_dot_count" -gt 0 ]; then
                    # Fix VCF: replace '.' REF with 'N' for SURVIVOR compatibility
                    fixed_vcf="$TMPDIR/$(basename "$vcf" .gz).fixed.vcf.gz"
                    echo "Info: Fixing REF column in $vcf ($ref_dot_count variants with '.' REF)" >> {log}
                    zcat "$vcf" | awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{if ($4 == ".") $4 = "N"; print}}' | bgzip -c > "$fixed_vcf"
                    tabix -f -p vcf "$fixed_vcf"
                    chosen_vcf="$fixed_vcf"
                fi

                # Additional sanitization for SURVIVOR stability
                # Output plain VCF (not .vcf.gz) — SURVIVOR 1.0.7 cannot read gzipped input
                sanitized_vcf="$(dirname {output.merge_list})/$(basename "$chosen_vcf" .vcf.gz | sed 's/.fixed$//')_sanitized.vcf"
                echo "Info: Sanitizing for SURVIVOR (single-valued INFO/SVLEN; ALT <DUP:TANDEM>-><DUP>): $chosen_vcf" >> {log}
                zcat "$chosen_vcf" | awk 'BEGIN{{OFS="\t"}}
                    /^#/ {{print; next}}
                    {{
                        # Normalize ALT
                        gsub(/<DUP:TANDEM>/, "<DUP>", $5)
                        # Keep the first integer SVLEN value; drop non-integer SVLEN
                        info=$8
                        if (info != ".") {{
                            n = split(info, kv, ";")
                            out = ""
                            for (i = 1; i <= n; i++) {{
                                if (kv[i] == "") continue
                                if (kv[i] ~ /^SVLEN=/) {{
                                    v = substr(kv[i], 7)
                                    sub(/,.*/, "", v)
                                    if (v !~ /^-?[0-9]+$/) continue
                                    kv[i] = "SVLEN=" v
                                }}
                                out = (out == "" ? kv[i] : out ";" kv[i])
                            }}
                            info = (out == "" ? "." : out)
                        }}
                        $8 = info
                        print
                    }}' > "$sanitized_vcf"
                echo "$sanitized_vcf" >> {output.merge_list}
                basename "$(dirname "$vcf")" >> {output.caller_list}
            else
                echo "Warning: Skipping empty VCF file: $vcf (no variant records)" >> {log}
            fi
        done
        
        # Clean up temp fixed files
        rm -rf $TMPDIR
        
        # Check if we have at least one VCF with variants
        vcf_count=$(wc -l < {output.merge_list})
        if [ "$vcf_count" -eq 0 ]; then
            echo "Error: No VCF files with variants found for merging" >> {log} 2>&1
            exit 1
        fi
        
        touch {output.done_flag}
        """

rule ngs_sv_survivor_merge_run:
    """Run SURVIVOR merge on sanitized VCFs. Uses survivor container."""
    input:
        merge_list = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_files.txt",
        done_flag = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/.sanitize_done"
    output:
        merged_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.ngs_sv.merged.vcf"
    params:
        max_distance = SURV_MAX_DISTANCE,
        min_support = NGS_MERGE_MIN_CALLERS,
        min_length = SURV_MIN_LENGTH
    threads: 1
    resources:
        mem_mb = 20 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["survivor_module"]
    conda:
        ENV_SURVIVOR
    container:
        config["containers"]["survivor"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/survivor_merge/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/ngs_sv/survivor_merge/{{sample}}.log"
    shell:
        """
        # Log binary resolution to help diagnose potential PATH/library issues
        echo "SURVIVOR which: $(command -v SURVIVOR)" >> {log} 2>&1 || true
        SURVIVOR --version >> {log} 2>&1 || true
        
        SURVIVOR merge \
            {input.merge_list} \
            {params.max_distance} {params.min_support} 1 1 0 {params.min_length} \
            {output.merged_vcf} >> {log} 2>&1
        """

rule ngs_sv_survivor_merge_sample:
    """Convert SURVIVOR output to a clean single-sample VCF, compress and index it."""
    input:
        merged_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.ngs_sv.merged.vcf",
        caller_list = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/merge_callers.txt"
    output:
        merged_vcf_gz = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.ngs_sv.merged.vcf.gz",
        tbi = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.ngs_sv.merged.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 4 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["tabix_module"],
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BASIC
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/ngs_sv/survivor_merge/{{sample}}.compress.log"
    shell:
        """
        # Validate and sanitize SURVIVOR output (known to produce corrupted records)
        TMP_VCF_CLEAN=$(dirname {output.merged_vcf_gz})/{wildcards.sample}.ngs_sv.merged.clean.vcf
        SORT_TMP=$(dirname {output.merged_vcf_gz})/tmp_sort_$$
        mkdir -p "$SORT_TMP"
        echo "Info: Validating and cleaning SURVIVOR output" >> {log}
        # Keep valid records, expose distinct-caller support, and reduce SURVIVOR's
        # per-caller sample columns to one sample column for a stable downstream contract.
        awk -v sample={wildcards.sample} 'BEGIN{{OFS="\t"}}
            FNR == NR {{ caller[++nc] = $1; next }}
            /^##/ {{print; next}}
            /^#CHROM/ {{
                print "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\\"Distinct SV callers supporting this variant (from SURVIVOR SUPP_VEC)\\">"
                print "##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description=\\"Number of distinct SV callers supporting this variant (SURVIVOR SUPP)\\">"
                print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample
                next
            }}
            {{
                if (NF >= 10 && $1 ~ /^(chr)?[0-9XYM]+$/ && $2 ~ /^[0-9]+$/ && length($4) > 0 && length($5) > 0 && $4 != "." && $5 != ".") {{
                    n = split($8, fields, ";")
                    cleaned = ""
                    supp = 0
                    vec = ""
                    is_breakend = (index($5, "[") > 0 || index($5, "]") > 0)
                    for (i = 1; i <= n; i++) {{
                        if (fields[i] ~ /^SUPP=/) supp = substr(fields[i], 6)
                        else if (fields[i] ~ /^SUPP_VEC=/) vec = substr(fields[i], 10)
                        if (fields[i] ~ /^END=/) {{
                            split(fields[i], end_value, "=")
                            if ((end_value[2] + 0) < ($2 + 0)) {{
                                if (!is_breakend) {{
                                    invalid_end++
                                }} else {{
                                    removed_end++
                                }}
                                continue
                            }}
                        }}
                        cleaned = (cleaned == "" ? fields[i] : cleaned ";" fields[i])
                    }}
                    callers = ""
                    for (i = 1; i <= length(vec); i++)
                        if (substr(vec, i, 1) == "1")
                            callers = (callers == "" ? caller[i] : callers "," caller[i])
                    if (callers == "") callers = "unknown"
                    $8 = (cleaned == "" ? "." : cleaned) ";CALLERS=" callers ";NUM_CALLERS=" supp
                    col = 0
                    fallback = 0
                    for (i = 10; i <= NF; i++) {{
                        split($i, gt_fields, ":")
                        if (gt_fields[1] != "./." && gt_fields[1] != ".") {{
                            if (fallback == 0) fallback = i
                            if (gt_fields[1] != "0/0" && gt_fields[1] != "0|0") {{
                                col = i
                                break
                            }}
                        }}
                    }}
                    if (col == 0) col = (fallback == 0 ? 10 : fallback)
                    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $col
                }}
            }}
            END {{
                if (removed_end > 0) print "Warning: Removed END<POS from " removed_end " breakend records; ALT/CHR2 retain the remote breakpoint." > "/dev/stderr"
                if (invalid_end > 0) {{
                    print "Error: Found " invalid_end " non-breakend records with END<POS." > "/dev/stderr"
                    exit 2
                }}
            }}' "{input.caller_list}" "{input.merged_vcf}" > "$TMP_VCF_CLEAN" 2>> {log}
        
        # Check if we have valid variants after cleaning (use || true to avoid strict mode failure)
        VALID_VARIANTS=$(grep -v "^#" "$TMP_VCF_CLEAN" 2>/dev/null | wc -l || echo "0")
        echo "Info: Found $VALID_VARIANTS valid variants after sanitization" >> {log}
        
        if [ "$VALID_VARIANTS" -eq 0 ]; then
            echo "Warning: No valid variants after SURVIVOR merge and sanitization. Creating minimal VCF." >> {log}
            # Preserve metadata but emit a single-sample column header.
            grep "^##" "{input.merged_vcf}" > "$TMP_VCF_CLEAN" || true
            printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' \
                "{wildcards.sample}" >> "$TMP_VCF_CLEAN"
        fi
        
        bcftools sort -T "$SORT_TMP" -Oz -o {output.merged_vcf_gz} "$TMP_VCF_CLEAN"
        tabix -f -p vcf {output.merged_vcf_gz}
        
        echo "Info: SURVIVOR output SVTYPE counts:" >> {log}
        bcftools query -f '%INFO/SVTYPE\\n' {output.merged_vcf_gz} | sort | uniq -c >> {log}
        DEL_INS_VARIANTS=$(bcftools view -H -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="INS"' {output.merged_vcf_gz} | wc -l)
        if [ "$DEL_INS_VARIANTS" -eq 0 ]; then
            echo "Warning: SURVIVOR output contains no DEL or INS records; inspect caller compatibility before downstream DEL/INS analysis." >> {log}
        fi

        # Clean up temporary files
        rm -f "$TMP_VCF_CLEAN"
        rm -rf "$SORT_TMP"
        
        # Clean up sanitized VCFs from previous step (they are intermediate files)
        rm -f $(dirname {output.merged_vcf_gz})/*_sanitized.vcf 2>/dev/null || true
        rm -f $(dirname {output.merged_vcf_gz})/*_sanitized.vcf.gz* 2>/dev/null || true
        """

# S06-alt: truvari collapse merge for NGS SVs (ngs_sv.merge_method: truvari)
# Mirrors the TGS collapse chain (rules/tgs_sv.smk): tag INFO/CALLER -> concat -> collapse ->
# distinct-caller support (NUM_CALLERS) -> filter. Unlike SURVIVOR it keeps one original
# record per cluster instead of rebuilding it. Split by container: cyvcf2 (python3),
# bcftools, truvari. Produces the same {sample}.ngs_sv.merged.vcf.gz as the SURVIVOR chain.
if NGS_SV_MERGE_METHOD == "truvari":
    ruleorder: ngs_sv_collapse_merge_sample > ngs_sv_survivor_merge_sample

    rule ngs_sv_collapse_tag_callers:
        """Stamp INFO/CALLER on each caller VCF (caller = parent directory name)."""
        input:
            vcfs = _ngs_sv_per_sample_vcfs
        output:
            tagged_dir = directory(f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/tagged"),
            tagged_list = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/tagged_files.txt"
        params:
            script = os.path.join(workflow.basedir, "scripts/sv_pipeline_scripts/tag_caller_source.py")
        threads: 1
        resources:
            mem_mb = 4 * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        conda:
            ENV_BASIC
        container:
            config["containers"]["cyvcf2"]
        log:
            f"{LOGS_DIR}/ngs_sv/collapse_merge/{{sample}}.tag.log"
        shell:
            """
            set -o pipefail
            : > {log}
            rm -rf {output.tagged_dir}; mkdir -p {output.tagged_dir}
            : > {output.tagged_list}
            for vcf in {input.vcfs}; do
                caller=$(basename $(dirname "$vcf"))
                out={output.tagged_dir}/$caller.vcf
                python3 {params.script:q} --caller "$caller" --input "$vcf" --output "$out" 2>> {log}
                echo "$out" >> {output.tagged_list}
            done
            echo "Tagged caller VCFs: $(wc -l < {output.tagged_list})" >> {log}
            """

    rule ngs_sv_collapse_concat:
        """Unify sample name/contigs, set missing QUAL to 0, sort and concatenate (bcftools)."""
        input:
            tagged_list = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/tagged_files.txt",
            tagged_dir = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/tagged",
            fai = f"{REFERENCE_FA}.fai"
        output:
            concat_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/{{sample}}.concat.vcf.gz",
            concat_tbi = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/{{sample}}.concat.vcf.gz.tbi"
        threads: 1
        resources:
            mem_mb = 8 * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["tabix_module"],
            config["tool_modules"]["bcftools_module"]
        conda:
            ENV_BASIC
        container:
            config["containers"]["bcftools"]
        log:
            f"{LOGS_DIR}/ngs_sv/collapse_merge/{{sample}}.concat.log"
        shell:
            """
            set -o pipefail
            : > {log}
            WORK=$(dirname {output.concat_vcf})/concat_tmp
            rm -rf $WORK; mkdir -p $WORK
            echo "{wildcards.sample}" > $WORK/sample_name.txt
            : > $WORK/concat_list.txt
            while read -r vcf; do
                [ -z "$vcf" ] && continue
                base=$(basename "$vcf" .vcf)
                # Sites-only VCFs (e.g. CNVpytor) get a GT=./. column so all inputs share one sample.
                # Missing QUAL is set to 0 for truvari, as in the TGS chain.
                awk 'BEGIN{{OFS="\t"}}
                    /^##FORMAT=<ID=GT,/ {{has_gt=1}}
                    /^##/ {{print; next}}
                    /^#CHROM/ {{
                        sites_only = (NF == 8)
                        if (sites_only) {{
                            if (!has_gt) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"
                            print $0, "FORMAT", "SAMPLE"
                        }} else print
                        next
                    }}
                    {{
                        if ($6 == ".") $6 = 0
                        if (sites_only) print $0, "GT", "./."
                        else print
                    }}' "$vcf" > $WORK/$base.fmt.vcf
                if [ "$(bcftools query -l $WORK/$base.fmt.vcf | wc -l)" -ne 1 ]; then
                    echo "Error: $vcf must contain exactly one sample" >> {log}
                    exit 1
                fi
                # Same sample name and contig dictionary for every caller, required by bcftools concat.
                bcftools reheader --samples $WORK/sample_name.txt --fai {input.fai} $WORK/$base.fmt.vcf 2>> {log} \
                    | bcftools sort -T $WORK/sort_$base -Oz -o $WORK/$base.norm.vcf.gz 2>> {log}
                tabix -f -p vcf $WORK/$base.norm.vcf.gz 2>> {log}
                echo "$WORK/$base.norm.vcf.gz" >> $WORK/concat_list.txt
            done < {input.tagged_list}

            bcftools concat -a -f $WORK/concat_list.txt -Oz -o {output.concat_vcf} 2>> {log}
            tabix -f -p vcf {output.concat_vcf} 2>> {log}
            echo "Concatenated caller records: $(zcat {output.concat_vcf} | grep -vc '^#' || true)" >> {log}
            rm -rf $WORK
            """

    rule ngs_sv_collapse_run:
        """Collapse redundant multi-caller calls and recover distinct-caller support (truvari)."""
        input:
            concat_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/{{sample}}.concat.vcf.gz",
            concat_tbi = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/{{sample}}.concat.vcf.gz.tbi",
            ref = REFERENCE_FA
        output:
            collapsed_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/{{sample}}.collapsed.vcf",
            removed_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/{{sample}}.collapsed.removed.vcf",
            support_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/{{sample}}.collapsed.support.vcf"
        params:
            pctseq = NGS_COLLAPSE_PCTSEQ,
            pctsize = NGS_COLLAPSE_PCTSIZE,
            refdist = NGS_COLLAPSE_REFDIST,
            sizemin = NGS_COLLAPSE_MIN_LENGTH,
            sizemax = NGS_COLLAPSE_SIZEMAX,
            keep = NGS_COLLAPSE_KEEP,
            script = os.path.join(workflow.basedir, "scripts/sv_pipeline_scripts/annotate_caller_support.py")
        threads: 1
        resources:
            mem_mb = 8 * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        conda:
            "../envs/long_read/ont_truvari.yaml"
        container:
            config["containers"].get("truvari", "docker://quay.io/biocontainers/truvari:5.4.0--pyhdfd78af_0")
        benchmark:
            f"{LOGS_DIR}/benchmarks/ngs_sv/collapse_merge/{{sample}}.tsv"
        log:
            f"{LOGS_DIR}/ngs_sv/collapse_merge/{{sample}}.collapse.log"
        shell:
            """
            set -o pipefail
            : > {log}
            # truvari refuses to overwrite existing outputs; clear them for reruns.
            rm -f {output.collapsed_vcf} {output.removed_vcf}
            truvari collapse -i {input.concat_vcf} -o {output.collapsed_vcf} -c {output.removed_vcf} \
                -f {input.ref} -k {params.keep} \
                --pctseq {params.pctseq} --pctsize {params.pctsize} --refdist {params.refdist} \
                --sizemin {params.sizemin} --sizemax {params.sizemax} 2>> {log}
            python3 {params.script:q} \
                --collapsed {output.collapsed_vcf} --removed {output.removed_vcf} \
                --output {output.support_vcf} 2>> {log}
            echo "NUM_CALLERS distribution:" >> {log}
            grep -v "^#" {output.support_vcf} \
                | sed -n 's/.*NUM_CALLERS=\\([0-9]*\\).*/\\1/p' \
                | sort -n | uniq -c >> {log} 2>&1 || true
            """

    rule ngs_sv_collapse_merge_sample:
        """Filter by distinct-caller support and length, keep standard chromosomes, compress (bcftools)."""
        input:
            support_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/collapse/{{sample}}.collapsed.support.vcf"
        output:
            merged_vcf_gz = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.ngs_sv.merged.vcf.gz",
            tbi = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.ngs_sv.merged.vcf.gz.tbi"
        params:
            min_callers = max(1, NGS_MERGE_MIN_CALLERS),
            min_length = NGS_COLLAPSE_MIN_LENGTH
        threads: 1
        resources:
            mem_mb = 4 * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["tabix_module"],
            config["tool_modules"]["bcftools_module"]
        conda:
            ENV_BASIC
        container:
            config["containers"]["bcftools"]
        log:
            f"{LOGS_DIR}/ngs_sv/collapse_merge/{{sample}}.merge.log"
        shell:
            """
            set -o pipefail
            : > {log}
            SORT_TMP=$(dirname {output.merged_vcf_gz})/tmp_sort_$$
            mkdir -p "$SORT_TMP"
            # Breakends have no SVLEN and are kept, as in the SURVIVOR chain.
            FILTER_EXPR='INFO/NUM_CALLERS >= {params.min_callers} && (INFO/SVLEN=="." || abs(INFO/SVLEN) >= {params.min_length})'
            awk 'BEGIN{{OFS="\t"}} /^#/{{print;next}} $1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/ {{print}}' \
                {input.support_vcf} \
                | bcftools view -i "$FILTER_EXPR" -Ou 2>> {log} \
                | bcftools sort -T "$SORT_TMP" -Oz -o {output.merged_vcf_gz} 2>> {log}
            tabix -f -p vcf {output.merged_vcf_gz} 2>> {log}
            rm -rf "$SORT_TMP"
            echo "Collapsed representatives: $(grep -vc '^#' {input.support_vcf} || true)" >> {log}
            echo "Min distinct callers: {params.min_callers}; min length: {params.min_length}" >> {log}
            echo "SVTYPE counts after filtering:" >> {log}
            bcftools query -f '%INFO/SVTYPE\\n' {output.merged_vcf_gz} | sort | uniq -c >> {log}
            """

def _ngs_sv_population_inputs(wildcards):
    vcfs = []
    for s in get_ngs_samples():
        vcfs.append(f"{NGS_SV_SINGLE_DIR}/{s}/merged/{s}.ngs_sv.merged.vcf.gz")
    return vcfs

rule ngs_sv_survivor_population_prepare:
    """Prepare merge file list for population SURVIVOR merge."""
    input:
        vcfs = _ngs_sv_population_inputs
    output:
        merge_list = f"{NGS_SV_COHORT_DIR}/merged/merge_files.txt"
    threads: 1
    resources:
        mem_mb = 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    localrule: True
    log:
        f"{LOGS_DIR}/ngs_sv/survivor_population_merge/cohort.prepare.log"
    shell:
        """
        mkdir -p $(dirname {output.merge_list})
        DECOMPDIR=$(dirname {output.merge_list})/decompressed
        mkdir -p "$DECOMPDIR"
        : > {output.merge_list}
        for vcf in {input.vcfs}; do
            # SURVIVOR 1.0.7 cannot read gzipped VCF — decompress first
            if [[ "$vcf" == *.vcf.gz ]]; then
                plain_vcf="$DECOMPDIR/$(basename "$vcf" .gz)"
                zcat "$vcf" > "$plain_vcf"
                echo "$plain_vcf" >> {output.merge_list}
            else
                echo "$vcf" >> {output.merge_list}
            fi
        done
        echo "Prepared merge list with $(wc -l < {output.merge_list}) VCF files" >> {log}
        """

rule ngs_sv_survivor_population_run:
    """Run population-level SURVIVOR merge. Uses survivor container."""
    input:
        merge_list = f"{NGS_SV_COHORT_DIR}/merged/merge_files.txt",
        vcfs = _ngs_sv_population_inputs  # ensure VCFs exist
    output:
        merged_vcf = f"{NGS_SV_COHORT_DIR}/merged/cohort.ngs_sv.merged.vcf"
    params:
        max_distance = SURV_MAX_DISTANCE,
        min_support = SURV_POP_SUPPORT,
        min_length = SURV_MIN_LENGTH
    threads: 1
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["survivor_module"]
    conda:
        ENV_SURVIVOR
    container:
        config["containers"]["survivor"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/survivor_population_merge/cohort.tsv"
    log:
        f"{LOGS_DIR}/ngs_sv/survivor_population_merge/cohort.log"
    shell:
        """
        echo "SURVIVOR which: $(command -v SURVIVOR)" >> {log} 2>&1 || true
        SURVIVOR --version >> {log} 2>&1 || true
        
        SURVIVOR merge \
            {input.merge_list} \
            {params.max_distance} {params.min_support} 1 1 0 {params.min_length} \
            {output.merged_vcf} >> {log} 2>&1
        """

rule ngs_sv_survivor_population_merge:
    """Validate, clean, compress and index population SURVIVOR output. Uses bcftools container."""
    input:
        merged_vcf = f"{NGS_SV_COHORT_DIR}/merged/cohort.ngs_sv.merged.vcf"
    output:
        merged_vcf_gz = f"{NGS_SV_COHORT_DIR}/merged/cohort.ngs_sv.merged.vcf.gz",
        tbi = f"{NGS_SV_COHORT_DIR}/merged/cohort.ngs_sv.merged.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 4 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["tabix_module"],
        config["tool_modules"]["bcftools_module"]
    conda:
        ENV_BASIC
    container:
        config["containers"]["bcftools"]
    log:
        f"{LOGS_DIR}/ngs_sv/survivor_population_merge/cohort.compress.log"
    shell:
        """
        # Validate and sanitize SURVIVOR output (known to produce corrupted records)
        TMP_VCF_CLEAN=$(dirname {output.merged_vcf_gz})/cohort.ngs_sv.merged.clean.vcf
        echo "Info: Validating and cleaning SURVIVOR output" >> {log}
        # Keep only valid records: proper CHROM, numeric POS, valid REF/ALT, minimum fields
        awk 'BEGIN{{OFS="\t"}} 
            /^##/ {{print; next}} 
            /^#CHROM/ {{print; next}} 
            {{if (NF >= 10 && $1 ~ /^(chr)?[0-9XYM]+$/ && $2 ~ /^[0-9]+$/ && length($4) > 0 && length($5) > 0 && $4 != "." && $5 != ".") print}}' "{input.merged_vcf}" > "$TMP_VCF_CLEAN" || true
        
        # Check if we have valid variants after cleaning (use || true to avoid strict mode failure)
        VALID_VARIANTS=$(grep -v "^#" "$TMP_VCF_CLEAN" 2>/dev/null | wc -l || echo "0")
        echo "Info: Found $VALID_VARIANTS valid variants after sanitization" >> {log}
        
        if [ "$VALID_VARIANTS" -eq 0 ]; then
            echo "Warning: No valid variants after SURVIVOR merge and sanitization. Creating minimal VCF." >> {log}
            # Create a minimal valid VCF with headers only
            grep "^#" "{input.merged_vcf}" > "$TMP_VCF_CLEAN" || true
        fi
        
        bgzip -f -c "$TMP_VCF_CLEAN" > {output.merged_vcf_gz}
        tabix -f -p vcf {output.merged_vcf_gz}
        rm -f "$TMP_VCF_CLEAN"

        # Clean up decompressed VCFs created by population_prepare
        rm -rf $(dirname {output.merged_vcf_gz})/decompressed 2>/dev/null || true
        """

# Aggregate per-sample SV targets for convenience

def ngs_sv_final_outputs():
    outs = []
    for s in get_ngs_samples():
        # Manta outputs if enabled
        if config.get("ngs_sv", {}).get("manta", {}).get("enabled", True):
            outs.extend([
                f"{NGS_SV_SINGLE_DIR}/{s}/manta/{s}.manta.vcf.gz",
                f"{NGS_SV_SINGLE_DIR}/{s}/manta/{s}.manta.doctored.vcf.gz",
                f"{NGS_SV_SINGLE_DIR}/{s}/manta/{s}.manta.counts.1.txt",
            ])
        # Smoove/Lumpy outputs if enabled
        if config.get("ngs_sv", {}).get("smoove", {}).get("enabled", True):
            outs.extend([
                f"{NGS_SV_SINGLE_DIR}/{s}/smoove/{s}-smoove.genotyped.vcf.gz",
                f"{NGS_SV_SINGLE_DIR}/{s}/lumpy/count_lumpy/{s}.lumpy.counts.1.txt",
            ])
        # Standalone callers do not require merge or annotation.
        if config.get("ngs_sv", {}).get("merge_per_sample", True):
            outs.extend([
                f"{NGS_SV_SINGLE_DIR}/{s}/merged/{s}.ngs_sv.merged.vcf.gz",
                f"{NGS_SV_SINGLE_DIR}/{s}/merged/{s}.ngs_sv.merged.vcf.gz.tbi",
            ])
            if config.get("annotsv", {}).get("enabled", True):
                outs.append(f"{NGS_SV_SINGLE_DIR}/{s}/annotation/{s}.ngs.annotsv.tsv")
        # Add CNVpytor outputs if enabled
        if config.get("ngs_sv", {}).get("cnvpytor", {}).get("enabled", False):
            outs.extend([
                f"{NGS_SV_SINGLE_DIR}/{s}/cnvpytor/{s}.summary.txt",
            ])
            # Add CNVpytor VCF if converter provided
            if str(config.get("ngs_sv", {}).get("cnvpytor", {}).get("to_vcf_script", "")).strip() != "":
                outs.extend([
                    f"{NGS_SV_SINGLE_DIR}/{s}/cnvpytor/{s}.cnvpytor.vcf.gz",
                    f"{NGS_SV_SINGLE_DIR}/{s}/cnvpytor/{s}.cnvpytor.vcf.gz.tbi",
                ])
    # Add population/cohort-level merged outputs only if enabled and we have at least one sample
    if SURV_ENABLE_POP_MERGE and config.get("ngs_sv", {}).get("merge_per_sample", True) and len(get_ngs_samples()) > 0:
        outs.extend([
            f"{NGS_SV_COHORT_DIR}/merged/cohort.ngs_sv.merged.vcf.gz",
            f"{NGS_SV_COHORT_DIR}/merged/cohort.ngs_sv.merged.vcf.gz.tbi",
        ])
    return outs
