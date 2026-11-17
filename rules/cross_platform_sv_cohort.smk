# Cohort-level cross-platform SV integration (opt-in).
#
# Purpose: integrate a short-read (NGS) cohort SV VCF with a long-read (TGS) cohort SV VCF
# into one site-level catalogue that records, for every SV site, which platform(s) detected
# it and the allele frequency in each cohort. The two cohorts do NOT need to share
# individuals (e.g. a large short-read cohort plus a smaller long-read cohort).
#
# Steps:
#   1. sites:    per platform, recompute AN/AC/AF from genotypes, drop genotypes, drop
#                non-carrier sites (AC < min_ac), give every site a unique ID and keep the
#                source ID in INFO/SRC_ID.
#   2. merge:    SURVIVOR merge (default) or truvari collapse (method: truvari). Matching is
#                by breakpoint distance and size only, because short-read SVs carry no
#                sequence (truvari runs with --pctseq 0).
#   3. annotate: write a sorted site-only VCF with INFO/PLATFORM (both / short_only /
#                long_only), INFO/SUPP_VEC (NGS, TGS) and per-platform AF/AC/AN, plus a
#                summary table.
#
# Inputs default to this run's cohort outputs; set cross_platform_sv_cohort.ngs_vcf /
# tgs_vcf to integrate cohort VCFs produced elsewhere (for example in separate runs).

import os

_CPC_CFG = config.get("cross_platform_sv_cohort", {})
CPC_METHOD = str(_CPC_CFG.get("method", "survivor")).strip().lower()
if CPC_METHOD not in ("survivor", "truvari"):
    raise ValueError(f"cross_platform_sv_cohort.method must be 'survivor' or 'truvari', got '{CPC_METHOD}'")
CPC_MAX_DISTANCE = _CPC_CFG.get("max_distance", 1000)
CPC_MIN_LENGTH = _CPC_CFG.get("min_length", 50)
CPC_MIN_AC = _CPC_CFG.get("min_ac", 1)
CPC_PASS_ONLY = _CPC_CFG.get("pass_only", True)
CPC_REQUIRE_SAME_TYPE = _CPC_CFG.get("require_same_type", False)
CPC_REPORT_MAF = _CPC_CFG.get("report_maf", 0.01)
_CPC_TRUVARI = _CPC_CFG.get("truvari", {})
CPC_TRUVARI_PCTSIZE = _CPC_TRUVARI.get("pctsize", 0.7)
CPC_TRUVARI_SIZEMAX = _CPC_TRUVARI.get("sizemax", 100000000)

CPC_DIR = (
    config["output"].get("cross_platform", {}).get("sv_cohort")
    if isinstance(config["output"].get("cross_platform"), dict)
    else None
) or os.path.join(config["output"]["root"], "cross_platform", "cohort_sv")
CPC_LOG_DIR = f"{config['output']['logs']}/cross_platform_sv_cohort"
CPC_SCRIPT_DIR = os.path.join(workflow.basedir, "scripts", "cross_platform")
CPC_ENV_BASIC = "../envs/basic.yaml"
CPC_ENV_SURVIVOR = "../envs/ngs_sv/ngs_sv.survivor.yaml"
CPC_ENV_TRUVARI = "../envs/long_read/ont_truvari.yaml"


def _cpc_input_vcf(platform):
    """Cohort SV VCF for one platform: explicit path, or this run's cohort output."""
    path = str(_CPC_CFG.get(f"{platform}_vcf", "") or "").strip()
    if path:
        return path
    if platform == "ngs":
        if not config["workflows"].get("enable_ngs_sv_detection", False):
            raise ValueError("cross_platform_sv_cohort: set ngs_vcf or enable workflows.enable_ngs_sv_detection")
        return f"{config['output']['ngs']['sv']['cohort']}/merged/cohort.ngs_sv.merged.vcf.gz"
    if not config["workflows"].get("enable_tgs_sv_detection", False):
        raise ValueError("cross_platform_sv_cohort: set tgs_vcf or enable workflows.enable_tgs_sv_detection")
    return f"{config['output']['tgs']['sv']['cohort']}/population.merged.vcf"


wildcard_constraints:
    cpc_platform = "ngs|tgs"


rule cross_platform_sv_cohort_sites:
    """Per-platform site-only VCF with recomputed AN/AC/AF and unique site IDs."""
    input:
        vcf = lambda w: _cpc_input_vcf(w.cpc_platform)
    output:
        sites = f"{CPC_DIR}/sites/{{cpc_platform}}.sites.vcf"
    params:
        label = lambda w: w.cpc_platform.upper(),
        filter_opt = "-f PASS,." if CPC_PASS_ONLY else "",
        min_ac = CPC_MIN_AC,
        awk = os.path.join(CPC_SCRIPT_DIR, "prepare_sv_sites.awk")
    threads: 2
    resources:
        mem_mb = 8 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["bcftools_module"]
    conda:
        CPC_ENV_BASIC
    container:
        config["containers"]["bcftools"]
    benchmark:
        f"{config['output']['logs']}/benchmarks/cross_platform_sv_cohort/{{cpc_platform}}.sites.tsv"
    log:
        f"{CPC_LOG_DIR}/{{cpc_platform}}.sites.log"
    shell:
        """
        set -o pipefail
        mkdir -p $(dirname {output.sites})
        TMP=$(dirname {output.sites})/{wildcards.cpc_platform}.sort_tmp
        rm -rf "$TMP"; mkdir -p "$TMP"
        # Genotypes are dropped right after AF is computed, so the sort and every later
        # step only handle site-level records, whatever the cohort size.
        bcftools view {params.filter_opt} --threads {threads} -Ou {input.vcf} 2> {log} \
            | bcftools +fill-tags -Ou -- -t AN,AC,AF 2>> {log} \
            | bcftools view -i 'INFO/AC>={params.min_ac}' -G -Ou 2>> {log} \
            | bcftools sort -T "$TMP" -Ov 2>> {log} \
            | awk -v label={params.label} -f {params.awk:q} > {output.sites}
        rm -rf "$TMP"
        echo "Sites kept: $(grep -vc '^#' {output.sites} || true)" >> {log}
        """


if CPC_METHOD == "survivor":
    rule cross_platform_sv_cohort_merge:
        """Match NGS and TGS cohort sites with SURVIVOR (union; NGS first, TGS second)."""
        input:
            ngs = f"{CPC_DIR}/sites/ngs.sites.vcf",
            tgs = f"{CPC_DIR}/sites/tgs.sites.vcf"
        output:
            merged = f"{CPC_DIR}/merge/cohort.survivor.vcf"
        params:
            max_distance = CPC_MAX_DISTANCE,
            type_flag = 1 if CPC_REQUIRE_SAME_TYPE else 0,
            min_length = CPC_MIN_LENGTH
        threads: 1
        resources:
            mem_mb = 16 * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["survivor_module"]
        conda:
            CPC_ENV_SURVIVOR
        container:
            config["containers"]["survivor"]
        benchmark:
            f"{config['output']['logs']}/benchmarks/cross_platform_sv_cohort/merge.tsv"
        log:
            f"{CPC_LOG_DIR}/merge.log"
        shell:
            """
            mkdir -p $(dirname {output.merged})
            LIST=$(dirname {output.merged})/merge_files.txt
            printf '%s\\n' {input.ngs} {input.tgs} > "$LIST"
            # Arguments: max distance, min support (1 = union), same type, same strand,
            # estimate distance from SV size (0 = off), min SV length.
            SURVIVOR merge "$LIST" {params.max_distance} 1 {params.type_flag} 0 0 {params.min_length} \
                {output.merged} > {log} 2>&1
            """
else:
    rule cross_platform_sv_cohort_concat:
        """Concatenate both site files for truvari; QUAL ranks TGS above NGS so that the
        sequence-resolved long-read record is kept as the representative."""
        input:
            ngs = f"{CPC_DIR}/sites/ngs.sites.vcf",
            tgs = f"{CPC_DIR}/sites/tgs.sites.vcf"
        output:
            vcf = f"{CPC_DIR}/merge/cohort.concat.vcf.gz",
            tbi = f"{CPC_DIR}/merge/cohort.concat.vcf.gz.tbi"
        threads: 1
        resources:
            mem_mb = 8 * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        envmodules:
            config["tool_modules"]["tabix_module"],
            config["tool_modules"]["bcftools_module"]
        conda:
            CPC_ENV_BASIC
        container:
            config["containers"]["bcftools"]
        log:
            f"{CPC_LOG_DIR}/concat.log"
        shell:
            """
            set -o pipefail
            OUT=$(dirname {output.vcf}); mkdir -p "$OUT"
            : > {log}
            for pair in "ngs:1:{input.ngs}" "tgs:2:{input.tgs}"; do
                label=${{pair%%:*}}; rest=${{pair#*:}}; qual=${{rest%%:*}}; src=${{rest#*:}}
                # Same placeholder sample name in both files, as bcftools concat requires.
                awk -v q=$qual 'BEGIN{{OFS="\\t"}} /^#CHROM/{{$10="SITES"}} /^#/{{print; next}} {{$6=q; print}}' "$src" \
                    | bcftools view -Oz -o "$OUT/$label.concat_in.vcf.gz" 2>> {log}
                tabix -f -p vcf "$OUT/$label.concat_in.vcf.gz"
            done
            bcftools concat -a "$OUT/ngs.concat_in.vcf.gz" "$OUT/tgs.concat_in.vcf.gz" -Oz -o {output.vcf} 2>> {log}
            tabix -f -p vcf {output.vcf}
            rm -f "$OUT"/ngs.concat_in.vcf.gz* "$OUT"/tgs.concat_in.vcf.gz*
            """

    rule cross_platform_sv_cohort_merge:
        """Match NGS and TGS cohort sites with truvari collapse (no sequence comparison)."""
        input:
            vcf = f"{CPC_DIR}/merge/cohort.concat.vcf.gz",
            tbi = f"{CPC_DIR}/merge/cohort.concat.vcf.gz.tbi",
            ref = config["reference"]["genome"]
        output:
            merged = f"{CPC_DIR}/merge/cohort.collapsed.vcf",
            removed = f"{CPC_DIR}/merge/cohort.collapsed.removed.vcf"
        params:
            refdist = CPC_MAX_DISTANCE,
            pctsize = CPC_TRUVARI_PCTSIZE,
            sizemin = CPC_MIN_LENGTH,
            sizemax = CPC_TRUVARI_SIZEMAX,
            type_opt = "" if CPC_REQUIRE_SAME_TYPE else "--typeignore"
        threads: 1
        resources:
            mem_mb = 16 * 1024,
            partition = config["cluster"]["partition"],
            time = config["cluster"]["time"]
        conda:
            CPC_ENV_TRUVARI
        container:
            config["containers"].get("truvari", "docker://quay.io/biocontainers/truvari:5.4.0--pyhdfd78af_0")
        benchmark:
            f"{config['output']['logs']}/benchmarks/cross_platform_sv_cohort/merge.tsv"
        log:
            f"{CPC_LOG_DIR}/merge.log"
        shell:
            """
            rm -f {output.merged} {output.removed}
            truvari collapse -i {input.vcf} -o {output.merged} -c {output.removed} \
                -f {input.ref} -k maxqual {params.type_opt} --pctseq 0 \
                --pctsize {params.pctsize} --refdist {params.refdist} \
                --sizemin {params.sizemin} --sizemax {params.sizemax} > {log} 2>&1
            """


rule cross_platform_sv_cohort_annotate:
    """Write the integrated site-only VCF (PLATFORM, SUPP_VEC, per-platform AF) and summary."""
    input:
        merged = f"{CPC_DIR}/merge/cohort.survivor.vcf" if CPC_METHOD == "survivor" else f"{CPC_DIR}/merge/cohort.collapsed.vcf",
        removed = [] if CPC_METHOD == "survivor" else f"{CPC_DIR}/merge/cohort.collapsed.removed.vcf",
        ngs = f"{CPC_DIR}/sites/ngs.sites.vcf",
        tgs = f"{CPC_DIR}/sites/tgs.sites.vcf",
        fai = f"{config['reference']['genome']}.fai"
    output:
        vcf = f"{CPC_DIR}/cohort.cross_sv.vcf.gz",
        tbi = f"{CPC_DIR}/cohort.cross_sv.vcf.gz.tbi",
        summary = f"{CPC_DIR}/cohort.cross_sv.summary.tsv"
    params:
        script = os.path.join(CPC_SCRIPT_DIR, "annotate_cohort_cross_sv.py"),
        method = CPC_METHOD,
        removed_opt = lambda w, input: f"--removed {input.removed}" if input.removed else "",
        report_maf = CPC_REPORT_MAF
    threads: 1
    resources:
        mem_mb = 8 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["tabix_module"]
    conda:
        CPC_ENV_BASIC
    container:
        config["containers"]["cyvcf2"]
    log:
        f"{CPC_LOG_DIR}/annotate.log"
    shell:
        """
        set -o pipefail
        python3 {params.script:q} --method {params.method} --merged {input.merged} {params.removed_opt} \
            --ngs-sites {input.ngs} --tgs-sites {input.tgs} --report-maf {params.report_maf} \
            --fai {input.fai} --summary {output.summary} 2> {log} | bgzip -c > {output.vcf}
        tabix -f -p vcf {output.vcf}
        """
