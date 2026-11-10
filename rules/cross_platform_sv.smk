# Cross-platform (NGS + TGS) Structural Variant integration
#
# Purpose: for samples that were sequenced on BOTH platforms, merge the per-sample
# short-read (NGS) and long-read (TGS) SV call sets with SURVIVOR to produce:
#   1. a merged VCF carrying SURVIVOR's INFO/SUPP_VEC (which platform(s) support each SV)
#   2. a concordance summary (NGS-only / TGS-only / shared counts, per SV type)
#
# This is an OPT-IN validation module (default off). It only CONSUMES the outputs of
# the existing NGS and TGS SV workflows and never modifies them. The set of samples to
# integrate is defined by the user-maintained mapping table loaded in
# workflow/samples_cross_platform.smk (see get_cross_platform_pairs()).
#
# Reuses the same containers and the proven "decompress -> sanitize -> SURVIVOR merge
# -> clean/compress" pattern as rules/ngs_sv.smk and rules/tgs_sv*.smk.

import os

# Per-rule conda environments (reuse existing definitions)
ENV_BASIC = "../envs/basic.yaml"
ENV_SURVIVOR = "../envs/ngs_sv/ngs_sv.survivor.yaml"

# Input SV call sets produced by the existing per-platform workflows
NGS_SV_SINGLE_DIR = config["output"]["ngs"]["sv"]["single_sample"]
TGS_SV_SINGLE_DIR = config["output"]["tgs"]["sv"]["single_sample"]
LOGS_DIR = config["output"]["logs"]

# Output directory for cross-platform results. Falls back to results/cross_platform/...
# so the module works even when the config file does not declare output.cross_platform.
CROSS_PLATFORM_SV_DIR = (
    config["output"].get("cross_platform", {}).get("sv")
    if isinstance(config["output"].get("cross_platform"), dict)
    else None
) or os.path.join(config["output"]["root"], "cross_platform", "structural_variants")

# Absolute path to the concordance script, anchored to the workflow directory so it
# resolves correctly regardless of the current working directory (per-subfolder runs).
CROSS_SV_SCRIPT = os.path.join(
    workflow.basedir, "scripts", "cross_platform", "sv_concordance.py"
)

# SURVIVOR merge parameters for cross-platform integration
_CROSS_SV_CFG = config.get("cross_platform_sv", {})
CROSS_SV_MAX_DISTANCE = _CROSS_SV_CFG.get("max_distance", 1000)
CROSS_SV_MIN_LENGTH = _CROSS_SV_CFG.get("min_length", 30)
# Union merge: keep every site so NGS-only / TGS-only calls remain visible for concordance
CROSS_SV_MIN_SUPPORT = 1
# Short-read and long-read callers often disagree on SVTYPE representation, so by default
# we do NOT require identical type when matching across platforms.
CROSS_SV_TYPE_FLAG = 1 if _CROSS_SV_CFG.get("require_same_type", False) else 0
CROSS_SV_STRAND_FLAG = 1 if _CROSS_SV_CFG.get("require_same_strand", False) else 0

# SUPP_VEC positions follow the order the inputs are written into the merge list below.
CROSS_SV_NGS_INDEX = 0
CROSS_SV_TGS_INDEX = 1


def _cross_ngs_vcf(pair_id):
    pair = get_cross_platform_pair(pair_id)
    return f"{NGS_SV_SINGLE_DIR}/{pair['ngs_sample']}/merged/{pair['ngs_sample']}.ngs_sv.merged.vcf.gz"


def _cross_tgs_vcf(pair_id):
    pair = get_cross_platform_pair(pair_id)
    return f"{TGS_SV_SINGLE_DIR}/{pair['tgs_sample']}/final/{pair['tgs_sample']}.sv.final.vcf.gz"


# Constrain the wildcard to known pair_ids to avoid ambiguous matches.
wildcard_constraints:
    pair="|".join([p["pair_id"] for p in get_cross_platform_pairs()]) or "__no_cross_platform_pair__"


rule cross_platform_sv_sanitize:
    """Decompress and sanitize the NGS and TGS SV VCFs for SURVIVOR.

    SURVIVOR 1.0.7 cannot read gzipped input and chokes on REF='.', multi-valued
    INFO/SVLEN (reduced to its first value), and '<DUP:TANDEM>' ALTs, so we normalize them first (same treatment
    as the per-platform SURVIVOR rules). The merge list is written NGS-first, TGS-second
    so INFO/SUPP_VEC position 0 = NGS and position 1 = TGS."""
    input:
        ngs_vcf = lambda w: _cross_ngs_vcf(w.pair),
        tgs_vcf = lambda w: _cross_tgs_vcf(w.pair),
    output:
        merge_list = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/merge_files.txt",
        done_flag = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/.sanitize_done",
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
        f"{LOGS_DIR}/cross_platform_sv/{{pair}}.sanitize.log"
    shell:
        """
        OUTDIR=$(dirname {output.merge_list})
        mkdir -p "$OUTDIR"
        : > {output.merge_list}

        # Sanitize a single VCF and append the plain-text result to the merge list.
        # $1 = source vcf(.gz), $2 = label used for the sanitized filename (ngs/tgs).
        sanitize_one() {{
            local src="$1"
            local label="$2"
            local out="$OUTDIR/${{label}}_sanitized.vcf"
            echo "Info: sanitizing $src -> $out" >> {log}
            zcat "$src" | awk 'BEGIN{{OFS="\t"}}
                /^#/ {{print; next}}
                {{
                    if ($4 == ".") $4 = "N"
                    gsub(/<DUP:TANDEM>/, "<DUP>", $5)
                    # Keep the first integer SVLEN value: SURVIVOR derives DEL/INS/DUP
                    # length from SVLEN and drops them via min_length if it is missing.
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
                }}' > "$out"
            echo "$out" >> {output.merge_list}
        }}

        # Order matters: NGS first (SUPP_VEC index 0), TGS second (index 1).
        sanitize_one "{input.ngs_vcf}" ngs
        sanitize_one "{input.tgs_vcf}" tgs

        if [ "$(wc -l < {output.merge_list})" -lt 2 ]; then
            echo "Error: expected 2 sanitized VCFs for cross-platform merge" >> {log}
            exit 1
        fi
        touch {output.done_flag}
        """


rule cross_platform_sv_merge_run:
    """Run SURVIVOR merge on the sanitized NGS+TGS VCFs (union; support>=1)."""
    input:
        merge_list = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/merge_files.txt",
        done_flag = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/.sanitize_done",
    output:
        merged_vcf = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/{{pair}}.cross_sv.merged.vcf",
    params:
        max_distance = CROSS_SV_MAX_DISTANCE,
        min_support = CROSS_SV_MIN_SUPPORT,
        type_flag = CROSS_SV_TYPE_FLAG,
        strand_flag = CROSS_SV_STRAND_FLAG,
        min_length = CROSS_SV_MIN_LENGTH,
    threads: 1
    resources:
        mem_mb = 8 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["survivor_module"]
    conda:
        ENV_SURVIVOR
    container:
        config["containers"]["survivor"]
    log:
        f"{LOGS_DIR}/cross_platform_sv/{{pair}}.merge.log"
    shell:
        """
        echo "SURVIVOR which: $(command -v SURVIVOR)" >> {log} 2>&1 || true
        SURVIVOR --version >> {log} 2>&1 || true

        SURVIVOR merge \
            {input.merge_list} \
            {params.max_distance} {params.min_support} {params.type_flag} {params.strand_flag} 0 {params.min_length} \
            {output.merged_vcf} >> {log} 2>&1
        """


rule cross_platform_sv_finalize:
    """Validate, clean, compress and index the SURVIVOR output; preserves SUPP_VEC."""
    input:
        merged_vcf = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/{{pair}}.cross_sv.merged.vcf",
    output:
        merged_vcf_gz = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/{{pair}}.cross_sv.merged.vcf.gz",
        tbi = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/{{pair}}.cross_sv.merged.vcf.gz.tbi",
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
        f"{LOGS_DIR}/cross_platform_sv/{{pair}}.finalize.log"
    shell:
        """
        TMP_VCF_CLEAN=$(dirname {output.merged_vcf_gz})/{wildcards.pair}.cross_sv.merged.clean.vcf
        echo "Info: validating and cleaning SURVIVOR output" >> {log}
        awk 'BEGIN{{OFS="\t"}}
            /^##/ {{print; next}}
            /^#CHROM/ {{print; next}}
            {{if (NF >= 10 && $1 ~ /^(chr)?[0-9XYM]+$/ && $2 ~ /^[0-9]+$/ && length($4) > 0 && length($5) > 0 && $4 != "." && $5 != ".") print}}' "{input.merged_vcf}" > "$TMP_VCF_CLEAN" || true

        VALID_VARIANTS=$(grep -v "^#" "$TMP_VCF_CLEAN" 2>/dev/null | wc -l || echo "0")
        echo "Info: found $VALID_VARIANTS valid variants after sanitization" >> {log}
        if [ "$VALID_VARIANTS" -eq 0 ]; then
            echo "Warning: no valid variants; creating header-only VCF" >> {log}
            grep "^#" "{input.merged_vcf}" > "$TMP_VCF_CLEAN" || true
        fi

        bgzip -f -c "$TMP_VCF_CLEAN" > {output.merged_vcf_gz}
        tabix -f -p vcf {output.merged_vcf_gz}
        rm -f "$TMP_VCF_CLEAN"

        # Remove intermediate sanitized VCFs
        rm -f $(dirname {output.merged_vcf_gz})/ngs_sanitized.vcf \
              $(dirname {output.merged_vcf_gz})/tgs_sanitized.vcf 2>/dev/null || true
        """


rule cross_platform_sv_concordance:
    """Summarize NGS/TGS SV concordance from the merged VCF's SUPP_VEC field."""
    input:
        merged_vcf_gz = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/{{pair}}.cross_sv.merged.vcf.gz",
    output:
        report = f"{CROSS_PLATFORM_SV_DIR}/{{pair}}/{{pair}}.concordance.tsv",
    params:
        script = CROSS_SV_SCRIPT,
        pair_id = "{pair}",
        ngs_index = CROSS_SV_NGS_INDEX,
        tgs_index = CROSS_SV_TGS_INDEX,
    threads: 1
    resources:
        mem_mb = 2 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        ENV_BASIC
    container:
        config["containers"]["cyvcf2"]
    log:
        f"{LOGS_DIR}/cross_platform_sv/{{pair}}.concordance.log"
    shell:
        """
        mkdir -p $(dirname {output.report})
        python {params.script} \
            --input {input.merged_vcf_gz} \
            --output {output.report} \
            --pair-id {params.pair_id} \
            --ngs-index {params.ngs_index} \
            --tgs-index {params.tgs_index} > {log} 2>&1
        """
