# NGS cohort-level SV integration with svtools (ngs_sv.cohort_merge_method: svtools).
# Follows the svtools/LUMPY cohort workflow: lsort -> lmerge -> per-sample re-genotyping
# of the merged sites (svtools genotype, i.e. SVTyper) -> vcfpaste -> afreq -> prune.
# svtools only accepts LUMPY-style VCFs (STRANDS, PRPOS/PREND), so the input is the
# per-sample Smoove VCF rather than the multi-caller merged VCF used by SURVIVOR.
# svtools and SVTyper ship in the Smoove container and Conda environment.
# Produces the same cohort.ngs_sv.merged.vcf.gz as the SURVIVOR population chain.
# Included from rules/ngs_sv.smk, which defines ENV_SMOOVE, ENV_BASIC and get_ngs_samples().

_SVTOOLS_CFG = config.get("ngs_sv", {}).get("svtools", {})
SVTOOLS_LSORT_BATCH = _SVTOOLS_CFG.get("lsort_batch_size", 200)
SVTOOLS_LMERGE_SLOP = _SVTOOLS_CFG.get("lmerge_slop", 20)
SVTOOLS_GENOTYPE_PAIRS = _SVTOOLS_CFG.get("genotype_insert_pairs", 1000000)
SVTOOLS_GENOTYPE_MEM_MB = _SVTOOLS_CFG.get("genotype_mem_mb", 8192)
SVTOOLS_PRUNE_DISTANCE = _SVTOOLS_CFG.get("prune_distance", 100)

SVTOOLS_DIR = f"{NGS_SV_COHORT_DIR}/svtools"
SVTOOLS_LOG_DIR = f"{LOGS_DIR}/ngs_sv/svtools_cohort"

ruleorder: ngs_sv_svtools_cohort_finalize > ngs_sv_survivor_population_merge


def _svtools_smoove_vcfs(wildcards):
    return [
        f"{NGS_SV_SINGLE_DIR}/{s}/smoove/{s}-smoove.genotyped.vcf.gz"
        for s in get_ngs_samples()
    ]


def _svtools_genotyped_vcfs(wildcards):
    return [f"{SVTOOLS_DIR}/genotyped/{s}.vcf.gz" for s in get_ngs_samples()]


rule ngs_sv_svtools_cohort_merge:
    """Sort and merge per-sample Smoove calls into cohort SV sites (svtools lsort + lmerge)."""
    input:
        vcfs = _svtools_smoove_vcfs
    output:
        sites = f"{SVTOOLS_DIR}/cohort.sites.vcf.gz"
    params:
        batch = SVTOOLS_LSORT_BATCH,
        slop = SVTOOLS_LMERGE_SLOP
    threads: 1
    resources:
        mem_mb = 16 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["smoove_module"]
    conda:
        ENV_SMOOVE
    container:
        config["containers"]["smoove"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/svtools_cohort/merge.tsv"
    log:
        f"{SVTOOLS_LOG_DIR}/merge.log"
    shell:
        """
        set -o pipefail
        export PYTHONNOUSERSITE=1
        OUTDIR=$(dirname {output.sites})
        TMP="$OUTDIR/tmp_lsort_$$"
        mkdir -p "$TMP"
        LIST="$OUTDIR/smoove_vcfs.txt"
        printf '%s\\n' {input.vcfs} > "$LIST"
        svtools --version > {log} 2>&1
        echo "Input Smoove VCFs: $(wc -l < "$LIST")" >> {log}
        svtools lsort -f "$LIST" -t "$TMP" -b {params.batch} > "$TMP/sorted.vcf" 2>> {log}
        svtools lmerge -i "$TMP/sorted.vcf" -f {params.slop} 2>> {log} | bgzip -c > {output.sites}
        echo "Sorted calls: $(grep -vc '^#' "$TMP/sorted.vcf" || true)" >> {log}
        echo "Merged sites: $(zcat {output.sites} | grep -vc '^#' || true)" >> {log}
        rm -rf "$TMP"
        """


rule ngs_sv_svtools_cohort_genotype:
    """Re-genotype the cohort SV sites in one sample (svtools genotype / SVTyper)."""
    input:
        sites = f"{SVTOOLS_DIR}/cohort.sites.vcf.gz",
        bam = f"{NGS_ALIGNMENT_DIR}/{{sample}}/bam/{{sample}}.dup.sorted.recal.bam",
        bai = f"{NGS_ALIGNMENT_DIR}/{{sample}}/bam/{{sample}}.dup.sorted.recal.bai"
    output:
        vcf = f"{SVTOOLS_DIR}/genotyped/{{sample}}.vcf.gz",
        lib_info = f"{SVTOOLS_DIR}/genotyped/{{sample}}.lib_info.json"
    params:
        pairs = SVTOOLS_GENOTYPE_PAIRS
    threads: 1
    resources:
        mem_mb = SVTOOLS_GENOTYPE_MEM_MB,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["smoove_module"],
        config["tool_modules"]["svtyper_module"]
    conda:
        ENV_SMOOVE
    container:
        config["containers"]["smoove"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/svtools_cohort/genotype/{{sample}}.tsv"
    log:
        f"{SVTOOLS_LOG_DIR}/genotype/{{sample}}.log"
    shell:
        """
        set -o pipefail
        export PYTHONNOUSERSITE=1
        rm -f {output.lib_info}
        # QUAL is overwritten by svtools genotype with this sample's evidence.
        # SVTyper writes FORMAT "GT" alone for sites it cannot score; vcfpaste takes FORMAT
        # from one sample only, so pad those records to the full SVTyper FORMAT.
        zcat {input.sites} \
            | svtools genotype -B {input.bam} -l {output.lib_info} -n {params.pairs} 2> {log} \
            | awk 'BEGIN{{OFS="\\t"; full="GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB"; n=split(full, f, ":")}}
                   /^#/ {{print; next}}
                   $9 == "GT" {{$9 = full; for (i = 2; i <= n; i++) $10 = $10 ":."; padded++}}
                   {{print}}
                   END {{print "Padded GT-only records: " padded + 0 > "/dev/stderr"}}' 2>> {log} \
            | bgzip -c > {output.vcf}
        """


rule ngs_sv_svtools_cohort_paste:
    """Paste per-sample genotypes, add allele frequencies and prune redundant sites."""
    input:
        sites = f"{SVTOOLS_DIR}/cohort.sites.vcf.gz",
        vcfs = _svtools_genotyped_vcfs
    output:
        vcf = f"{SVTOOLS_DIR}/cohort.pruned.vcf.gz"
    params:
        prune_distance = SVTOOLS_PRUNE_DISTANCE
    threads: 1
    resources:
        mem_mb = 16 * 1024,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    envmodules:
        config["tool_modules"]["smoove_module"]
    conda:
        ENV_SMOOVE
    container:
        config["containers"]["smoove"]
    benchmark:
        f"{LOGS_DIR}/benchmarks/ngs_sv/svtools_cohort/paste.tsv"
    log:
        f"{SVTOOLS_LOG_DIR}/paste.log"
    shell:
        """
        set -o pipefail
        export PYTHONNOUSERSITE=1
        OUTDIR=$(dirname {output.vcf})
        LIST="$OUTDIR/genotyped_vcfs.txt"
        printf '%s\\n' {input.vcfs} > "$LIST"
        BEDPE="$OUTDIR/cohort.pasted.sorted.bedpe"
        : > {log}
        svtools vcfpaste -m {input.sites} -f "$LIST" -q 2>> {log} \
            | svtools afreq 2>> {log} \
            | svtools vcftobedpe 2>> {log} \
            | svtools bedpesort > "$BEDPE" 2>> {log}
        svtools prune -s -d {params.prune_distance} -e AF "$BEDPE" 2>> {log} \
            | svtools bedpetovcf 2>> {log} \
            | svtools vcfsort 2>> {log} \
            | bgzip -c > {output.vcf}
        rm -f "$BEDPE"
        """


rule ngs_sv_svtools_cohort_finalize:
    """Drop bulky LUMPY breakpoint probabilities, sort, compress and index the cohort VCF."""
    input:
        vcf = f"{SVTOOLS_DIR}/cohort.pruned.vcf.gz",
        fai = f"{REFERENCE_FA}.fai"
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
        f"{SVTOOLS_LOG_DIR}/finalize.log"
    shell:
        """
        set -o pipefail
        SORT_TMP=$(dirname {output.merged_vcf_gz})/tmp_sort_$$
        mkdir -p "$SORT_TMP"
        # svtools writes ##contig lines in reverse order; rebuild them in reference order
        # (bcftools reheader --fai keeps the existing order) so bcftools sort follows it.
        bcftools annotate -x INFO/PRPOS,INFO/PREND -Ov {input.vcf} 2> {log} \
            | awk 'BEGIN{{OFS="\\t"}}
                   NR == FNR {{contigs = contigs "##contig=<ID=" $1 ",length=" $2 ">\\n"; next}}
                   /^##contig=/ {{next}}
                   /^#CHROM/ {{printf "%s", contigs}}
                   {{print}}' {input.fai} - \
            | bcftools sort -T "$SORT_TMP" -Oz -o {output.merged_vcf_gz} 2>> {log}
        tabix -f -p vcf {output.merged_vcf_gz} 2>> {log}
        rm -rf "$SORT_TMP"
        echo "Cohort SV records: $(bcftools index -n {output.merged_vcf_gz})" >> {log}
        echo "SVTYPE counts:" >> {log}
        bcftools query -f '%INFO/SVTYPE\\n' {output.merged_vcf_gz} | sort | uniq -c >> {log}
        """
