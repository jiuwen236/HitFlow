# Optional HitSV 2.0.4 SRS, LRS and matched-individual Hybrid calling.
import hashlib
import json
import re
import shlex

_HITSV_CFG = config.get("hitsv", {})
_HITSV_INVALID_POLICY = _HITSV_CFG.get("invalid_allele_policy", "error")
if _HITSV_INVALID_POLICY not in ("error", "quarantine"):
    raise ValueError("hitsv.invalid_allele_policy must be error or quarantine")
_HITSV_SRS = config["workflows"].get("enable_ngs_sv_detection", False) and config.get("ngs_sv", {}).get("hitsv", {}).get("enabled", False)
_HITSV_LRS = config["workflows"].get("enable_tgs_sv_detection", False) and config.get("tgs_config", {}).get("sv_detection", {}).get("hitsv", {}).get("enabled", False)
_HITSV_HYBRID = _HITSV_CFG.get("hybrid", {}).get("enabled", False)
_HITSV_SOURCE = str(workflow.basedir)
_HITSV_SCRIPT = os.path.join(_HITSV_SOURCE, "scripts/sv_pipeline_scripts/hitsv.py")
_HITSV_BIN = _HITSV_CFG.get("binary", os.path.join(_HITSV_SOURCE, "bin/HitSV/2.0.4/HitSV"))
_HITSV_REF = config["resolved_paths"]["reference_genome"]
_HITSV_WORK = _HITSV_CFG.get("work_dir", os.path.join(config["output"]["root"], "hitsv_work"))
_HITSV_SRS_WORK = _HITSV_CFG.get("srs_work_dir", _HITSV_WORK)
_HITSV_CHROMS = [str(c) for c in _HITSV_CFG.get("chromosomes", config["chromosomes"])]
_HITSV_ENV = "../envs/hitsv.yaml"
_HITSV_CONTAINER = config["containers"].get("hitsv", config["containers"].get("truvari", "docker://quay.io/biocontainers/truvari:5.4.0--pyhdfd78af_0"))
_HITSV_PRESET = _HITSV_CFG.get("lrs_preset", "HIFI" if config["tgs_config"]["platform"] == "PACBIO_HIFI" else "ERR_PRONE")
_HITSV_ANNOTATE = _HITSV_CFG.get("annotate", False) and config.get("annotsv", {}).get("enabled", True)
_HITSV_ANNOT_DIR = os.path.join(config["output"]["root"], "annotations", "hitsv")
_HITSV_ERROR_PROFILE = _HITSV_CFG.get("hybrid", {}).get("error_profile", os.path.join(_HITSV_SOURCE, "bin/HitSV/2.0.4/9mer_profile"))
_HITSV_HYBRID_DIR = _HITSV_CFG.get("hybrid", {}).get("output_dir", os.path.join(config["output"]["root"], "cross_platform", "hitsv"))
if not _HITSV_CHROMS or len(set(_HITSV_CHROMS)) != len(_HITSV_CHROMS):
    raise ValueError("hitsv.chromosomes must be a nonempty unique list")
if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", c) for c in _HITSV_CHROMS):
    raise ValueError("HitSV chromosome names must be safe path components")
# Cache identity includes reference metadata; Snakemake also tracks FASTA/FAI inputs.
_ref_stat = os.stat(_HITSV_REF)
_ref_key = hashlib.sha256((os.path.realpath(_HITSV_REF) + str(_ref_stat.st_size) + str(_ref_stat.st_mtime_ns) + "2.0.4").encode()).hexdigest()[:20]
_HITSV_FA_STAT = os.path.join(_HITSV_CFG.get("cache_dir", _HITSV_WORK + "/cache"), _ref_key + ".stat.txt")
_HITSV_PAIRS = {}
if _HITSV_HYBRID:
    for pair in CROSS_PLATFORM_PAIRS:
        if pair["ngs_sample"] not in get_all_sample_ids() or pair["tgs_sample"] not in get_ont_samples():
            raise ValueError("Unresolved HitSV Hybrid pair: " + pair["pair_id"])
        _HITSV_PAIRS[pair["pair_id"]] = pair
    if not _HITSV_PAIRS:
        raise ValueError("HitSV Hybrid requires samples_pairs_csv with matched samples")


def _hitsv_bam(sample, platform):
    supplied = _HITSV_CFG.get(platform + "_bams", {}).get(sample)
    if supplied:
        return supplied
    if platform == "srs":
        return f"{NGS_ALIGNMENT_DIR}/{sample}/bam/{sample}.dup.sorted.recal.bam"
    suffix = "ont" if config["tgs_config"]["platform"] == "ONT" else "pb"
    return f"{TGS_PREPROCESSING_DIR}/{sample}/alignment/{sample}.{suffix}.sorted.bam"


def _hitsv_index(bam):
    candidates = [bam + ".bai", bam[:-4] + ".bai"]
    return next((p for p in candidates if os.path.exists(p)), candidates[1] if ".recal.bam" in bam else candidates[0])


def _hitsv_samples(mode, sample):
    if mode == "hybrid":
        pair = _HITSV_PAIRS[sample]
        return [("srs", pair["ngs_sample"]), ("lrs", pair["tgs_sample"])]
    if mode not in ("srs", "lrs"):
        raise ValueError("Unknown HitSV mode: " + mode)
    return [(mode, sample)]


def _hitsv_bam_inputs(w):
    paths = {}
    for mode, sample in _hitsv_samples(w.mode, w.sample):
        paths[mode] = _hitsv_bam(sample, mode)
        paths[mode + "_index"] = _hitsv_index(paths[mode])
    return paths


def _hitsv_preflight_args(w):
    args = []
    for mode, sample in _hitsv_samples(w.mode, w.sample):
        bam = _hitsv_bam(sample, mode)
        args += ["--bam", bam, "--index", _hitsv_index(bam)]
    return shlex.join(args)


def _hitsv_srs_sample(w):
    return _HITSV_PAIRS[w.sample]["ngs_sample"] if w.mode == "hybrid" else w.sample


def _hitsv_call_inputs(w):
    inputs = _hitsv_bam_inputs(w)
    if w.mode != "lrs":
        sample = _hitsv_srs_sample(w)
        base = f"{_HITSV_SRS_WORK}/srs/{sample}"
        inputs.update(fa_stat=_HITSV_FA_STAT, bam_stat=base + "/bam.stat.json",
                      tl=base + "/TL.sorted.bam", tl_index=base + "/TL.sorted.bam.bai")
    if w.mode == "hybrid" and _HITSV_ERROR_PROFILE:
        inputs["error_profile"] = _HITSV_ERROR_PROFILE
    return inputs


def _hitsv_call_args(w):
    args = []
    inputs = _hitsv_call_inputs(w)
    for key in ("srs", "lrs", "fa_stat", "bam_stat", "tl", "error_profile"):
        if key in inputs:
            args += ["--" + key.replace("_", "-"), inputs[key]]
    return shlex.join(args)


def hitsv_outputs():
    outputs = []
    groups = []
    if _HITSV_SRS:
        groups.append((NGS_SV_SINGLE_DIR, get_all_sample_ids(), True, "srs"))
    if _HITSV_LRS:
        groups.append((TGS_SV_SINGLE_DIR, get_ont_samples(), True, "lrs"))
    if _HITSV_HYBRID:
        groups.append((_HITSV_HYBRID_DIR, _HITSV_PAIRS, False, "hybrid"))
    for root, samples, nested, mode in groups:
        for sample in samples:
            path = f"{root}/{sample}/" + ("hitsv/" if nested else "") + f"{sample}.hitsv.vcf.gz"
            outputs.extend([path, path + ".tbi"])
            if _HITSV_ANNOTATE:
                outputs.append(f"{_HITSV_ANNOT_DIR}/{mode}/{sample}.hitsv.annotsv.tsv")
    return outputs


rule hitsv_preflight:
    input:
        unpack(_hitsv_bam_inputs),
        ref=_HITSV_REF,
        fai=_HITSV_REF + ".fai",
        script=_HITSV_SCRIPT
    output:
        manifest=f"{_HITSV_WORK}/{{mode}}/{{sample}}/manifest.json"
    wildcard_constraints:
        mode="srs|lrs|hybrid"
    params:
        bams=_hitsv_preflight_args,
        chromosomes=",".join(_HITSV_CHROMS),
        identities=json.dumps(_HITSV_CFG.get("hybrid", {}).get("sample_identities", {}))
    resources:
        mem_mb=2048,
        hitsv_io=1,
        partition=config["cluster"]["partition"],
        time="00:15:00"
    conda: _HITSV_ENV
    container: _HITSV_CONTAINER
    log: f"{LOGS_DIR}/hitsv/{{mode}}/{{sample}}.preflight.log"
    shell:
        """
        python3 {input.script:q} preflight {params.bams} --reference {input.ref:q} \
          --chromosomes {params.chromosomes:q} --identities {params.identities:q} \
          --output {output.manifest:q} > {log:q} 2>&1
        """


rule hitsv_reference_stat:
    input:
        ref=_HITSV_REF,
        fai=_HITSV_REF + ".fai",
        binary=_HITSV_BIN
    output: stat=_HITSV_FA_STAT
    threads: 1
    resources:
        mem_mb=int(_HITSV_CFG.get("reference_stat_mem_mb", _HITSV_CFG.get("preprocess_mem_mb", 1024))),
        hitsv_io=1,
        partition=config["cluster"]["partition"],
        time=config["cluster"]["time"]
    conda: _HITSV_ENV
    container: _HITSV_CONTAINER
    benchmark: f"{LOGS_DIR}/benchmarks/hitsv/reference.{_ref_key}.tsv"
    log: f"{LOGS_DIR}/hitsv/reference.{_ref_key}.log"
    shell:
        """
        HITSV_STAT_TMP=$(mktemp {output.stat:q}.partial.XXXXXX)
        trap 'rm -f "$HITSV_STAT_TMP"' EXIT
        {input.binary:q} fa_stat {input.ref:q} > "$HITSV_STAT_TMP" 2> {log:q}
        test -s "$HITSV_STAT_TMP"
        mv "$HITSV_STAT_TMP" {output.stat:q}
        """


rule hitsv_srs_bam_stat:
    input:
        bam=lambda w: _hitsv_bam(w.sample, "srs"),
        bai=lambda w: _hitsv_index(_hitsv_bam(w.sample, "srs")),
        manifest=f"{_HITSV_WORK}/srs/{{sample}}/manifest.json",
        ref=_HITSV_REF,
        fai=_HITSV_REF + ".fai",
        binary=_HITSV_BIN
    output: stat=f"{_HITSV_WORK}/srs/{{sample}}/bam.stat.json"
    threads: 1
    resources:
        mem_mb=int(_HITSV_CFG.get("bam_stat_mem_mb", _HITSV_CFG.get("preprocess_mem_mb", 1024))),
        hitsv_io=1,
        partition=config["cluster"]["partition"],
        time=config["cluster"]["time"]
    conda: _HITSV_ENV
    container: _HITSV_CONTAINER
    benchmark: f"{LOGS_DIR}/benchmarks/hitsv/srs/{{sample}}.bam_stat.tsv"
    log: f"{LOGS_DIR}/hitsv/srs/{{sample}}.bam_stat.log"
    shell: "{input.binary:q} srs_bam_stat {input.ref:q} {input.bam:q} > {output.stat:q} 2> {log:q}"


rule hitsv_srs_trans_reads:
    input:
        bam=lambda w: _hitsv_bam(w.sample, "srs"),
        bai=lambda w: _hitsv_index(_hitsv_bam(w.sample, "srs")),
        manifest=f"{_HITSV_WORK}/srs/{{sample}}/manifest.json",
        ref=_HITSV_REF,
        fai=_HITSV_REF + ".fai",
        binary=_HITSV_BIN
    output: bam=temp(f"{_HITSV_WORK}/srs/{{sample}}/TL.bam")
    threads: 1
    resources:
        mem_mb=int(_HITSV_CFG.get("trans_reads_mem_mb", _HITSV_CFG.get("preprocess_mem_mb", 1024))),
        hitsv_io=1,
        partition=config["cluster"]["partition"],
        time=config["cluster"]["time"]
    conda: _HITSV_ENV
    container: _HITSV_CONTAINER
    benchmark: f"{LOGS_DIR}/benchmarks/hitsv/srs/{{sample}}.trans_reads.tsv"
    log: f"{LOGS_DIR}/hitsv/srs/{{sample}}.trans_reads.log"
    shell: "{input.binary:q} srs_trans_reads {input.ref:q} {input.bam:q} {output.bam:q} > {log:q} 2>&1"


rule hitsv_sort_trans_reads:
    input: bam=f"{_HITSV_WORK}/srs/{{sample}}/TL.bam"
    output:
        bam=f"{_HITSV_WORK}/srs/{{sample}}/TL.sorted.bam",
        bai=f"{_HITSV_WORK}/srs/{{sample}}/TL.sorted.bam.bai"
    threads: int(_HITSV_CFG.get("sort_threads", 2))
    resources:
        mem_mb=lambda wildcards, threads: samtools_sort_memory_mb("1G", threads),
        hitsv_io=1,
        partition=config["cluster"]["partition"],
        time=config["cluster"]["time"]
    conda: _HITSV_ENV
    container: config["containers"]["samtools"]
    benchmark: f"{LOGS_DIR}/benchmarks/hitsv/srs/{{sample}}.sort_tl.tsv"
    log: f"{LOGS_DIR}/hitsv/srs/{{sample}}.sort_tl.log"
    shell:
        """
        samtools sort -@ {threads} -m 1G -T {output.bam:q}.tmp -o {output.bam:q} {input.bam:q} 2> {log:q}
        samtools index -@ {threads} {output.bam:q} {output.bai:q} 2>> {log:q}
        """


rule hitsv_call_chromosome:
    input:
        unpack(_hitsv_call_inputs),
        manifest=f"{_HITSV_WORK}/{{mode}}/{{sample}}/manifest.json",
        ref=_HITSV_REF,
        fai=_HITSV_REF + ".fai",
        binary=_HITSV_BIN,
        script=_HITSV_SCRIPT
    output: vcf=f"{_HITSV_WORK}/{{mode}}/{{sample}}/raw/{{chromosome}}.vcf"
    wildcard_constraints:
        mode="srs|lrs|hybrid",
        chromosome="|".join(re.escape(c) for c in _HITSV_CHROMS)
    params:
        options=_hitsv_call_args,
        preset=_HITSV_PRESET,
        alignment_preset=_HITSV_CFG.get("alignment_preset", "asm5"),
        min_size=int(_HITSV_CFG.get("min_sv_length", 50))
    threads: 1
    resources:
        mem_mb=int(_HITSV_CFG.get("call_mem_mb", 16384)),
        hitsv_io=1,
        partition=config["cluster"]["partition"],
        time=config["cluster"]["time"]
    conda: _HITSV_ENV
    container: _HITSV_CONTAINER
    benchmark: f"{LOGS_DIR}/benchmarks/hitsv/{{mode}}/{{sample}}.{{chromosome}}.call.tsv"
    log: f"{LOGS_DIR}/hitsv/{{mode}}/{{sample}}.{{chromosome}}.call.log"
    shell:
        """
        python3 {input.script:q} call --binary {input.binary:q} --manifest {input.manifest:q} \
          --chromosome {wildcards.chromosome:q} --reference {input.ref:q} --output {output.vcf:q} \
          --preset {params.preset:q} --alignment-preset {params.alignment_preset:q} \
          --min-size {params.min_size} {params.options} > {log:q} 2>&1
        """


rule hitsv_normalize:
    input:
        vcfs=lambda w: expand(f"{_HITSV_WORK}/{{mode}}/{{sample}}/raw/{{chromosome}}.vcf", mode=w.mode, sample=w.sample, chromosome=_HITSV_CHROMS),
        ref=_HITSV_REF,
        fai=_HITSV_REF + ".fai",
        script=_HITSV_SCRIPT
    output:
        vcf=f"{_HITSV_WORK}/{{mode}}/{{sample}}/{{sample}}.normalized.vcf",
        audit=f"{_HITSV_WORK}/{{mode}}/{{sample}}/vcf_audit.json",
        rejected=f"{_HITSV_WORK}/{{mode}}/{{sample}}/quarantined.jsonl"
    params:
        min_size=int(_HITSV_CFG.get("min_sv_length", 50)),
        repair="--repair-padding" if _HITSV_CFG.get("repair_padding", False) else "",
        policy=_HITSV_INVALID_POLICY
    resources:
        mem_mb=4096,
        partition=config["cluster"]["partition"],
        time="01:00:00"
    conda: _HITSV_ENV
    container: _HITSV_CONTAINER
    log: f"{LOGS_DIR}/hitsv/{{mode}}/{{sample}}.normalize.log"
    shell:
        """
        python3 {input.script:q} normalize --inputs {input.vcfs:q} --sample {wildcards.sample:q} \
          --output {output.vcf:q} --audit {output.audit:q} --reference {input.ref:q} --min-size {params.min_size} \
          {params.repair} --invalid-allele-policy {params.policy:q} --rejected {output.rejected:q} > {log:q} 2>&1
        """


rule hitsv_compress:
    input: vcf=f"{_HITSV_WORK}/{{mode}}/{{sample}}/{{sample}}.normalized.vcf"
    output:
        vcf=f"{_HITSV_WORK}/{{mode}}/{{sample}}/{{sample}}.hitsv.vcf.gz",
        tbi=f"{_HITSV_WORK}/{{mode}}/{{sample}}/{{sample}}.hitsv.vcf.gz.tbi"
    resources:
        mem_mb=4096,
        partition=config["cluster"]["partition"],
        time="01:00:00"
    conda: _HITSV_ENV
    container: config["containers"]["bcftools"]
    log: f"{LOGS_DIR}/hitsv/{{mode}}/{{sample}}.compress.log"
    shell:
        """
        bcftools sort -m 1G -T {output.vcf:q}.sort -Oz -o {output.vcf:q} {input.vcf:q} 2> {log:q}
        bcftools index --tbi --force {output.vcf:q} 2>> {log:q}
        """


rule hitsv_publish_srs:
    input:
        vcf=f"{_HITSV_WORK}/srs/{{sample}}/{{sample}}.hitsv.vcf.gz",
        tbi=f"{_HITSV_WORK}/srs/{{sample}}/{{sample}}.hitsv.vcf.gz.tbi"
    output:
        vcf=f"{NGS_SV_SINGLE_DIR}/{{sample}}/hitsv/{{sample}}.hitsv.vcf.gz",
        tbi=f"{NGS_SV_SINGLE_DIR}/{{sample}}/hitsv/{{sample}}.hitsv.vcf.gz.tbi"
    resources:
        mem_mb=1024,
        partition=config["cluster"]["partition"],
        time="00:10:00"
    shell:
        """
        cp {input.vcf:q} {output.vcf:q}
        cp {input.tbi:q} {output.tbi:q}
        """


rule hitsv_publish_lrs:
    input:
        vcf=f"{_HITSV_WORK}/lrs/{{sample}}/{{sample}}.hitsv.vcf.gz",
        tbi=f"{_HITSV_WORK}/lrs/{{sample}}/{{sample}}.hitsv.vcf.gz.tbi"
    output:
        vcf=f"{TGS_SV_SINGLE_DIR}/{{sample}}/hitsv/{{sample}}.hitsv.vcf.gz",
        tbi=f"{TGS_SV_SINGLE_DIR}/{{sample}}/hitsv/{{sample}}.hitsv.vcf.gz.tbi"
    resources:
        mem_mb=1024,
        partition=config["cluster"]["partition"],
        time="00:10:00"
    shell:
        """
        cp {input.vcf:q} {output.vcf:q}
        cp {input.tbi:q} {output.tbi:q}
        """


rule hitsv_publish_hybrid:
    input:
        vcf=f"{_HITSV_WORK}/hybrid/{{sample}}/{{sample}}.hitsv.vcf.gz",
        tbi=f"{_HITSV_WORK}/hybrid/{{sample}}/{{sample}}.hitsv.vcf.gz.tbi"
    output:
        vcf=f"{_HITSV_HYBRID_DIR}/{{sample}}/{{sample}}.hitsv.vcf.gz",
        tbi=f"{_HITSV_HYBRID_DIR}/{{sample}}/{{sample}}.hitsv.vcf.gz.tbi"
    resources:
        mem_mb=1024,
        partition=config["cluster"]["partition"],
        time="00:10:00"
    shell:
        """
        cp {input.vcf:q} {output.vcf:q}
        cp {input.tbi:q} {output.tbi:q}
        """


# Optional caller-level annotation, including standalone Hybrid runs.
rule hitsv_annotsv:
    input:
        vcf=f"{_HITSV_WORK}/{{mode}}/{{sample}}/{{sample}}.normalized.vcf"
    output:
        tsv=f"{_HITSV_ANNOT_DIR}/{{mode}}/{{sample}}.hitsv.annotsv.tsv"
    wildcard_constraints:
        mode="srs|lrs|hybrid"
    params:
        annotations=config.get("annotsv", {}).get("annotations_dir", "ref_data/AnnotSV_annotations"),
        genome=config.get("annotsv", {}).get("genome_build", "GRCh38"),
        min_size=int(_HITSV_CFG.get("min_sv_length", 50))
    threads: 1
    resources:
        mem_mb=8192,
        partition=config["cluster"]["partition"],
        time=config["cluster"]["time"]
    conda: "../envs/sv_annotation/annotsv.yaml"
    container: config.get("containers", {}).get("annotsv", "")
    benchmark: f"{LOGS_DIR}/benchmarks/hitsv/{{mode}}/{{sample}}.annotsv.tsv"
    log: f"{LOGS_DIR}/hitsv/{{mode}}/{{sample}}.annotsv.log"
    shell:
        """
        export LANG=C LC_ALL=C
        ANNOTSV_TASK_DIR=$(mktemp -d {output.tsv:q}.work.XXXXXX)
        trap 'rm -rf "$ANNOTSV_TASK_DIR"' EXIT
        bcftools view -e 'INFO/SVTYPE="BND" || INFO/SVTYPE="TRA"' -Ov \
          -o "$ANNOTSV_TASK_DIR/input.vcf" {input.vcf:q} 2> {log:q}
        if ! grep -qv '^#' "$ANNOTSV_TASK_DIR/input.vcf"; then
            printf 'AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\n' > {output.tsv:q}
        else
            AnnotSV -SVinputFile "$ANNOTSV_TASK_DIR/input.vcf" \
              -outputDir "$ANNOTSV_TASK_DIR" -outputFile result \
              -annotationsDir {params.annotations:q} -genomeBuild {params.genome:q} \
              -SVminSize {params.min_size} >> {log:q} 2>&1
            test -s "$ANNOTSV_TASK_DIR/result.tsv"
            cp "$ANNOTSV_TASK_DIR/result.tsv" {output.tsv:q}
        fi
        """
