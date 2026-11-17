# Aggregate QC and statistics outputs with MultiQC.

import os
import shlex

ENV_MULTIQC = "../envs/multiqc.yaml"
CONTAINER_MULTIQC = config.get("containers", {}).get("multiqc", "")

MULTIQC_CFG = config.get("multiqc", {})
MULTIQC_OUTDIR = MULTIQC_CFG.get(
    "outdir",
    config.get("output", {}).get("multiqc", os.path.join(config["output"]["root"], "multiqc")),
)
MULTIQC_FILENAME = MULTIQC_CFG.get("filename", "multiqc_report.html")
MULTIQC_BASENAME = os.path.splitext(MULTIQC_FILENAME)[0]
MULTIQC_TITLE = MULTIQC_CFG.get("title", "WGS Pipeline QC Report")
MULTIQC_EXTRA_ARGS = MULTIQC_CFG.get("extra_args", "")
MULTIQC_EXPORT_PLOTS = MULTIQC_CFG.get("export_plots", True)
MULTIQC_MODULES = MULTIQC_CFG.get("modules", [])


def get_multiqc_input_files(wildcards=None):
    """Return concrete upstream files that should exist before MultiQC runs."""
    if not config["workflows"].get("enable_multiqc", False):
        return []

    inputs = []

    if config["workflows"].get("enable_basic_ngs_workflows", False):
        sample_ids = get_all_sample_ids()
        sample_ids_with_fastq = get_sample_ids_with_fastq()

        if config.get("fastp", {}).get("enabled", False):
            inputs.extend([
                f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastp/{sample_id}.fastp.json"
                for sample_id in sample_ids_with_fastq
            ])

        inputs.extend([
            f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastqc/{sample_id}_{read}_fastqc.zip"
            for sample_id in sample_ids_with_fastq
            for read in ["1", "2"]
        ])

        for sample_id in sample_ids:
            inputs.extend([
                f"{NGS_ALIGNMENT_DIR}/{sample_id}/qc/qualimap/genome_results.txt",
                f"{NGS_ALIGNMENT_DIR}/{sample_id}/qc/qualimap/qualimapReport.html",
            ])
            if config.get("verifybamid2", {}).get("enabled", True):
                inputs.append(f"{NGS_ALIGNMENT_DIR}/{sample_id}/qc/verifybamid/{sample_id}.selfSM")

        if config["workflows"].get("enable_ngs_single_sample_filtering", False):
            inputs.extend([
                f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/variant_calling/normalized/{sample_id}.filtered.normalized.stats.txt"
                for sample_id in sample_ids
            ])

        if (
            config["workflows"].get("enable_ngs_single_sample_filtering", False)
            and config.get("ngs_single_sample_variants", {}).get("annotation", {}).get("tool", "funcotator") == "snpeff"
        ):
            inputs.extend([
                f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/variant_calling/annotated/{sample_id}.snpeff_summary.html"
                for sample_id in sample_ids
            ])

    if config["workflows"].get("enable_ngs_joint_calling", False):
        inputs.append(f"{NGS_VARIANTS_COHORT_DIR}/filtered/cohort.filtered.recode.stats.txt")

    if config["workflows"].get("enable_tgs_sv_detection", False):
        inputs.extend([
            f"{TGS_PREPROCESSING_DIR}/{sample}/qc/nanoplot/{sample}NanoPlot-report.html"
            for sample in get_ont_samples()
        ])

    if config["workflows"].get("enable_tgs_variant_calling", False):
        selected_tool = get_tgs_variant_calling_tool()
        if selected_tool == "clair3":
            inputs.extend([
                f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/clair3/{sample}.clair3.stats.txt"
                for sample in get_ont_samples()
            ])
        elif selected_tool == "deepvariant":
            inputs.extend([
                f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/deepvariant/{sample}.deepvariant.stats.txt"
                for sample in get_ont_samples()
            ])

    if config["workflows"].get("enable_tgs_variant_annotation", False):
        tgs_annotation_tool = config.get("tgs_variant_annotation", {}).get("tool", "snpeff").strip().lower()
        if tgs_annotation_tool == "snpeff":
            selected_tool = get_tgs_variant_calling_tool()
            inputs.extend([
                f"{TGS_VARIANTS_SINGLE_DIR}/{sample}/{selected_tool}/annotated/{sample}.snpeff_summary.html"
                for sample in get_ont_samples()
            ])

    return inputs


def get_multiqc_search_paths():
    """Return directories that MultiQC should scan for sidecar files and plots."""
    paths = []

    if config["workflows"].get("enable_basic_ngs_workflows", False):
        paths.extend([NGS_PREPROCESSING_DIR, NGS_ALIGNMENT_DIR])

    if config["workflows"].get("enable_tgs_sv_detection", False):
        paths.append(TGS_PREPROCESSING_DIR)

    if config["workflows"].get("enable_tgs_variant_calling", False):
        paths.append(TGS_VARIANTS_SINGLE_DIR)

    if (
        config["workflows"].get("enable_ngs_single_sample_filtering", False)
        or config["workflows"].get("enable_tgs_variant_annotation", False)
    ):
        paths.extend([NGS_VARIANTS_SINGLE_DIR, TGS_VARIANTS_SINGLE_DIR])

    return paths


def get_multiqc_outputs():
    """Return MultiQC outputs for rule all."""
    if not config["workflows"].get("enable_multiqc", False):
        return []
    if not get_multiqc_input_files():
        return []

    outputs = [
        f"{MULTIQC_OUTDIR}/{MULTIQC_FILENAME}",
        f"{MULTIQC_OUTDIR}/{MULTIQC_BASENAME}_data",
        f"{MULTIQC_OUTDIR}/multiqc_inputs.txt",
    ]
    if MULTIQC_EXPORT_PLOTS:
        outputs.append(f"{MULTIQC_OUTDIR}/{MULTIQC_BASENAME}_plots")
    return outputs


def _multiqc_module_args():
    if not MULTIQC_MODULES:
        return ""
    return " ".join(f"--module {shlex.quote(str(module))}" for module in MULTIQC_MODULES)


rule multiqc_report:
    input:
        get_multiqc_input_files
    output:
        report = f"{MULTIQC_OUTDIR}/{MULTIQC_FILENAME}",
        data_dir = directory(f"{MULTIQC_OUTDIR}/{MULTIQC_BASENAME}_data"),
        # MultiQC only writes the plot directory when --export is used.
        **({"plots_dir": directory(f"{MULTIQC_OUTDIR}/{MULTIQC_BASENAME}_plots")}
           if MULTIQC_EXPORT_PLOTS else {}),
        manifest = f"{MULTIQC_OUTDIR}/multiqc_inputs.txt"
    params:
        outdir = MULTIQC_OUTDIR,
        title = MULTIQC_TITLE,
        extra_args = MULTIQC_EXTRA_ARGS,
        module_args = _multiqc_module_args(),
        export_plots = MULTIQC_EXPORT_PLOTS
    threads: 2
    resources:
        mem_mb = 4096,
        partition = config["cluster"]["partition"],
        time = "04:00:00"
    envmodules:
        config["tool_modules"].get("multiqc_module", "")
    conda:
        ENV_MULTIQC
    container:
        CONTAINER_MULTIQC
    log:
        f"{LOGS_DIR}/multiqc/multiqc.log"
    benchmark:
        f"{LOGS_DIR}/benchmarks/multiqc/multiqc.txt"
    run:
        os.makedirs(params.outdir, exist_ok=True)
        os.makedirs(os.path.dirname(str(log)), exist_ok=True)

        manifest_paths = [str(path) for path in input]
        manifest_paths.extend([
            path for path in get_multiqc_search_paths()
            if os.path.exists(path)
        ])

        seen = set()
        unique_paths = []
        for path in manifest_paths:
            if path not in seen:
                seen.add(path)
                unique_paths.append(path)

        with open(output.manifest, "w") as manifest:
            for path in unique_paths:
                manifest.write(f"{path}\n")

        export_arg = "--export" if params.export_plots else ""
        command = " ".join([
            "multiqc",
            "--force",
            "--data-dir",
            export_arg,
            "--file-list", shlex.quote(output.manifest),
            "--outdir", shlex.quote(params.outdir),
            "--filename", shlex.quote(MULTIQC_FILENAME),
            "--title", shlex.quote(params.title),
            params.module_args,
            params.extra_args,
        ]).strip()

        shell(
            f"{command} > {log} 2>&1"
        )
