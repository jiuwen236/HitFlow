# WGS Pipeline Snakemake workflow
# Modular structure for better maintainability

import os

# Include common configuration and helper modules
include: "workflow/common.smk"       # Configuration and directory definitions
include: "workflow/helpers.smk"      # Helper functions
include: "workflow/samples_ngs.smk"  # NGS sample discovery and management
include: "workflow/samples_tgs.smk"  # TGS sample discovery and management
include: "workflow/samples_cross_platform.smk"  # Cross-platform (NGS+TGS) sample pairing


# ---------- Resolved input file paths (setup before rule includes) ----------
RESOLVED_PATHS = {
    "reference_genome": resolve_input_file_path(config["reference"]["genome"]),
}

# VerifyBamID2 SVD prefix (only resolve if enabled)
_verifybamid2_cfg = config.get("verifybamid2", {})
if _verifybamid2_cfg.get("enabled", True):
    _svd_prefix = _verifybamid2_cfg.get("svd_prefix", config.get("svd_prefix", ""))
    if _svd_prefix:
        RESOLVED_PATHS["svd_prefix"] = resolve_input_file_path(_svd_prefix + ".bed")[:-4]

_known_sites_cfg = str(config.get("gatk", {}).get("known_sites_vcf", "")).strip()
RESOLVED_PATHS["known_sites_vcf"] = (
    resolve_input_file_path(_known_sites_cfg)
    if _known_sites_cfg else ""
)

if config["vqsr"]["enabled"]:
    RESOLVED_PATHS["vqsr_resources"] = {
        "hapmap": resolve_input_file_path(config["vqsr"]["resources"]["hapmap"]),
        "omni": resolve_input_file_path(config["vqsr"]["resources"]["omni"]),
        "1000g_snps": resolve_input_file_path(config["vqsr"]["resources"]["1000g_snps"]),
        "dbsnp": resolve_input_file_path(config["vqsr"]["resources"]["dbsnp"]),
        "mills_indels": resolve_input_file_path(config["vqsr"]["resources"]["mills_indels"]),
    }

# Resolve WES targets if provided
if "wes" in config and str(config["wes"].get("targets_bed", "")).strip() != "":
    RESOLVED_PATHS["wes"] = {
        "targets_bed": resolve_input_file_path(config["wes"]["targets_bed"]),
        "interval_padding": int(config["wes"].get("interval_padding", 0)),
    }

if config["phasing"]["enabled"]:
    RESOLVED_PATHS["phasing"] = {
        "genetic_map": config["phasing"]["shapeit"]["genetic_map"]
    }

tgs_phasing_cfg = config.get("tgs_joint_calling", {}).get("phasing", {})
if (
    config["workflows"].get("enable_tgs_phasing", False)
    and tgs_phasing_cfg.get("enabled", False)
    and tgs_phasing_cfg.get("shapeit", {}).get("genetic_map")
):
    RESOLVED_PATHS["tgs_phasing"] = {
        "genetic_map": tgs_phasing_cfg["shapeit"]["genetic_map"]
    }

config["resolved_paths"] = RESOLVED_PATHS

# HitSV is opt-in and can run without downstream merge/annotation targets.
if (config.get("ngs_sv", {}).get("hitsv", {}).get("enabled", False)
    or config.get("tgs_config", {}).get("sv_detection", {}).get("hitsv", {}).get("enabled", False)
    or config.get("hitsv", {}).get("hybrid", {}).get("enabled", False)):
    include: "rules/hitsv.smk"

# ---------- Include workflow rules based on configuration ----------
if config["workflows"]["enable_basic_ngs_workflows"]:
    if config["fastp"]["enabled"]:
        include: "rules/ngs_fastp.smk"
    include: "rules/ngs_align.smk"
    include: "rules/ngs_bqsr.smk"
    include: "rules/ngs_qc.smk"
    include: "rules/ngs_variant_calling.smk"

if config["workflows"].get("enable_ngs_single_sample_filtering", False):
    include: "rules/ngs_single_sample_filtering.smk"
    annotation_tool = config.get("ngs_single_sample_variants", {}).get("annotation", {}).get("tool", "funcotator")
    if annotation_tool == "funcotator":
        include: "rules/ngs_funcotator.smk"
    elif annotation_tool == "snpeff":
        include: "rules/ngs_snpeff.smk"

if config["workflows"]["enable_ngs_joint_calling"]:
    include: "rules/ngs_joint_calling.smk"

# Long-read alignment rules are shared by TGS SV and small-variant workflows.
if (
    config["workflows"]["enable_tgs_sv_detection"]
    or config["workflows"].get("enable_tgs_variant_calling", False)
):
    include: "rules/tgs_sv.smk"

if config["workflows"]["enable_tgs_sv_detection"]:
    include: "rules/tgs_sv_part2.smk"

if config["workflows"].get("enable_tgs_variant_calling", False):
    include: "rules/tgs_variant_calling.smk"
    include: "rules/tgs_deepvariant.smk"

if config["workflows"].get("enable_tgs_variant_annotation", False):
    tgs_annotation_tool = config.get("tgs_variant_annotation", {}).get("tool", "snpeff").strip().lower()
    if tgs_annotation_tool == "snpeff":
        include: "rules/tgs_snpeff.smk"
    elif tgs_annotation_tool == "annovar":
        include: "rules/tgs_annovar.smk"

if config["workflows"].get("enable_tgs_joint_calling", False):
    include: "rules/tgs_joint_calling.smk"
    include: "rules/tgs_joint_calling_part2.smk"
    if config["workflows"].get("enable_tgs_phasing", False) and tgs_phasing_cfg.get("enabled", False):
        include: "rules/tgs_phasing.smk"

if config["workflows"].get("enable_ngs_sv_detection", False):
    include: "rules/ngs_sv.smk"

# SV annotation with AnnotSV
if config["workflows"].get("enable_tgs_sv_detection", False) or config["workflows"].get("enable_ngs_sv_detection", False):
    include: "rules/annotsv.smk"

# Cross-platform (NGS + TGS) SV integration for matched samples (opt-in).
# Requires both NGS and TGS SV detection so their per-sample SV VCFs exist.
if config["workflows"].get("enable_cross_platform_sv", False):
    include: "rules/cross_platform_sv.smk"

if config["workflows"].get("enable_multiqc", False):
    include: "rules/multiqc.smk"

# Tool deployment check (Conda or containers) without input data.
if config["workflows"].get("enable_env_check", False):
    include: "rules/env_check.smk"

# Include output definitions
include: "workflow/outputs.smk"

# ---------- Target rule to define all output files ----------
rule all:
    input:
        # Conditionally include outputs based on workflow settings
        get_basic_workflow_outputs() if config["workflows"]["enable_basic_ngs_workflows"] else [],
        get_joint_calling_outputs() if config["workflows"]["enable_ngs_joint_calling"] else [],
        get_sv_detection_outputs() if config["workflows"]["enable_tgs_sv_detection"] else [],
        get_ngs_sv_outputs() if config["workflows"]["enable_ngs_sv_detection"] else [],
        get_tgs_variant_calling_outputs_wrapper() if config["workflows"].get("enable_tgs_variant_calling", False) else [],
        get_tgs_annotation_outputs_wrapper() if config["workflows"].get("enable_tgs_variant_annotation", False) else [],
        get_tgs_joint_calling_outputs() if config["workflows"].get("enable_tgs_joint_calling", False) else [],
        get_cross_platform_sv_outputs() if config["workflows"].get("enable_cross_platform_sv", False) else [],
        get_multiqc_outputs() if config["workflows"].get("enable_multiqc", False) else [],
        f"{ENV_CHECK_DIR}/report.tsv" if config["workflows"].get("enable_env_check", False) else [],
        get_hitsv_outputs()
       
       
