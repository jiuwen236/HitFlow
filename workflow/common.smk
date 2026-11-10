# Common configuration and directory definitions
# This file contains shared configuration loaded from config.yaml

import os

# Load configuration (support duplicate keys with deep-merge)
from scripts.snakemake.merge_yaml_duplicates import load_yaml_with_merged_duplicates, deep_merge_dicts

# Capture config provided via the Snakemake CLI (--config / --configfile) before
# we re-read the YAML ourselves to preserve duplicate top-level keys.
try:
    _cli_cfg = config if isinstance(config, dict) else {}
except NameError:
    _cli_cfg = {}

# Resolve which config file to (re-)read. Honor --configfile when given so the
# workflow can run from any working directory (e.g. a per-test subfolder run via
# `snakemake --directory <folder>`); otherwise fall back to "config.yaml" in the
# current working directory. All paths are resolved relative to the working dir.
_config_candidates = []
try:
    _config_candidates.extend(str(p) for p in workflow.configfiles)
except (NameError, AttributeError):
    pass
_config_candidates.append("config.yaml")

_config_path = next((p for p in _config_candidates if p and os.path.exists(p)), None)
if _config_path is None:
    raise FileNotFoundError(
        "HitFlow could not locate a config file. Pass one explicitly with "
        "'snakemake --configfile <path>', or run from a directory that contains "
        "'config.yaml' (e.g. 'snakemake --directory <test-folder>'). "
        f"Looked for: {', '.join(_config_candidates)}; cwd: {os.getcwd()}."
    )

_yaml_cfg = load_yaml_with_merged_duplicates(_config_path)
config = deep_merge_dicts(_yaml_cfg, _cli_cfg)

# Bundled helper scripts and binaries live next to the Snakefile, so runs work
# from any working directory. Relative values of the script-path settings
# below are first looked up in the working directory (previous behaviour),
# then in the workflow directory.
WORKFLOW_DIR = os.path.abspath(workflow.basedir)


def workflow_path(*parts):
    return os.path.join(WORKFLOW_DIR, *parts)


def _resolve_workflow_path(value):
    value = str(value or "").strip()
    if not value or os.path.isabs(value) or os.path.exists(value):
        return value
    candidate = workflow_path(value)
    return candidate if os.path.exists(candidate) else value


for _key, _default in (
    ("ngs_scripts", "scripts/ngs_scripts"),
    ("sv_pipeline_scripts", "scripts/sv_pipeline_scripts"),
    ("configure_path", "scripts/ngs_scripts/configure.v1"),
):
    config[_key] = _resolve_workflow_path(config.get(_key) or _default)
_cnvpytor_cfg = config.get("ngs_sv", {}).get("cnvpytor", {})
if isinstance(_cnvpytor_cfg, dict) and _cnvpytor_cfg.get("to_vcf_script"):
    _cnvpytor_cfg["to_vcf_script"] = _resolve_workflow_path(_cnvpytor_cfg["to_vcf_script"])
IMPORT_SHAPEIT5 = workflow_path("scripts/pipeline/import_shapeit5.sh")

# Define common directories using new hierarchical structure
NGS_PREPROCESSING_DIR = config["output"]["ngs"]["preprocessing"]
NGS_ALIGNMENT_DIR = config["output"]["ngs"]["alignment"]
NGS_VARIANTS_SINGLE_DIR = config["output"]["ngs"]["variants"]["single_sample"]
NGS_VARIANTS_COHORT_DIR = config["output"]["ngs"]["variants"]["cohort"]
NGS_SV_SINGLE_DIR = config["output"]["ngs"]["sv"]["single_sample"]
NGS_SV_COHORT_DIR = config["output"]["ngs"]["sv"]["cohort"]
NGS_QC_REPORTS_DIR = config["output"]["ngs"]["qc_reports"]

TGS_PREPROCESSING_DIR = config["output"]["tgs"]["preprocessing"]
TGS_VARIANTS_SINGLE_DIR = config["output"]["tgs"]["variants"]["single_sample"]
TGS_VARIANTS_COHORT_DIR = config["output"]["tgs"]["variants"]["cohort"]
TGS_SV_SINGLE_DIR = config["output"]["tgs"]["sv"]["single_sample"]
TGS_SV_COHORT_DIR = config["output"]["tgs"]["sv"]["cohort"]
TGS_QC_REPORTS_DIR = config["output"]["tgs"]["qc_reports"]

LOGS_DIR = config["output"]["logs"]
TMP_DIR = config["output"]["tmp"]


# Prepare cluster logs before the scheduler submits its first job.
onstart:
    if not workflow.local_exec:
        os.makedirs("logs_cluster", exist_ok=True)

# Legacy variables for backward compatibility (deprecated - will be removed)
FINAL_DIR = NGS_VARIANTS_SINGLE_DIR  # Deprecated: use NGS_VARIANTS_SINGLE_DIR
ALIGNMENT_DIR = NGS_ALIGNMENT_DIR  # Deprecated: use NGS_ALIGNMENT_DIR
JOINT_CALLING_DIR = NGS_VARIANTS_COHORT_DIR  # Deprecated: use NGS_VARIANTS_COHORT_DIR

# Chromosome list for expansion
CHROMS = config["chromosomes"]

# Base directories for input data
BASE_DIRS = config["input_data"]["base_dir"]
SCAN_DIRS = config.get("scan_ngs_dirs", [])
SCAN_TGS_DIRS = config.get("scan_tgs_dirs", [])
