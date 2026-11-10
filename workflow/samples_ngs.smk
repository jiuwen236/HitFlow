# NGS sample discovery and management
# This file handles loading and processing of NGS (short-read) samples

import os
import csv
import glob
import re

# CSV file path for NGS samples
SAMPLES_CSV = resolve_input_file_path(config.get("samples_csv", "samples.csv"))


def _is_pattern_path(path_value):
    """Return True if the path contains glob or brace-style patterns."""
    if not path_value:
        return False
    return glob.has_magic(path_value) or "{" in path_value or "}" in path_value


def _expand_fastq_inputs(path_value):
    """Expand FASTQ entries supporting glob/brace patterns and base directories."""
    if not path_value:
        return []

    raw_value = path_value.strip()
    if not raw_value:
        return []

    candidate_patterns = []

    if os.path.isabs(raw_value):
        candidate_patterns.append(raw_value)
    else:
        candidate_patterns.append(raw_value)
        for base_dir in BASE_DIRS:
            base_dir = (base_dir or "").strip()
            if base_dir:
                candidate_patterns.append(os.path.join(base_dir, raw_value))

    expanded_paths = []
    for pattern in candidate_patterns:
        if _is_pattern_path(pattern):
            expanded_paths.extend(_expand_path_pattern(pattern))
        else:
            expanded_paths.append(resolve_input_file_path(pattern, BASE_DIRS))

    seen = set()
    normalized = []
    for path in expanded_paths:
        if not path:
            continue
        norm_path = os.path.normpath(path) if os.path.isabs(path) else path
        if norm_path not in seen:
            normalized.append(norm_path)
            seen.add(norm_path)
    return normalized


def _load_samples_csv(csv_path):
    """Load samples from CSV file. 
    
    Supports multiple rows for the same sample - FASTQ files will be merged.
    Each row adds a pair of FASTQ files to the sample's file list.
    
    Supports optional 'bam' column. If provided, alignment steps are skipped.
    When 'bam' is specified, 'fq1' and 'fq2' columns are optional.
    """
    samples = {}
    if os.path.exists(csv_path):
        with open(csv_path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sample = row["sample"].strip()
                fq1_entries = _expand_fastq_inputs(row.get("fq1"))
                fq2_entries = _expand_fastq_inputs(row.get("fq2"))
                bam = row["bam"].strip() if row.get("bam") else None

                if bam:
                    bam = resolve_input_file_path(bam, BASE_DIRS)

                # Initialize sample if not exists
                if sample not in samples:
                    samples[sample] = {
                        "fq1": [],
                        "fq2": [],
                        "bam": None
                    }
                
                # Append FASTQ files to lists (supporting multiple rows per sample)
                if fq1_entries:
                    samples[sample]["fq1"].extend(fq1_entries)
                if fq2_entries:
                    samples[sample]["fq2"].extend(fq2_entries)
                
                # Set BAM path (only first row's bam is used if multiple rows for same sample)
                if bam and samples[sample]["bam"] is None:
                    samples[sample]["bam"] = bam
    return samples


def _scan_fastq_dirs(dirs):
    """Scan directories for FASTQ pairs and create sample dictionary."""
    samples = {}
    fastq_regex = re.compile(
        r"(.+?)[._-](?:R|read)?([12])(?:[._-][0-9A-Za-z]+)?\.f(?:ast)?q(?:\.gz)?$",
        re.IGNORECASE,
    )
    
    fq_pairs = {}  # {base_name: {1: path1, 2: path2}}
    
    # Filter out empty strings to avoid scanning current directory
    dirs = [d for d in dirs if d and d.strip()]
    
    for d in dirs:
        for fq_path in glob.glob(os.path.join(d, "**", "*.f*q*"), recursive=True):
            fname = os.path.basename(fq_path)
            m = fastq_regex.match(fname)
            if not m:
                continue
            base_name = m.group(1)
            read_num = int(m.group(2))
            
            if base_name not in fq_pairs:
                fq_pairs[base_name] = {}
            fq_pairs[base_name][read_num] = fq_path
    
    # Convert to samples format (lists for consistency with CSV loader)
    for base_name, reads in fq_pairs.items():
        if 1 in reads and 2 in reads:  # Only complete pairs
            samples[base_name] = {
                "fq1": [reads[1]],
                "fq2": [reads[2]],
                "bam": None
            }
    
    return samples


def _merge_samples_dict(*dicts):
    """Merge multiple sample dictionaries. 
    
    For the same sample appearing in multiple dictionaries, merge their FASTQ lists.
    BAM paths are taken from the first dictionary that provides one.
    """
    merged = {}
    for d in dicts:
        for sid, info in d.items():
            if sid not in merged:
                merged[sid] = {"fq1": [], "fq2": [], "bam": None}
            
            # Merge FASTQ lists
            if isinstance(info.get("fq1"), list):
                merged[sid]["fq1"].extend(info["fq1"])
            elif info.get("fq1"):
                merged[sid]["fq1"].append(info["fq1"])
            
            if isinstance(info.get("fq2"), list):
                merged[sid]["fq2"].extend(info["fq2"])
            elif info.get("fq2"):
                merged[sid]["fq2"].append(info["fq2"])
            
            # BAM path: first one wins
            if info.get("bam") and merged[sid]["bam"] is None:
                merged[sid]["bam"] = info["bam"]
    
    return merged


# Load samples from all sources
SAMPLES_FROM_CSV = _load_samples_csv(SAMPLES_CSV)
SAMPLES_FROM_SCAN = _scan_fastq_dirs(SCAN_DIRS)

# Create the final samples dictionary
ALL_SAMPLES = _merge_samples_dict(SAMPLES_FROM_CSV, SAMPLES_FROM_SCAN)

# Write NGS samples log
_scan_dir_label = ",".join([d for d in SCAN_DIRS if d and d.strip()])
_write_samples_log_csv(SAMPLES_FROM_CSV, f"{LOGS_DIR}/ngs_samples_from_csv.csv", "NGS", source_path=SAMPLES_CSV)
_write_samples_log_csv(SAMPLES_FROM_SCAN, f"{LOGS_DIR}/ngs_samples_from_scan.csv", "NGS", source_path=_scan_dir_label)
_write_samples_log_csv(
    ALL_SAMPLES,
    f"{LOGS_DIR}/ngs_samples_merged.csv",
    "NGS",
    source_path=f"csv:{SAMPLES_CSV or 'NA'};scan:{_scan_dir_label or 'NA'}"
)

# ---------- NGS sample accessor functions ----------
def get_fastq_path(wildcards, read_num):
    """Get FASTQ path(s) for a sample. Returns list of paths."""
    sample_id = wildcards.sample_id
    
    if sample_id not in ALL_SAMPLES:
        raise ValueError(f"Sample '{sample_id}' not found in samples")
    
    fq_key = f"fq{read_num}"
    paths = ALL_SAMPLES[sample_id].get(fq_key)
    
    if not paths:
        raise ValueError(f"FASTQ R{read_num} not found for sample '{sample_id}'")
    
    # Return list (already in list format after merge)
    return paths if isinstance(paths, list) else [paths]


def get_input_fastq_for_alignment(wildcards, read_num):
    """
    Get FASTQ path for alignment - returns fastp-cleaned FASTQ if fastp is enabled,
    otherwise returns raw FASTQ path. This ensures backward compatibility.
    """
    sample_id = wildcards.sample_id
    
    # Check if fastp is enabled
    if config["fastp"]["enabled"]:
        # Return fastp-cleaned FASTQ from preprocessing directory
        return f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastp/{sample_id}.clean.R{read_num}.fq.gz"
    else:
        # Return raw FASTQ (backward compatible)
        return get_fastq_path(wildcards, str(read_num))


def sample_has_bam(sample_id):
    """Check if sample has a pre-aligned BAM file specified."""
    if sample_id not in ALL_SAMPLES:
        return False
    return bool(ALL_SAMPLES[sample_id].get("bam"))


def sample_has_fastq(sample_id):
    """Check if sample has FASTQ files specified."""
    if sample_id not in ALL_SAMPLES:
        return False
    info = ALL_SAMPLES[sample_id]
    return bool(info.get("fq1") and info.get("fq2"))


def get_sample_bam_path(sample_id):
    """Get the user-provided BAM path for a sample."""
    if sample_id not in ALL_SAMPLES:
        raise ValueError(f"Sample '{sample_id}' not found in samples")
    bam = ALL_SAMPLES[sample_id].get("bam")
    if not bam:
        raise ValueError(f"BAM path not found for sample '{sample_id}'")
    return bam


def get_all_sample_ids():
    """Get list of all sample IDs that have FASTQ files or BAM files."""
    return [sid for sid, info in ALL_SAMPLES.items() 
            if (info.get("fq1") and info.get("fq2")) or info.get("bam")]


def get_sample_ids_with_fastq():
    """Get list of sample IDs that have FASTQ files (require alignment)."""
    return [sid for sid, info in ALL_SAMPLES.items() 
            if info.get("fq1") and info.get("fq2") and not info.get("bam")]


def get_sample_ids_with_bam():
    """Get list of sample IDs that have BAM files (skip alignment)."""
    return [sid for sid, info in ALL_SAMPLES.items() if info.get("bam")]


# Regex helper for rules limited to BAM-provided samples
_BAM_SAMPLE_IDS = get_sample_ids_with_bam()
SYMLINK_BAM_WILDCARD_REGEX = (
    "|".join([re.escape(sid) for sid in _BAM_SAMPLE_IDS])
    if _BAM_SAMPLE_IDS else r"(?!x)x"
)


def get_bam_for_variant_calling(wildcards):
    """
    Get BAM file for variant calling.
    Returns user-provided BAM if available, otherwise returns pipeline-generated recal BAM.
    """
    sample_id = wildcards.sample_id
    if sample_has_bam(sample_id):
        return get_sample_bam_path(sample_id)
    else:
        return f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.dup.sorted.recal.bam"


def get_bai_for_variant_calling(wildcards):
    """
    Get BAM index file for variant calling.
    Returns user-provided BAI if available, otherwise returns pipeline-generated recal BAI.
    Assumes BAI is at <bam>.bai or <bam_without_extension>.bai.
    """
    sample_id = wildcards.sample_id
    if sample_has_bam(sample_id):
        bam_path = get_sample_bam_path(sample_id)
        # Try common BAI naming conventions
        bai_path = bam_path + ".bai"
        if not os.path.exists(bai_path):
            # Try replacing .bam with .bai
            bai_path = bam_path.replace(".bam", ".bai")
        return bai_path
    else:
        return f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.dup.sorted.recal.bai"


def get_bam_for_bqsr(wildcards):
    """
    Get BAM file for BQSR.
    Returns user-provided BAM if available, otherwise returns pipeline-generated dup.sorted BAM.
    """
    sample_id = wildcards.sample_id
    if sample_has_bam(sample_id):
        return get_sample_bam_path(sample_id)
    else:
        return f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.dup.sorted.bam"


def get_bai_for_bqsr(wildcards):
    """
    Get BAM index for BQSR.
    """
    sample_id = wildcards.sample_id
    if sample_has_bam(sample_id):
        bam_path = get_sample_bam_path(sample_id)
        bai_path = bam_path + ".bai"
        if not os.path.exists(bai_path):
            bai_path = bam_path.replace(".bam", ".bai")
        return bai_path
    else:
        return f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.dup.sorted.bai"
