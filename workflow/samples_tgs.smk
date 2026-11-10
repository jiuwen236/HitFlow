# TGS/ONT sample discovery and management
# This file handles loading and processing of TGS (long-read) samples

import os
import csv
import glob
import re
import subprocess

# CSV file paths for TGS samples
SAMPLES_ONT_CSV = resolve_input_file_path(config.get("samples_ont_csv", "samples_ont.csv"))


def _load_ont_samples_csv(csv_path):
    """Load TGS samples from CSV file. 
    
    Supports multiple rows for the same sample - files will be merged.
    
    The file_path can be:
    - A single FASTQ/BAM file
    - A directory containing FASTQ files
    - A path pattern with brace expansion: /path/{1,2,3}/files/
    - A path pattern with wildcards: /path/*/barcode01/
    - A combination: /path/{1,2}/*/fastq_pass/barcode03/
    
    For multiple FASTQ files, they will be automatically merged during alignment.
    """
    samples = {}
    if os.path.exists(csv_path):
        with open(csv_path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sample = row["sample"].strip()
                file_path = row["file_path"].strip() if row.get("file_path") else None
                
                if not file_path:
                    continue
                
                # Convert relative paths to absolute using resolve_input_file_path
                file_path = resolve_input_file_path(file_path, BASE_DIRS)
                
                # Initialize sample if not exists
                if sample not in samples:
                    samples[sample] = {
                        "file_path": [],
                        "fastq_files": [],
                        "fastq_file": None,
                        "bam_file": None,
                        "unaligned_bam": None
                    }
                
                # Check if this is a path pattern (contains {,} or ends with /)
                is_pattern = '{' in file_path or file_path.endswith('/')
                
                if is_pattern:
                    # This is a pattern - find all matching FASTQ files
                    fastq_files = _find_fastq_files_in_pattern(file_path)
                    
                    if not fastq_files:
                        print(f"Warning: No FASTQ files found for pattern '{file_path}' (sample: {sample})")
                        continue
                    
                    # Append to existing file lists
                    samples[sample]["file_path"].append(file_path)
                    samples[sample]["fastq_files"].extend(fastq_files)
                else:
                    # Single file path
                    is_bam = file_path.endswith('.bam')
                    is_aligned_bam = False
                    
                    aligned_hint = (row.get("aligned") or "").strip().lower()
                    if aligned_hint not in ("", "true", "false"):
                        raise ValueError("samples_tgs.csv aligned must be true, false or empty")
                    if is_bam and aligned_hint:
                        # Explicit metadata avoids rescanning large BAMs on every worker.
                        # Downstream callers still validate the header/reference/index.
                        is_aligned_bam = aligned_hint == "true"
                    elif is_bam and os.path.exists(file_path):
                        # Check if BAM contains aligned reads
                        try:
                            result = subprocess.run(
                                ['samtools', 'view', '-c', '-F', '4', file_path],
                                capture_output=True, text=True, timeout=30
                            )
                            aligned_count = int(result.stdout.strip())
                            is_aligned_bam = aligned_count > 0
                        except Exception as e:
                            print(f"Warning: Could not check alignment status for {file_path}: {e}")
                            is_aligned_bam = False
                    
                    # Append to file lists
                    samples[sample]["file_path"].append(file_path)
                    
                    if not is_bam:
                        samples[sample]["fastq_files"].append(file_path)
                        if not samples[sample]["fastq_file"]:
                            samples[sample]["fastq_file"] = file_path
                    elif is_aligned_bam:
                        if not samples[sample]["bam_file"]:
                            samples[sample]["bam_file"] = file_path
                    else:
                        if not samples[sample]["unaligned_bam"]:
                            samples[sample]["unaligned_bam"] = file_path
    
    return samples


def _scan_tgs_dirs(dirs):
    """Scan directories for TGS FASTQ and BAM files and create sample dictionary.
    
    Automatically detects:
    - FASTQ files (.fq, .fastq, .fq.gz, .fastq.gz)
    - BAM files (.bam) and checks alignment status
    
    Sample names are derived from base filenames (without extensions).
    """
    samples = {}
    
    # Filter out empty strings to avoid scanning current directory
    dirs = [d for d in dirs if d and d.strip()]
    
    # Regex patterns for TGS files (single-end, no R1/R2 distinction)
    # Match: sample.fq.gz, sample.fastq, sample_nanopore.fq.gz, etc.
    fastq_regex = re.compile(r"^(.+?)\.f(?:ast)?q(?:\.gz)?$", re.IGNORECASE)
    bam_regex = re.compile(r"^(.+?)\.bam$", re.IGNORECASE)
    
    for d in dirs:
        if not os.path.exists(d):
            print(f"Warning: TGS scan directory does not exist: {d}")
            continue
            
        # Scan for FASTQ files
        for fq_path in glob.glob(os.path.join(d, "**", "*.f*q*"), recursive=True):
            fname = os.path.basename(fq_path)
            m = fastq_regex.match(fname)
            if not m:
                continue
            
            base_name = m.group(1)
            
            # Skip if already found (first occurrence wins)
            if base_name in samples:
                continue
            
            samples[base_name] = {
                "file_path": [fq_path],  # Store as list for consistency
                "fastq_files": [fq_path],  # Store as list for consistency
                "fastq_file": fq_path,
                "bam_file": None,
                "unaligned_bam": None
            }
        
        # Scan for BAM files
        for bam_path in glob.glob(os.path.join(d, "**", "*.bam"), recursive=True):
            fname = os.path.basename(bam_path)
            m = bam_regex.match(fname)
            if not m:
                continue
            
            base_name = m.group(1)
            
            # Skip if already found (FASTQ takes priority)
            if base_name in samples:
                continue
            
            # Check if BAM is aligned
            is_aligned_bam = False
            if os.path.exists(bam_path):
                try:
                    result = subprocess.run(
                        ['samtools', 'view', '-c', '-F', '4', bam_path],
                        capture_output=True, text=True, timeout=30
                    )
                    aligned_count = int(result.stdout.strip())
                    is_aligned_bam = aligned_count > 0
                except Exception as e:
                    print(f"Warning: Could not check alignment status for {bam_path}: {e}")
                    is_aligned_bam = False
            
            samples[base_name] = {
                "file_path": [bam_path],  # Store as list for consistency
                "fastq_files": [],  # No FASTQ files for BAM
                "fastq_file": None,
                "bam_file": bam_path if is_aligned_bam else None,
                "unaligned_bam": bam_path if not is_aligned_bam else None
            }
    
    return samples


def _merge_tgs_samples_dict(*dicts):
    """Merge multiple TGS sample dictionaries.
    
    For the same sample appearing in multiple dictionaries, merge their file lists.
    """
    merged = {}
    for d in dicts:
        for sid, info in d.items():
            if sid not in merged:
                merged[sid] = {
                    "file_path": [],
                    "fastq_files": [],
                    "fastq_file": None,
                    "bam_file": None,
                    "unaligned_bam": None
                }
            
            # Merge file_path lists
            file_paths = info.get("file_path", [])
            if isinstance(file_paths, list):
                merged[sid]["file_path"].extend(file_paths)
            elif file_paths:
                merged[sid]["file_path"].append(file_paths)
            
            # Merge fastq_files lists
            fastq_files = info.get("fastq_files", [])
            if isinstance(fastq_files, list):
                merged[sid]["fastq_files"].extend(fastq_files)
            elif fastq_files:
                merged[sid]["fastq_files"].append(fastq_files)
            
            # Keep first non-None values for single file fields
            if not merged[sid]["fastq_file"] and info.get("fastq_file"):
                merged[sid]["fastq_file"] = info["fastq_file"]
            if not merged[sid]["bam_file"] and info.get("bam_file"):
                merged[sid]["bam_file"] = info["bam_file"]
            if not merged[sid]["unaligned_bam"] and info.get("unaligned_bam"):
                merged[sid]["unaligned_bam"] = info["unaligned_bam"]
    
    return merged


# Load ONT/TGS samples from CSV and directory scanning
ONT_SAMPLES_FROM_CSV = _load_ont_samples_csv(SAMPLES_ONT_CSV)
ONT_SAMPLES_FROM_SCAN = _scan_tgs_dirs(SCAN_TGS_DIRS)

# Merge TGS samples using dedicated merge function
ONT_SAMPLES = _merge_tgs_samples_dict(ONT_SAMPLES_FROM_CSV, ONT_SAMPLES_FROM_SCAN)

# Write TGS samples log
_write_samples_log_csv(ONT_SAMPLES_FROM_CSV, f"{LOGS_DIR}/tgs_samples_from_csv.csv", "TGS", source_path=SAMPLES_ONT_CSV)
_write_samples_log_csv(ONT_SAMPLES_FROM_SCAN, f"{LOGS_DIR}/tgs_samples_from_scan.csv", "TGS")
_write_samples_log_csv(ONT_SAMPLES, f"{LOGS_DIR}/tgs_samples_merged.csv", "TGS")


# ---------- TGS sample accessor functions ----------
def get_ont_fastq_path(wildcards):
    """Get TGS FASTQ path(s) for a sample.
    
    Returns:
        - For single FASTQ: returns the file path
        - For multiple FASTQs: returns a list of file paths (will be merged in alignment rule)
    """
    sample = wildcards.sample
    
    if sample not in ONT_SAMPLES:
        return f"MISSING_SAMPLE_{sample}"
    
    # Check if we have multiple FASTQ files
    fastq_files = ONT_SAMPLES[sample].get("fastq_files", [])
    
    if not fastq_files:
        # Try the old single file path for backward compatibility
        path = ONT_SAMPLES[sample].get("fastq_file")
        if not path:
            return f"MISSING_FASTQ_FOR_{sample}"
        return path
    
    # Return list of files if multiple, single file if only one
    if len(fastq_files) == 1:
        return fastq_files[0]
    else:
        return fastq_files


def get_ont_unaligned_bam_path(wildcards):
    """Get TGS unaligned BAM path for a sample (for conversion to FASTQ)."""
    sample = wildcards.sample
    
    if sample not in ONT_SAMPLES:
        return f"MISSING_SAMPLE_{sample}"
    
    path = ONT_SAMPLES[sample].get("unaligned_bam")
    
    if not path:
        return f"MISSING_UNALIGNED_BAM_FOR_{sample}"
    
    return path


def get_ont_bam_path(wildcards):
    """Get TGS BAM path for a sample (pre-aligned BAM)."""
    sample = wildcards.sample
    
    if sample not in ONT_SAMPLES:
        return f"MISSING_SAMPLE_{sample}"
    
    path = ONT_SAMPLES[sample].get("bam_file")
    
    if not path:
        return f"MISSING_BAM_FOR_{sample}"
    
    return path


def has_ont_bam(wildcards):
    """Check if sample has a pre-aligned BAM file."""
    sample = wildcards.sample
    return sample in ONT_SAMPLES and ONT_SAMPLES[sample].get("bam_file") is not None


def has_ont_unaligned_bam(wildcards):
    """Check if sample has an unaligned BAM file."""
    sample = wildcards.sample
    return sample in ONT_SAMPLES and ONT_SAMPLES[sample].get("unaligned_bam") is not None


def get_ont_samples():
    """Get list of all TGS sample IDs."""
    return [sid for sid, info in ONT_SAMPLES.items() if info.get("file_path")]
