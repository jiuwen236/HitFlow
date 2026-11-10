# Helper functions for file path resolution and sample processing
# This file contains utility functions used across the workflow

import os
import csv
import glob
import re
import subprocess

from scripts.snakemake.resources import (
    memory_mb,
    retry_memory_mb,
    samtools_sort_memory_mb,
    sambamba_sort_memory_mb,
    deepvariant_memory_mb,
)


def _retry_mem_mb(base_mb, attempt=1):
    return retry_memory_mb(
        base_mb, attempt, config.get("cluster", {}).get("retry_max_mem_mb", 128 * 1024)
    )


# ---------- File path resolution function ----------
def resolve_input_file_path(file_path, base_dirs=None):
    """
    Resolve input file path by checking if it exists, and if not, try with base directories.
    
    Args:
        file_path (str): The file path to resolve
        base_dirs (list): List of base directories to try if file doesn't exist
    
    Returns:
        str: The resolved file path
    """
    if not file_path:
        # Allow optional/disabled inputs (e.g., empty known-sites)
        return ""

    if base_dirs is None:
        base_dirs = config["input_data"]["base_dir"]
    
    # If it's an absolute path, return as is
    if os.path.isabs(file_path):
        return file_path
    
    # If the relative path exists, return as is
    if os.path.exists(file_path):
        return file_path
    
    # Try with each base directory
    for base_dir in base_dirs:
        full_path = os.path.join(base_dir, file_path)
        if os.path.exists(full_path):
            return full_path
    
    # If file not found, return original path (let snakemake handle the error)
    return file_path


# ---------- CSV logging function ----------
def _write_samples_log_csv(samples_dict, log_path, sample_type="NGS", source_path=""):
    """Write samples dictionary to CSV log file for debugging.
    
    Args:
        samples_dict: Dictionary of samples
        log_path: Path to output CSV file
        sample_type: "NGS" or "TGS" for different formats
    """
    from datetime import datetime
    
    # Create logs directory if not exists
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    
    with open(log_path, 'w', newline='') as f:
        if sample_type == "NGS":
            writer = csv.writer(f)
            writer.writerow(['# Generated at', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
            if source_path:
                writer.writerow(['# Source', source_path])
            writer.writerow(['sample', 'fq1_count', 'fq2_count', 'fq1_files', 'fq2_files'])
            
            for sample, info in sorted(samples_dict.items()):
                fq1_list = info.get('fq1', [])
                fq2_list = info.get('fq2', [])
                fq1_count = len(fq1_list) if isinstance(fq1_list, list) else (1 if fq1_list else 0)
                fq2_count = len(fq2_list) if isinstance(fq2_list, list) else (1 if fq2_list else 0)
                fq1_str = ';'.join(fq1_list) if isinstance(fq1_list, list) else (fq1_list or '')
                fq2_str = ';'.join(fq2_list) if isinstance(fq2_list, list) else (fq2_list or '')
                
                writer.writerow([sample, fq1_count, fq2_count, fq1_str, fq2_str])
        
        else:  # TGS
            writer = csv.writer(f)
            writer.writerow(['# Generated at', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
            if source_path:
                writer.writerow(['# Source', source_path])
            writer.writerow(['sample', 'file_count', 'fastq_count', 'bam_file', 'file_paths', 'fastq_files'])
            
            for sample, info in sorted(samples_dict.items()):
                file_paths = info.get('file_path', [])
                fastq_files = info.get('fastq_files', [])
                bam_file = info.get('bam_file') or info.get('unaligned_bam') or ''
                
                file_count = len(file_paths) if isinstance(file_paths, list) else (1 if file_paths else 0)
                fastq_count = len(fastq_files) if isinstance(fastq_files, list) else 0
                
                file_paths_str = ';'.join(file_paths) if isinstance(file_paths, list) else (file_paths or '')
                fastq_files_str = ';'.join(fastq_files) if isinstance(fastq_files, list) else ''
                
                writer.writerow([sample, file_count, fastq_count, bam_file, file_paths_str, fastq_files_str])


# ---------- Path pattern expansion functions for TGS samples ----------
def _expand_path_pattern(path_pattern):
    """
    Expand path patterns containing {a,b,c} brace expansion and * wildcards.
    Returns a list of all matching files.
    
    Examples:
        /path/{1,2}/*/file.txt -> [/path/1/subdir1/file.txt, /path/1/subdir2/file.txt, /path/2/subdir1/file.txt, ...]
        /path/to/*.fastq.gz -> [/path/to/file1.fastq.gz, /path/to/file2.fastq.gz, ...]
    """
    # Use bash to expand the pattern (handles both brace expansion and wildcards)
    try:
        result = subprocess.run(
            ['bash', '-c', f'echo {path_pattern}'],
            capture_output=True, text=True, timeout=10
        )
        expanded = result.stdout.strip()
        
        # If expansion produced multiple paths, split them
        if expanded and expanded != path_pattern:
            potential_paths = expanded.split()
            # Filter out paths that don't exist or are the original pattern (no match)
            matched_files = []
            for p in potential_paths:
                # Use glob to get all files matching the pattern
                matched = glob.glob(p, recursive=False)
                matched_files.extend(matched)
            
            if matched_files:
                return sorted(set(matched_files))  # Remove duplicates and sort
        
        # Fallback: try direct glob expansion
        matched = glob.glob(path_pattern, recursive=False)
        if matched:
            return sorted(matched)
        
        # If no matches found, return the original pattern
        return [path_pattern]
    
    except Exception as e:
        print(f"Warning: Could not expand pattern '{path_pattern}': {e}")
        return [path_pattern]


# ---------- TGS variant calling tool resolution ----------
def get_tgs_variant_calling_tool():
    """
    Resolve TGS variant calling tool based on configuration.
    
    If tool is "Default" (case-insensitive), automatically select based on platform:
      - ONT -> Clair3 (recommended for Oxford Nanopore)
      - PACBIO_HIFI -> DeepVariant (recommended for PacBio HiFi)
    
    Returns:
        str: The resolved tool name in lowercase ("clair3" or "deepvariant")
    """
    selected_tool = config.get("tgs_variant_calling", {}).get("tool", "Default").strip()
    
    if selected_tool.lower() == "default":
        platform = config.get("tgs_config", {}).get("platform", "ONT").strip().upper()
        if platform == "PACBIO_HIFI":
            return "deepvariant"
        else:
            # ONT or any other platform defaults to Clair3
            return "clair3"
    else:
        return selected_tool.lower()


def _find_fastq_files_in_pattern(path_pattern):
    """
    Find all FASTQ files matching the given path pattern.
    Supports brace expansion {a,b,c} and wildcards *.
    
    If path_pattern is a directory or directory pattern, recursively search for FASTQ files.
    If path_pattern points to specific files, return those files.
    """
    expanded_paths = _expand_path_pattern(path_pattern)
    
    all_fastq_files = []
    fastq_extensions = ('.fq', '.fastq', '.fq.gz', '.fastq.gz')
    
    for path in expanded_paths:
        if os.path.isfile(path):
            # If it's a file, check if it's a FASTQ file
            if path.lower().endswith(fastq_extensions):
                all_fastq_files.append(path)
        elif os.path.isdir(path):
            # If it's a directory, recursively find all FASTQ files
            for root, dirs, files in os.walk(path):
                for f in files:
                    if f.lower().endswith(fastq_extensions):
                        all_fastq_files.append(os.path.join(root, f))
    
    return sorted(set(all_fastq_files))  # Remove duplicates and sort
