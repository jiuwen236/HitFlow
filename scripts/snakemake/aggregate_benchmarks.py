#!/usr/bin/env python3
"""
Aggregate Benchmark Statistics from Snakemake Workflow

This script aggregates benchmark data from all Snakemake benchmark files,
providing comprehensive statistics on runtime, memory usage, and I/O.

Usage:
    python scripts/aggregate_benchmarks.py [benchmark_dir] [--format {table|json|csv}]
    
Examples:
    # Basic usage with default directory
    python scripts/aggregate_benchmarks.py
    
    # Specify benchmark directory
    python scripts/aggregate_benchmarks.py logs/benchmarks/
    
    # Output as JSON
    python scripts/aggregate_benchmarks.py --format json
    
    # Output as CSV
    python scripts/aggregate_benchmarks.py --format csv > benchmarks_summary.csv
"""

import os
import sys
import glob
import argparse
import json
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple
import statistics


def parse_benchmark_file(filepath: str) -> Dict:
    """Parse a single Snakemake benchmark file (TSV format)."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            return None
        
        # Parse header and data
        header = lines[0].strip().split('\t')
        data = lines[1].strip().split('\t')
        
        if len(header) != len(data):
            return None
        
        # Create dictionary
        result = {}
        for key, value in zip(header, data):
            try:
                # Try to convert to float
                result[key] = float(value)
            except ValueError:
                # Keep as string if not numeric
                result[key] = value
        
        # Add filepath for reference
        result['_filepath'] = filepath
        result['_filename'] = os.path.basename(filepath)
        
        return result
    except Exception as e:
        print(f"Warning: Failed to parse {filepath}: {e}", file=sys.stderr)
        return None


def get_rule_name_from_path(filepath: str) -> str:
    """Extract rule name from benchmark file path."""
    # Path structure: logs/benchmarks/{rule_name}/{sample}.txt
    parts = Path(filepath).parts
    if 'benchmarks' in parts:
        idx = parts.index('benchmarks')
        if idx + 1 < len(parts):
            return parts[idx + 1]
    return "unknown"


def aggregate_benchmarks(benchmark_dir: str) -> Tuple[Dict, List[Dict]]:
    """
    Aggregate all benchmark files in the directory.
    
    Returns:
        Tuple of (aggregated_stats, all_records)
    """
    all_records = []
    by_rule = defaultdict(list)
    
    # Find all benchmark files
    benchmark_files = glob.glob(os.path.join(benchmark_dir, '**/*.txt'), recursive=True)
    
    if not benchmark_files:
        print(f"Warning: No benchmark files found in {benchmark_dir}", file=sys.stderr)
        return {}, []
    
    print(f"Found {len(benchmark_files)} benchmark files", file=sys.stderr)
    
    # Parse all files
    for filepath in benchmark_files:
        data = parse_benchmark_file(filepath)
        if data:
            rule_name = get_rule_name_from_path(filepath)
            data['_rule'] = rule_name
            all_records.append(data)
            by_rule[rule_name].append(data)
    
    if not all_records:
        print("Warning: No valid benchmark data found", file=sys.stderr)
        return {}, []
    
    # Aggregate statistics
    stats = {
        'total_files': len(all_records),
        'total_rules': len(by_rule),
        'by_rule': {},
        'overall': {}
    }
    
    # Calculate overall statistics
    all_runtime = [r['s'] for r in all_records if 's' in r]
    all_max_rss = [r.get('max_rss', 0) for r in all_records]
    all_max_vms = [r.get('max_vms', 0) for r in all_records]
    all_io_in = [r.get('io_in', 0) for r in all_records]
    all_io_out = [r.get('io_out', 0) for r in all_records]
    all_cpu_time = [r.get('cpu_time', 0) for r in all_records]
    
    stats['overall'] = {
        'total_runtime_seconds': sum(all_runtime),
        'total_runtime_hours': sum(all_runtime) / 3600,
        'total_cpu_time_seconds': sum(all_cpu_time),
        'total_cpu_time_hours': sum(all_cpu_time) / 3600,
        'max_memory_mb': max(all_max_rss) if all_max_rss else 0,
        'total_memory_mb': sum(all_max_rss),
        'avg_memory_mb': statistics.mean(all_max_rss) if all_max_rss else 0,
        'total_io_in_mb': sum(all_io_in),
        'total_io_out_mb': sum(all_io_out),
        'sample_count': len(all_records)
    }
    
    # Calculate per-rule statistics
    for rule_name, records in by_rule.items():
        runtime = [r['s'] for r in records if 's' in r]
        max_rss = [r.get('max_rss', 0) for r in records]
        cpu_time = [r.get('cpu_time', 0) for r in records]
        
        stats['by_rule'][rule_name] = {
            'count': len(records),
            'total_runtime_seconds': sum(runtime),
            'total_runtime_hours': sum(runtime) / 3600,
            'avg_runtime_seconds': statistics.mean(runtime) if runtime else 0,
            'min_runtime_seconds': min(runtime) if runtime else 0,
            'max_runtime_seconds': max(runtime) if runtime else 0,
            'total_cpu_time_seconds': sum(cpu_time),
            'avg_memory_mb': statistics.mean(max_rss) if max_rss else 0,
            'max_memory_mb': max(max_rss) if max_rss else 0,
            'total_io_in_mb': sum(r.get('io_in', 0) for r in records),
            'total_io_out_mb': sum(r.get('io_out', 0) for r in records),
        }
    
    return stats, all_records


def format_time(seconds: float) -> str:
    """Format seconds to human-readable time."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    
    if hours > 0:
        return f"{hours}h {minutes}m {secs}s"
    elif minutes > 0:
        return f"{minutes}m {secs}s"
    else:
        return f"{secs}s"


def format_memory(mb: float) -> str:
    """Format memory in MB to human-readable format."""
    if mb > 1024:
        return f"{mb/1024:.2f} GB"
    else:
        return f"{mb:.2f} MB"


def print_table_format(stats: Dict, all_records: List[Dict]):
    """Print statistics in table format."""
    print("\n" + "="*80)
    print("OVERALL BENCHMARK STATISTICS")
    print("="*80)
    
    overall = stats['overall']
    print(f"\nTotal Samples Analyzed: {overall['sample_count']}")
    print(f"Total Rules: {stats['total_rules']}")
    print(f"\nRuntime Summary:")
    print(f"  Total Runtime:     {format_time(overall['total_runtime_seconds'])} ({overall['total_runtime_hours']:.2f} hours)")
    print(f"  Total CPU Time:    {format_time(overall['total_cpu_time_seconds'])} ({overall['total_cpu_time_hours']:.2f} hours)")
    print(f"\nMemory Summary:")
    print(f"  Peak Memory:       {format_memory(overall['max_memory_mb'])}")
    print(f"  Average Memory:    {format_memory(overall['avg_memory_mb'])}")
    print(f"  Total Memory Used: {format_memory(overall['total_memory_mb'])}")
    print(f"\nI/O Summary:")
    print(f"  Total Input:       {format_memory(overall['total_io_in_mb'])}")
    print(f"  Total Output:      {format_memory(overall['total_io_out_mb'])}")
    
    # Sort rules by total runtime
    sorted_rules = sorted(stats['by_rule'].items(), 
                         key=lambda x: x[1]['total_runtime_seconds'], 
                         reverse=True)
    
    print("\n" + "="*80)
    print("PER-RULE STATISTICS (sorted by total runtime)")
    print("="*80)
    print(f"\n{'Rule Name':<30} {'Count':>8} {'Total Time':>15} {'Avg Time':>12} {'Avg Memory':>12}")
    print("-"*80)
    
    for rule_name, rule_stats in sorted_rules:
        print(f"{rule_name:<30} {rule_stats['count']:>8} "
              f"{format_time(rule_stats['total_runtime_seconds']):>15} "
              f"{format_time(rule_stats['avg_runtime_seconds']):>12} "
              f"{format_memory(rule_stats['avg_memory_mb']):>12}")
    
    # Top 10 longest running individual jobs
    print("\n" + "="*80)
    print("TOP 10 LONGEST RUNNING JOBS")
    print("="*80)
    
    sorted_records = sorted(all_records, key=lambda x: x.get('s', 0), reverse=True)[:10]
    print(f"\n{'Rule':<25} {'Sample':<25} {'Runtime':>12} {'Memory':>12}")
    print("-"*80)
    
    for record in sorted_records:
        sample_name = os.path.splitext(record['_filename'])[0]
        print(f"{record['_rule']:<25} {sample_name:<25} "
              f"{format_time(record.get('s', 0)):>12} "
              f"{format_memory(record.get('max_rss', 0)):>12}")


def print_json_format(stats: Dict, all_records: List[Dict]):
    """Print statistics in JSON format."""
    output = {
        'summary': stats,
        'details': all_records
    }
    print(json.dumps(output, indent=2))


def print_csv_format(stats: Dict, all_records: List[Dict]):
    """Print detailed records in CSV format."""
    if not all_records:
        return
    
    # Print header
    keys = ['_rule', '_filename', 's', 'max_rss', 'max_vms', 'io_in', 'io_out', 'cpu_time', 'mean_load']
    print(','.join(keys))
    
    # Print records
    for record in all_records:
        values = [str(record.get(k, '')) for k in keys]
        print(','.join(values))


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate Snakemake benchmark statistics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        'benchmark_dir',
        nargs='?',
        default='logs/benchmarks',
        help='Directory containing benchmark files (default: logs/benchmarks)'
    )
    parser.add_argument(
        '--format',
        choices=['table', 'json', 'csv'],
        default='table',
        help='Output format (default: table)'
    )
    parser.add_argument(
        '--top',
        type=int,
        default=10,
        help='Number of top jobs to show (default: 10)'
    )
    
    args = parser.parse_args()
    
    # Check if benchmark directory exists
    if not os.path.exists(args.benchmark_dir):
        print(f"Error: Benchmark directory '{args.benchmark_dir}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Aggregate benchmarks
    stats, all_records = aggregate_benchmarks(args.benchmark_dir)
    
    if not stats:
        print("No benchmark data to display", file=sys.stderr)
        sys.exit(1)
    
    # Output in requested format
    if args.format == 'json':
        print_json_format(stats, all_records)
    elif args.format == 'csv':
        print_csv_format(stats, all_records)
    else:
        print_table_format(stats, all_records)


if __name__ == '__main__':
    main()
