#!/usr/bin/env python3
"""
Generate sample map for GATK GenomicsDBImport from TGS samples.
"""

import argparse
import ast
import sys


def parse_samples(samples_str):
    """Parse samples string into a list."""
    samples = []
    if samples_str.strip():
        try:
            parsed = ast.literal_eval(samples_str)
        except (ValueError, SyntaxError):
            parsed = samples_str

        if isinstance(parsed, (list, tuple, set)):
            samples = [str(s).strip() for s in parsed if str(s).strip()]
        elif isinstance(parsed, str):
            samples = [tok for tok in parsed.split() if tok]
        elif parsed is not None:
            parsed_str = str(parsed).strip()
            if parsed_str:
                samples = [parsed_str]
    return samples


def main():
    parser = argparse.ArgumentParser(description="Generate sample map for GenomicsDBImport")
    parser.add_argument("--chromosome", required=True, help="Chromosome name")
    parser.add_argument("--output", required=True, help="Output sample map file")
    parser.add_argument("--final-dir", required=True, help="Final directory path")
    parser.add_argument("--tool", required=True, help="Variant calling tool (clair3 or deepvariant)")
    parser.add_argument("--samples", required=True, help="Samples list")
    
    args = parser.parse_args()
    
    samples = parse_samples(args.samples)
    
    if not samples:
        print("Error: No samples provided", file=sys.stderr)
        sys.exit(1)
    
    with open(args.output, "w") as out_f:
        for s in samples:
            if args.tool == "clair3":
                path = f"{args.final_dir}/{s}/clair3/{s}.clair3.gvcf.gz"
            elif args.tool == "deepvariant":
                path = f"{args.final_dir}/{s}/deepvariant/{s}.deepvariant.gatk.gvcf.gz"
            else:
                print(f"Error: Unsupported tool: {args.tool}", file=sys.stderr)
                sys.exit(1)
            out_f.write(f"{s}\t{path}\n")
    
    print(f"Generated sample map for {len(samples)} samples: {args.output}")


if __name__ == "__main__":
    main()
