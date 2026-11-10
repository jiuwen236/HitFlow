#!/usr/bin/env python3
"""
Fix DeepVariant GVCF for GATK compatibility by replacing <*> with <NON_REF>
and ensuring the NON_REF header is present.
"""

import gzip
import sys
import os
from datetime import datetime


def main():
    if len(sys.argv) != 4:
        print("Usage: fix_deepvariant_gvcf.py <input.gvcf.gz> <output.txt> <log.txt>")
        sys.exit(1)
    
    src = sys.argv[1]
    dst = sys.argv[2]
    log_path = sys.argv[3]
    
    non_ref_header = '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">'
    has_non_ref_header = False
    added_header = False
    modified_records = 0
    
    with gzip.open(src, "rt") as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("##"):
                if line.startswith("##ALT=<ID=NON_REF"):
                    has_non_ref_header = True
                fout.write(line)
                continue
            
            if line.startswith("#"):
                if not has_non_ref_header:
                    fout.write(non_ref_header + '\n')
                    has_non_ref_header = True
                    added_header = True
                fout.write(line)
                continue
            
            if "<*>" in line:
                line = line.replace("<*>", "<NON_REF>")
                modified_records += 1
            
            fout.write(line)
    
    with open(log_path, "w") as log_f:
        timestamp = datetime.now().isoformat()
        log_f.write(f"[{timestamp}] Processed {src}\n")
        log_f.write(f"Added NON_REF header: {added_header}\n")
        log_f.write(f"Records modified: {modified_records}\n")
    
    print(f"Successfully processed {src}")
    print(f"Added NON_REF header: {added_header}")
    print(f"Records modified: {modified_records}")


if __name__ == "__main__":
    main()
