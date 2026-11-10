#!/usr/bin/env python3
"""Tag every VCF record with its source caller via INFO/CALLER.

truvari collapse annotates INFO/NumCollapsed (a folded-record count) but does not
preserve which caller produced each record. To recover true distinct-caller support
after collapsing, we stamp each per-caller VCF with INFO/CALLER=<name> BEFORE the
records are concatenated and collapsed. annotate_caller_support.py then reads this
tag back through the CollapseId/MatchId link to compute NUM_CALLERS.

Standard library only (runs in the bcftools/cyvcf2 biocontainers). Streams line by
line so it scales to whole-genome, many-sample cohorts. Only one INFO field is added;
all other columns are passed through unchanged.
"""
import argparse
import gzip
import sys

CALLER_INFO_HEADER = (
    '##INFO=<ID=CALLER,Number=1,Type=String,'
    'Description="Source SV caller for this record (added before multi-caller collapse)">'
)


def open_text(path, mode):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", errors="surrogateescape")
    return open(path, mode, errors="surrogateescape")


def tag_records(in_fh, out_fh, caller):
    """Copy VCF from in_fh to out_fh, adding INFO/CALLER=<caller> to each record."""
    header_done = False
    caller_field = "CALLER=" + caller
    for line in in_fh:
        if line.startswith("##"):
            out_fh.write(line)
            continue
        if line.startswith("#CHROM"):
            # Emit the CALLER INFO definition once, right before the column header.
            out_fh.write(CALLER_INFO_HEADER + "\n")
            out_fh.write(line)
            header_done = True
            continue
        if not header_done:
            # Malformed VCF (record before #CHROM); still ensure the header exists.
            out_fh.write(CALLER_INFO_HEADER + "\n")
            out_fh.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            )
            header_done = True
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            # Not a well-formed data line; pass through untouched.
            out_fh.write(line)
            continue
        info = fields[7]
        if info in (".", ""):
            fields[7] = caller_field
        else:
            fields[7] = info + ";" + caller_field
        out_fh.write("\t".join(fields) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--caller", required=True, help="Caller name to stamp (e.g. sniffles)")
    ap.add_argument("--input", default="-", help="Input VCF (default stdin)")
    ap.add_argument("--output", default="-", help="Output VCF (default stdout)")
    args = ap.parse_args()

    in_fh = open_text(args.input, "r")
    out_fh = open_text(args.output, "w")
    try:
        tag_records(in_fh, out_fh, args.caller)
    finally:
        if in_fh is not sys.stdin:
            in_fh.close()
        if out_fh is not sys.stdout:
            out_fh.close()


if __name__ == "__main__":
    main()
