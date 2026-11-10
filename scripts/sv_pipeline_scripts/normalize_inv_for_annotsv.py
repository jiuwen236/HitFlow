#!/usr/bin/env python3
"""Normalize sequence-resolved inversions so AnnotSV can size them.

Some long-read callers (SVIM in particular) emit inversions as an explicit
sequence swap: REF = the forward strand, ALT = its reverse complement, with equal
lengths and NO INFO/END or INFO/SVLEN. AnnotSV computes variant length from
END-POS (or SVLEN); with neither present it derives length 0, decides
`variantLength (0) < SVminSize (30)`, and aborts the whole annotation run.

This script rewrites ONLY those records into symbolic-inversion form that AnnotSV
understands, leaving everything else byte-for-byte unchanged:
    ALT   -> <INV>
    END   =  POS + len(REF) - 1
    SVLEN =  len(REF)
It is intentionally conservative: a record is rewritten only when ALL hold:
    * INFO/SVTYPE == INV
    * INFO/END and INFO/SVLEN are both absent
    * len(REF) == len(ALT)
    * ALT == reverse_complement(REF)
This guarantees we never mislabel a non-inversion or clobber records that already
carry proper END/SVLEN. Standard library only; streams line by line.
"""
import argparse
import gzip
import sys

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

END_HEADER = (
    '##INFO=<ID=END,Number=1,Type=Integer,'
    'Description="End position of the variant described in this record">'
)
SVLEN_HEADER = (
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,'
    'Description="Difference in length between REF and ALT alleles">'
)


def open_text(path, mode):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def parse_info(info):
    d = {}
    order = []
    if info in (".", ""):
        return d, order
    for token in info.split(";"):
        if not token:
            continue
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
        else:
            k = token
            d[k] = None  # flag
        order.append(k)
    return d, order


def build_info(info_dict, order):
    parts = []
    for k in order:
        v = info_dict[k]
        parts.append(k if v is None else "%s=%s" % (k, v))
    return ";".join(parts) if parts else "."


def needs_normalization(ref, alt, info_dict):
    if info_dict.get("SVTYPE") != "INV":
        return False
    if "END" in info_dict or "SVLEN" in info_dict:
        return False
    # Symbolic ALTs (already <INV>) have '<' and won't equal a revcomp.
    if not ref or not alt or "<" in alt or "[" in alt or "]" in alt:
        return False
    if len(ref) != len(alt):
        return False
    return alt.upper() == reverse_complement(ref).upper()


def normalize(in_fh, out_fh):
    header_emitted = False
    rewritten = 0
    for line in in_fh:
        if line.startswith("##"):
            out_fh.write(line)
            continue
        if line.startswith("#CHROM"):
            # Make sure END/SVLEN are declared so AnnotSV parses the rewritten records.
            out_fh.write(END_HEADER + "\n")
            out_fh.write(SVLEN_HEADER + "\n")
            out_fh.write(line)
            header_emitted = True
            continue
        if not header_emitted:
            out_fh.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            out_fh.write(line)
            continue
        ref, alt = fields[3], fields[4]
        info_dict, order = parse_info(fields[7])
        if needs_normalization(ref, alt, info_dict):
            pos = int(fields[1])
            end = pos + len(ref) - 1
            info_dict["END"] = str(end)
            info_dict["SVLEN"] = str(len(ref))
            if "END" not in order:
                order.append("END")
            if "SVLEN" not in order:
                order.append("SVLEN")
            fields[4] = "<INV>"
            fields[7] = build_info(info_dict, order)
            rewritten += 1
        out_fh.write("\t".join(fields) + "\n")
    return rewritten


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="-", help="Input VCF (default stdin)")
    ap.add_argument("--output", default="-", help="Output VCF (default stdout)")
    ap.add_argument("--report", default=None, help="Optional path to write a rewrite count")
    args = ap.parse_args()

    in_fh = open_text(args.input, "r")
    out_fh = open_text(args.output, "w")
    try:
        rewritten = normalize(in_fh, out_fh)
    finally:
        if in_fh is not sys.stdin:
            in_fh.close()
        if out_fh is not sys.stdout:
            out_fh.close()
    if args.report:
        with open(args.report, "w") as rep:
            rep.write("normalized_inv_records=%d\n" % rewritten)


if __name__ == "__main__":
    main()
