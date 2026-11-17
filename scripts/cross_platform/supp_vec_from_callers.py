#!/usr/bin/env python3
"""Derive INFO/SUPP_VEC from INFO/CALLERS in a cross-platform truvari merge.

The truvari branch of rules/cross_platform_sv.smk tags every input record with
INFO/CALLER=ngs|tgs before collapsing, and scripts/sv_pipeline_scripts/annotate_caller_support.py
then writes INFO/CALLERS (the distinct tags folded into each representative). This script turns
that into the same INFO/SUPP_VEC string the SURVIVOR branch produces, so the downstream
concordance summary and documented output format are identical for both methods:

    "11" -> both platforms, "10" -> NGS only, "01" -> TGS only

Records are written sorted in the header's contig order (truvari collapse output is unsorted
and tabix needs sorted input). Standard library only, so it runs in the truvari container
(python3, no bcftools); one sample's call set fits in memory.
"""
import argparse
import gzip
import sys

SUPP_VEC_HEADER = (
    '##INFO=<ID=SUPP_VEC,Number=1,Type=String,'
    'Description="Support vector over the input call sets in NGS,TGS order (1 = detected)">'
)
SUPP_HEADER = (
    '##INFO=<ID=SUPP,Number=1,Type=Integer,'
    'Description="Number of input call sets supporting this variant">'
)
LABELS = ("ngs", "tgs")


def open_text(path, mode="r"):
    if path == "-":
        return sys.stdin if mode == "r" else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", errors="surrogateescape")
    return open(path, mode, errors="surrogateescape")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", default="-")
    args = ap.parse_args()

    counts = {"11": 0, "10": 0, "01": 0, "00": 0}
    headers, records, contigs = [], [], {}
    with open_text(args.input) as fin:
        for line in fin:
            if line.startswith("##"):
                if "ID=SUPP_VEC," in line or "ID=SUPP," in line:
                    continue
                if line.startswith("##contig=<ID="):
                    name = line[len("##contig=<ID="):].split(",", 1)[0].rstrip(">\n")
                    contigs.setdefault(name, len(contigs))
                headers.append(line)
                continue
            if line.startswith("#CHROM"):
                headers.append(SUPP_VEC_HEADER + "\n" + SUPP_HEADER + "\n" + line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            callers = set()
            kept = []
            for token in fields[7].split(";"):
                if not token:
                    continue
                key, _, value = token.partition("=")
                if key == "CALLERS":
                    callers = {c.strip() for c in value.split(",") if c.strip()}
                if key in ("SUPP_VEC", "SUPP"):
                    continue
                kept.append(token)
            supp_vec = "".join("1" if lab in callers else "0" for lab in LABELS)
            counts[supp_vec] = counts.get(supp_vec, 0) + 1
            kept.append("SUPP=%d" % supp_vec.count("1"))
            kept.append("SUPP_VEC=%s" % supp_vec)
            fields[7] = ";".join(kept) if kept else "."
            records.append(fields)

    records.sort(key=lambda f: (contigs.get(f[0], len(contigs)), f[0], int(f[1]) if f[1].isdigit() else 0))
    with open_text(args.output, "w") as fout:
        fout.writelines(headers)
        for fields in records:
            fout.write("\t".join(fields) + "\n")

    sys.stderr.write("SUPP_VEC counts: both %d, NGS-only %d, TGS-only %d, unlabelled %d\n" % (
        counts["11"], counts["10"], counts["01"], counts["00"]))
    if counts["00"]:
        sys.stderr.write(
            "Warning: %d records carry no ngs/tgs tag; INFO/CALLER was expected on every input record\n"
            % counts["00"])


if __name__ == "__main__":
    main()
