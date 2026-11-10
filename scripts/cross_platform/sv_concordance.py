#!/usr/bin/env python3
"""Summarize NGS/TGS concordance from a SURVIVOR-merged SV VCF.

SURVIVOR annotates every merged record with an INFO/SUPP_VEC field: a string of
'0'/'1' flags, one per input VCF, in the exact order the inputs were listed in the
merge file. For a cross-platform merge whose merge list is [NGS, TGS], SUPP_VEC is:
    "11" -> supported by both platforms (concordant)
    "10" -> NGS-only
    "01" -> TGS-only

This script parses SUPP_VEC (and SVTYPE) and reports counts / percentages broken
down by SV type. Uses only the Python standard library so it can run in any
container that provides python3.
"""

import argparse
import gzip
import sys
from collections import defaultdict


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _parse_info(info_str):
    info = {}
    for field in info_str.split(";"):
        if not field:
            continue
        if "=" in field:
            key, _, value = field.partition("=")
            info[key] = value
        else:
            info[field] = True
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="SURVIVOR-merged VCF (.vcf or .vcf.gz)")
    ap.add_argument("--output", required=True, help="Output concordance TSV")
    ap.add_argument("--pair-id", default="", help="Cross-platform pair identifier")
    ap.add_argument("--ngs-index", type=int, default=0,
                    help="0-based position of the NGS input in SUPP_VEC (default 0)")
    ap.add_argument("--tgs-index", type=int, default=1,
                    help="0-based position of the TGS input in SUPP_VEC (default 1)")
    args = ap.parse_args()

    # counts[svtype] = [both, ngs_only, tgs_only]
    counts = defaultdict(lambda: [0, 0, 0])
    missing_supp_vec = 0

    with _open(args.input) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            info = _parse_info(cols[7])
            supp_vec = info.get("SUPP_VEC")
            svtype = info.get("SVTYPE", "UNKNOWN")
            if not isinstance(supp_vec, str):
                missing_supp_vec += 1
                continue
            # Guard against SUPP_VEC shorter than expected
            ngs = args.ngs_index < len(supp_vec) and supp_vec[args.ngs_index] == "1"
            tgs = args.tgs_index < len(supp_vec) and supp_vec[args.tgs_index] == "1"
            if ngs and tgs:
                counts[svtype][0] += 1
            elif ngs:
                counts[svtype][1] += 1
            elif tgs:
                counts[svtype][2] += 1

    def _pct(part, whole):
        return f"{(100.0 * part / whole):.2f}" if whole else "0.00"

    # Totals across all SV types
    tot_both = sum(v[0] for v in counts.values())
    tot_ngs = sum(v[1] for v in counts.values())
    tot_tgs = sum(v[2] for v in counts.values())
    grand_total = tot_both + tot_ngs + tot_tgs

    with open(args.output, "w") as out:
        out.write("pair_id\tsvtype\tboth\tngs_only\ttgs_only\ttotal\t"
                  "concordance_pct\tngs_only_pct\ttgs_only_pct\n")

        def _row(svtype, both, ngs_only, tgs_only):
            total = both + ngs_only + tgs_only
            out.write(
                f"{args.pair_id}\t{svtype}\t{both}\t{ngs_only}\t{tgs_only}\t{total}\t"
                f"{_pct(both, total)}\t{_pct(ngs_only, total)}\t{_pct(tgs_only, total)}\n"
            )

        for svtype in sorted(counts):
            both, ngs_only, tgs_only = counts[svtype]
            _row(svtype, both, ngs_only, tgs_only)
        _row("ALL", tot_both, tot_ngs, tot_tgs)

    # Human-readable summary to stderr (captured in the rule log)
    sys.stderr.write(
        f"[concordance] pair={args.pair_id or 'NA'} total={grand_total} "
        f"both={tot_both} ngs_only={tot_ngs} tgs_only={tot_tgs} "
        f"concordance={_pct(tot_both, grand_total)}%\n"
    )
    if missing_supp_vec:
        sys.stderr.write(
            f"[concordance] warning: {missing_supp_vec} records had no SUPP_VEC and were skipped\n"
        )


if __name__ == "__main__":
    main()
