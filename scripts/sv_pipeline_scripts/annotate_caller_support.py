#!/usr/bin/env python3
"""Annotate a truvari-collapsed VCF with true distinct-caller support.

Background
----------
`truvari collapse` writes INFO/NumCollapsed = the number of *records* folded into a
representative. That is NOT the number of distinct callers: one caller can emit
several near-identical records for the same event, so NumCollapsed over-counts
support (observed max 22 with only 4 callers on HG002 ONT).

truvari does, however, preserve a join key:
  * each representative (in the collapsed VCF) carries INFO/CollapseId
  * each folded/removed record (in the *.collapsed.removed.vcf) carries INFO/MatchId,
    whose value(s) equal the CollapseId of the representative it was folded into.

If every input record was first tagged with INFO/CALLER (see tag_caller_source.py),
we can recover the real support: for each representative, the distinct callers are
its own CALLER plus the CALLER of every removed record whose MatchId == its CollapseId.

This script emits two new INFO fields on the representatives:
  * CALLERS     = comma-separated sorted distinct caller names
  * NUM_CALLERS = number of distinct callers

Standard library only. First pass indexes the (small) removed file by CollapseId;
second pass streams the representatives, so peak memory stays proportional to the
number of collapse groups, not the callset size.
"""
import argparse
import gzip
import sys

MATCH_KEYS = {
    "PctSeqSimilarity", "PctSizeSimilarity", "PctRecOverlap", "SizeDiff",
    "StartDistance", "EndDistance", "TruScore", "MatchId",
}
GROUP_KEYS = ("CollapseId", "NumCollapsed", "NumConsolidated")

CALLERS_HEADER = (
    '##INFO=<ID=CALLERS,Number=.,Type=String,'
    'Description="Distinct SV callers supporting this variant, recovered via '
    'truvari CollapseId/MatchId over per-record INFO/CALLER tags">'
)
NUM_CALLERS_HEADER = (
    '##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,'
    'Description="Number of distinct SV callers supporting this variant">'
)


def open_text(path, mode):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", errors="surrogateescape")
    return open(path, mode, errors="surrogateescape")


def parse_info(info):
    """Return a dict of INFO key->value (flags map to True)."""
    d = {}
    if info in (".", ""):
        return d
    for token in info.split(";"):
        if not token:
            continue
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
        else:
            d[token] = True
    return d


def caller_of(info_dict):
    val = info_dict.get("CALLER")
    if isinstance(val, str) and val:
        return val
    return None


def rank(fields, info, priority):
    caller = caller_of(info)
    order = priority.index(caller) if caller in priority else len(priority)
    try:
        qual = float(fields[5])
    except ValueError:
        qual = 0.0
    return (order, -qual)


def index_removed(removed_path, priority=None):
    """Map each CollapseId -> set of callers of records folded under it."""
    groups, best = {}, {}
    with open_text(removed_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            caller = caller_of(info)
            match_id = info.get("MatchId")
            if not match_id or caller is None:
                continue
            # MatchId is Number=. ; a record may reference several groups.
            for cid in str(match_id).split(","):
                cid = cid.strip()
                if not cid or cid == ".":
                    continue
                groups.setdefault(cid, set()).add(caller)
                if priority:
                    key = rank(fields, info, priority)
                    if cid not in best or key < best[cid][0]:
                        best[cid] = (key, fields, info)
    return groups, best


def as_representative(fields, rep_info):
    kept = [t for t in fields[7].split(";") if t and t.split("=", 1)[0] not in MATCH_KEYS]
    kept += ["%s=%s" % (k, rep_info[k]) for k in GROUP_KEYS if k in rep_info]
    fields = list(fields)
    fields[7] = ";".join(kept) if kept else "."
    return fields


def annotate(collapsed_path, groups, out_fh, best=None, priority=None):
    best = best or {}
    header_emitted = False
    with open_text(collapsed_path, "r") as fh:
        for line in fh:
            if line.startswith("##"):
                out_fh.write(line)
                continue
            if line.startswith("#CHROM"):
                out_fh.write(CALLERS_HEADER + "\n")
                out_fh.write(NUM_CALLERS_HEADER + "\n")
                out_fh.write(line)
                header_emitted = True
                continue
            if not header_emitted:
                out_fh.write(CALLERS_HEADER + "\n")
                out_fh.write(NUM_CALLERS_HEADER + "\n")
                out_fh.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                )
                header_emitted = True
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                out_fh.write(line)
                continue
            info = parse_info(fields[7])
            callers = set()
            own = caller_of(info)
            if own is not None:
                callers.add(own)
            cid = info.get("CollapseId")
            if cid and cid != ".":
                callers |= groups.get(cid, set())
                if cid in best and best[cid][0] < rank(fields, info, priority):
                    fields = as_representative(best[cid][1], info)
            # A representative with no CALLER tag and no folded records still counts
            # as one observed variant.
            num = len(callers) if callers else 1
            callers_str = ",".join(sorted(callers)) if callers else "unknown"
            extra = "CALLERS=%s;NUM_CALLERS=%d" % (callers_str, num)
            if fields[7] in (".", ""):
                fields[7] = extra
            else:
                fields[7] = fields[7] + ";" + extra
            out_fh.write("\t".join(fields) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--collapsed", required=True, help="Representative VCF (has CollapseId)")
    ap.add_argument("--removed", required=True, help="collapsed.removed VCF (has MatchId)")
    ap.add_argument("--output", default="-", help="Output VCF (default stdout)")
    ap.add_argument("--priority", default="",
                    help="Comma-separated caller order for choosing group representatives "
                         "(default: keep Truvari's representative)")
    args = ap.parse_args()

    priority = [c for c in args.priority.split(",") if c]
    groups, best = index_removed(args.removed, priority)
    out_fh = open_text(args.output, "w")
    try:
        annotate(args.collapsed, groups, out_fh, best, priority)
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()


if __name__ == "__main__":
    main()
