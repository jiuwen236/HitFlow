#!/usr/bin/env python3
"""Annotate a cross-platform cohort SV merge and summarize platform support.

Reads the merged output of SURVIVOR merge or truvari collapse over two site-only cohort
VCFs (short-read/NGS and long-read/TGS), and writes to stdout a sorted, site-only VCF in
which every record carries:

  PLATFORM   both | short_only | long_only
  SUPP_VEC   two flags in NGS,TGS order ("11", "10", "01"), as in the per-sample module
  NGS_AF / NGS_AC / NGS_AN, TGS_AF / TGS_AC / TGS_AN
             allele frequency and counts in each source cohort (missing when that platform
             did not detect the site)
  NGS_ID / TGS_ID
             the matched site IDs in the source cohort VCFs

The two cohorts need not share individuals, so there is no joint genotype matrix: each
site keeps the allele frequency of the cohort(s) in which it was observed.

How support is recovered:
  * SURVIVOR writes INFO/SUPP_VEC itself, and per-sample FORMAT/ID gives the matched input
    IDs; both are re-derived here from the ID columns so the two methods agree.
  * truvari collapse keeps one representative per cluster (INFO/CollapseId) and writes the
    folded records to a separate file (INFO/MatchId); joining them over the site IDs yields
    the same support vector.

Standard library only, so it runs in any container that provides python3. Memory stays
proportional to the number of sites, not to cohort size, because genotypes were already
dropped upstream.
"""
import argparse
import gzip
import sys
from collections import Counter, defaultdict

NGS, TGS = "NGS", "TGS"
NEW_HEADERS = [
    '##INFO=<ID=PLATFORM,Number=1,Type=String,Description="Platform support: both, short_only or long_only">',
    '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Support vector over the input cohorts in NGS,TGS order (1 = detected)">',
    '##INFO=<ID=NGS_ID,Number=.,Type=String,Description="Matched site IDs in the short-read cohort VCF">',
    '##INFO=<ID=TGS_ID,Number=.,Type=String,Description="Matched site IDs in the long-read cohort VCF">',
    '##INFO=<ID=NGS_AF,Number=1,Type=Float,Description="Allele frequency in the short-read cohort">',
    '##INFO=<ID=NGS_AC,Number=1,Type=Integer,Description="Allele count in the short-read cohort">',
    '##INFO=<ID=NGS_AN,Number=1,Type=Integer,Description="Allele number in the short-read cohort">',
    '##INFO=<ID=TGS_AF,Number=1,Type=Float,Description="Allele frequency in the long-read cohort">',
    '##INFO=<ID=TGS_AC,Number=1,Type=Integer,Description="Allele count in the long-read cohort">',
    '##INFO=<ID=TGS_AN,Number=1,Type=Integer,Description="Allele number in the long-read cohort">',
]
# INFO keys written by the merge tools that describe the merge, not the variant.
DROP_KEYS = {
    "SUPP", "SUPP_VEC", "SVMETHOD", "CollapseId", "NumCollapsed", "NumConsolidated",
    "MatchId", "PctSeqSimilarity", "PctSizeSimilarity", "PctRecOverlap", "SizeDiff",
    "StartDistance", "EndDistance", "TruScore", "AN", "AC", "AF", "SRC_ID",
}


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="surrogateescape")
    return open(path, "rt", errors="surrogateescape")


def parse_info(info):
    out = {}
    if info in (".", ""):
        return out
    for token in info.split(";"):
        if not token:
            continue
        key, _, value = token.partition("=")
        out[key] = value if _ else True
    return out


def platform_of(site_id):
    """Site IDs were stamped by prepare_sv_sites.awk as NGS_<n> / TGS_<n>."""
    if site_id.startswith(NGS + "_"):
        return NGS
    if site_id.startswith(TGS + "_"):
        return TGS
    return None


def load_sites(path, platform, stats):
    """Map site ID -> (AF, AC, AN, source ID) for one platform."""
    table = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            table[fields[2]] = (
                info.get("AF", "."), info.get("AC", "."), info.get("AN", "."),
                info.get("SRC_ID", fields[2]),
            )
    stats[platform] = len(table)
    return table


def collapse_groups(removed_path):
    """truvari: map a representative's CollapseId -> site IDs folded into it."""
    groups = defaultdict(set)
    if not removed_path:
        return groups
    with open_text(removed_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            match_id = info.get("MatchId")
            if not isinstance(match_id, str):
                continue
            for cid in match_id.split(","):
                cid = cid.strip()
                if cid and cid != ".":
                    groups[cid].add(fields[2])
    return groups


def member_ids(fields, info, method, groups):
    """Site IDs of every input record represented by this merged record."""
    ids = {fields[2]}
    if method == "survivor":
        # SURVIVOR carries the matched input IDs in the per-sample FORMAT/ID field.
        if len(fields) > 9 and "ID" in fields[8].split(":"):
            idx = fields[8].split(":").index("ID")
            for sample in fields[9:]:
                parts = sample.split(":")
                if len(parts) > idx:
                    value = parts[idx]
                    if value not in ("", ".", "NA", "NaN"):
                        ids.update(p for p in value.split(",") if platform_of(p))
    else:
        cid = info.get("CollapseId")
        if isinstance(cid, str):
            ids |= groups.get(cid, set())
    return {i for i in ids if platform_of(i)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged", required=True, help="SURVIVOR merge or truvari collapse output VCF")
    ap.add_argument("--removed", help="truvari collapse removed-records VCF (truvari only)")
    ap.add_argument("--ngs-sites", required=True)
    ap.add_argument("--tgs-sites", required=True)
    ap.add_argument("--method", choices=("survivor", "truvari"), default="survivor")
    ap.add_argument("--report-maf", type=float, default=0.01,
                    help="also count long-read-only sites above this allele frequency")
    ap.add_argument("--fai", help="reference .fai; output is sorted in reference contig order")
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    stats = {}
    sites = {
        NGS: load_sites(args.ngs_sites, NGS, stats),
        TGS: load_sites(args.tgs_sites, TGS, stats),
    }
    groups = collapse_groups(args.removed) if args.method == "truvari" else defaultdict(set)

    headers, records = [], []
    by_platform, by_type = Counter(), Counter()
    long_only_common = 0
    matched = {NGS: set(), TGS: set()}

    with open_text(args.merged) as fh:
        for line in fh:
            if line.startswith("##"):
                if not any(f"ID={k}," in line for k in ("PLATFORM", "SUPP_VEC", "NGS_AF", "TGS_AF")):
                    headers.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            ids = member_ids(fields, info, args.method, groups)
            per_platform = {p: sorted(i for i in ids if platform_of(i) == p) for p in (NGS, TGS)}
            if not per_platform[NGS] and not per_platform[TGS]:
                continue
            for p in (NGS, TGS):
                matched[p].update(per_platform[p])

            supp_vec = "".join("1" if per_platform[p] else "0" for p in (NGS, TGS))
            platform = {"11": "both", "10": "short_only", "01": "long_only"}[supp_vec]
            kept = [f"{k}={v}" if v is not True else k
                    for k, v in info.items() if k not in DROP_KEYS]
            kept.append(f"PLATFORM={platform}")
            kept.append(f"SUPP_VEC={supp_vec}")
            for p in (NGS, TGS):
                if not per_platform[p]:
                    continue
                # A merged site can fold several records of one cohort; report the most
                # frequent allele and list every matched source ID.
                best = max(per_platform[p],
                           key=lambda i: (float(sites[p][i][0]) if _is_float(sites[p][i][0]) else -1.0))
                af, ac, an, _src = sites[p][best]
                kept.append("%s_ID=%s" % (p, ",".join(sites[p][i][3] for i in per_platform[p])))
                for key, value in (("AF", af), ("AC", ac), ("AN", an)):
                    if value not in (".", ""):
                        kept.append(f"{p}_{key}={value}")
                if p == TGS and platform == "long_only" and _is_float(af) and float(af) > args.report_maf:
                    long_only_common += 1

            by_platform[platform] += 1
            sv_type = info.get("SVTYPE", "NA")
            by_type[(sv_type if isinstance(sv_type, str) else "NA", platform)] += 1
            fields[7] = ";".join(kept) if kept else "."
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            records.append((fields[0], pos, fields[:8]))

    out = sys.stdout
    for line in headers:
        out.write(line + "\n")
    for line in NEW_HEADERS:
        out.write(line + "\n")
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    contig_order = {}
    if args.fai:
        with open(args.fai) as fh:
            for line in fh:
                contig_order.setdefault(line.split("\t", 1)[0], len(contig_order))
    for line in headers:
        if line.startswith("##contig=<ID="):
            name = line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">")
            contig_order.setdefault(name, len(contig_order))
    records.sort(key=lambda r: (contig_order.get(r[0], len(contig_order)), r[0], r[1]))
    for _chrom, _pos, fields in records:
        out.write("\t".join(fields) + "\n")

    total = sum(by_platform.values())
    with open(args.summary, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"merge_method\t{args.method}\n")
        fh.write(f"short_read_cohort_sites\t{stats[NGS]}\n")
        fh.write(f"long_read_cohort_sites\t{stats[TGS]}\n")
        fh.write(f"integrated_sites\t{total}\n")
        for key in ("both", "short_only", "long_only"):
            pct = 100.0 * by_platform[key] / total if total else 0.0
            fh.write(f"{key}\t{by_platform[key]}\t{pct:.2f}%\n")
        fh.write(f"long_only_af_above_{args.report_maf}\t{long_only_common}\n")
        fh.write(f"short_read_sites_matched\t{len(matched[NGS])}\n")
        fh.write(f"long_read_sites_matched\t{len(matched[TGS])}\n")
        fh.write("\nsvtype\tboth\tshort_only\tlong_only\n")
        for sv_type in sorted({t for t, _ in by_type}):
            fh.write("%s\t%d\t%d\t%d\n" % (
                sv_type, by_type[(sv_type, "both")],
                by_type[(sv_type, "short_only")], by_type[(sv_type, "long_only")]))

    sys.stderr.write(
        "Integrated %d sites: both %d, short-only %d, long-only %d (long-only AF>%.3g: %d)\n" % (
            total, by_platform["both"], by_platform["short_only"],
            by_platform["long_only"], args.report_maf, long_only_common))
    # Every input site must end up in exactly one integrated record.
    for p, path in ((NGS, args.ngs_sites), (TGS, args.tgs_sites)):
        if len(matched[p]) != stats[p]:
            sys.stderr.write(
                "Note: %d of %d %s sites are listed in the output (%s); the rest were filtered by "
                "the merge tool (size/type) or, with SURVIVOR, folded into a nearby site of the "
                "same cohort without keeping their ID\n" % (len(matched[p]), stats[p], p, path))


def _is_float(value):
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    main()
