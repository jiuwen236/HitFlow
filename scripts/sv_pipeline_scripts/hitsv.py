#!/usr/bin/env python3
"""Validate HitSV inputs, dispatch one chromosome, and audit its VCF output."""
import argparse
import collections
import json
import os
import re
import subprocess
from pathlib import Path


def preflight(args):
    import pysam
    fai = {}
    with open(args.reference + ".fai") as handle:
        for line in handle:
            fields = line.split()
            fai[fields[0]] = int(fields[1])
    chromosomes = args.chromosomes.split(",")
    if len(set(chromosomes)) != len(chromosomes) or not chromosomes:
        raise ValueError("Chromosomes must be a nonempty unique list")
    dictionaries = []
    inputs = []
    for bam, index in zip(args.bam, args.index):
        with pysam.AlignmentFile(bam, "rb", index_filename=index) as handle:
            handle.check_index()
            header = handle.header.to_dict()
            if header.get("HD", {}).get("SO") != "coordinate":
                raise ValueError("BAM must be coordinate sorted: " + bam)
            dictionary = list(zip(handle.references, handle.lengths))
            for name, length in dictionary:
                if fai.get(name) != length:
                    raise ValueError("BAM/FAI dictionary mismatch: " + name)
            ids = {name: i for i, (name, _) in enumerate(dictionary)}
            if any(name not in ids for name in chromosomes):
                raise ValueError("Requested chromosome missing in BAM/FAI")
            samples = sorted({rg["SM"] for rg in header.get("RG", []) if rg.get("SM")})
            if not samples:
                raise ValueError("BAM requires an RG/SM sample identity: " + bam)
            inputs.append({"bam": os.path.realpath(bam), "index": os.path.realpath(index),
                           "samples": samples, "size": os.stat(bam).st_size,
                           "mtime_ns": os.stat(bam).st_mtime_ns})
            dictionaries.append(dictionary)
    if len(inputs) == 2:
        if dictionaries[0] != dictionaries[1]:
            raise ValueError("Hybrid BAM dictionaries must have identical order and lengths")
        # A declared biological identity is mandatory when BAM SM tags differ.
        aliases = json.loads(args.identities)
        identities = [{aliases.get(sm, sm) for sm in item["samples"]} for item in inputs]
        if len(identities[0]) != 1 or identities[0] != identities[1]:
            raise ValueError("Hybrid BAMs are not the same individual; verify SM tags/identity aliases")
    ids = {name: i for i, (name, _) in enumerate(dictionaries[0])}
    payload = {"reference": os.path.realpath(args.reference), "inputs": inputs,
               "chromosomes": {name: ids[name] for name in chromosomes}}
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n")


def call(args):
    manifest = json.loads(Path(args.manifest).read_text())
    chromosome_id = manifest["chromosomes"][args.chromosome]
    partial = args.output + ".partial"
    Path(partial).unlink(missing_ok=True)
    command = [args.binary, "call", "-S", str(chromosome_id), "-E", str(chromosome_id),
               "-r", args.reference, "-o", partial, "-m", str(args.min_size),
               "-v", "0", "-f", "0", "-p", args.preset, "-P", args.alignment_preset]
    if args.srs:
        for stat_path in (args.fa_stat, args.bam_stat):
            if Path(stat_path).stat().st_size == 0:
                raise ValueError("Empty SRS statistics file: " + stat_path)
        command += ["-n", args.srs, "-L", args.tl, "-I", args.fa_stat, "-T", args.bam_stat]
    if args.lrs:
        command += ["-l", args.lrs]
    if args.error_profile:
        # Required by HitSV for SRS polishing of LRS contigs in Hybrid mode.
        command += ["-k", args.error_profile]
    print(json.dumps(command), flush=True)
    subprocess.run(command, check=True)
    # A successful exit without a complete header is a failed task, not an empty callset.
    with open(partial) as handle:
        if not any(line.startswith("#CHROM\t") for line in handle):
            raise ValueError("HitSV did not emit a VCF header")
    os.replace(partial, args.output)


class AlleleValidationError(ValueError):
    """An allele is internally inconsistent or cannot be verified against FASTA."""


def _padding_evidence(fields, info):
    """Validate the local indel against the retained contig alignment."""
    cigar = info.get("CIGAR", "")
    contig = info.get("CONTIG", "").upper()
    location = info.get("CONTIG_POSITION", "")
    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not operations or "".join(n + op for n, op in operations) != cigar:
        return False
    if not location.startswith(fields[0] + ":"):
        return False
    try:
        refpos = int(location.rsplit(":", 1)[1])
        size = abs(int(info["SVLEN"]))
    except (ValueError, KeyError):
        return False
    querypos = 0
    candidates = []
    # The zero-based indel cursor equals the one-based preceding VCF anchor.
    for n, op in operations:
        length = int(n)
        if length <= 0:
            return False
        if refpos == int(fields[1]) and ((op == "I" and info["SVTYPE"] == "INS") or
                                         (op == "D" and info["SVTYPE"] == "DEL")):
            valid = length == size
            if op == "I":
                valid = valid and querypos + length <= len(contig) and contig[querypos:querypos + length] == fields[4][1:].upper()
            candidates.append(valid)
        if op in "M=XDN":
            refpos += length
        if op in "M=XIS":
            querypos += length
    return querypos == len(contig) and candidates == [True]


def _sequence_alleles(fields, info, reference, repair_padding):
    """Verify REF for every allele and sequence indels without committing failed changes."""
    counts = collections.Counter()
    if not re.fullmatch("[ACGTNacgtn]+", fields[3]):
        raise AlleleValidationError("REF is not a nucleotide sequence")
    fields = fields.copy()
    observed = fields[3].upper()
    expected = reference.fetch(fields[0], int(fields[1]) - 1,
                               int(fields[1]) - 1 + len(observed)).upper()
    if len(expected) != len(observed):
        raise AlleleValidationError("REF span exceeds reference")
    if info["SVTYPE"] not in {"DEL", "INS"} or not re.fullmatch("[ACGTNacgtn]+", fields[4]):
        if observed != expected:
            raise AlleleValidationError("REF disagrees with FASTA for non-sequence allele")
        return fields, counts
    pure_indel = ((info["SVTYPE"] == "INS" and len(observed) == 1 and len(fields[4]) > 1) or
                  (info["SVTYPE"] == "DEL" and len(fields[4]) == 1 and len(observed) > 1))
    if observed != expected:
        if pure_indel and observed[0] != fields[4][0].upper():
            raise AlleleValidationError("REF/ALT have no common padding base")
        n_encoding_only = all(a == b or (a, b) == ("A", "N") for a, b in zip(observed, expected))
        if n_encoding_only:
            if observed[0] != expected[0]:
                if not pure_indel:
                    raise AlleleValidationError("Cannot restore ambiguous REF anchor")
                fields[4] = expected[0] + fields[4][1:]
            fields[3] = expected
            fields[7] += ";HITSV_REF_N_RESTORED"
            counts["reference_n_restored"] += 1
        elif (repair_padding and pure_indel and observed[0] != expected[0] and
              observed[1:] == expected[1:] and _padding_evidence(fields, info)):
            fields[3] = expected
            fields[4] = expected[0] + fields[4][1:]
            fields[7] += ";HITSV_PADDING_REPAIRED;HITSV_ORIGINAL_REF=" + observed
            counts["padding_repaired"] += 1
        else:
            raise AlleleValidationError("REF disagrees with FASTA or padding evidence is insufficient")
    expected_length = len(fields[4]) - len(fields[3])
    try:
        declared_length = int(info["SVLEN"])
    except (ValueError, KeyError) as exc:
        raise AlleleValidationError("Invalid SVLEN") from exc
    if declared_length != expected_length:
        raise AlleleValidationError("SVLEN disagrees with sequence alleles")
    if (info["SVTYPE"] == "DEL" and expected_length >= 0) or (info["SVTYPE"] == "INS" and expected_length <= 0):
        raise AlleleValidationError("SVTYPE disagrees with sequence allele lengths")
    expected_end = int(fields[1]) + len(fields[3]) - 1
    if info.get("END") != str(expected_end):
        items = [item for item in fields[7].split(";") if not item.startswith("END=")]
        items.append("END=" + str(expected_end))
        if info.get("END") is not None:
            items.append("HITSV_ORIGINAL_END=" + info["END"])
        fields[7] = ";".join(items)
        counts["end_corrected"] += 1
    return fields, counts


def normalize(args):
    import pysam
    counts = collections.Counter()
    reasons = collections.Counter()
    repair_padding = getattr(args, "repair_padding", False)
    policy = getattr(args, "invalid_allele_policy", "error")
    rejected_path = getattr(args, "rejected", "") or args.audit + ".rejected.jsonl"
    renamed_info = {"CIGAR": "HITSV_CONTIG_CIGAR", "SOMATIC": "HITSV_SOMATIC"}
    header = None
    with pysam.FastaFile(args.reference) as reference, open(args.output, "w") as output, open(rejected_path, "w") as rejected:
        for path in args.inputs:
            saw_header = False
            with open(path) as handle:
                for line in handle:
                    if line.startswith("##"):
                        if header is None:
                            if line.startswith("##fileDate="):
                                continue
                            for old, new in renamed_info.items():
                                line = line.replace("##INFO=<ID=" + old + ",", "##INFO=<ID=" + new + ",")
                            output.write(line)
                        continue
                    if line.startswith("#CHROM"):
                        fields = line.rstrip().split("\t")
                        if len(fields) != 10:
                            raise ValueError("Expected one VCF sample: " + path)
                        if header is None:
                            output.write('##INFO=<ID=HITSV_ORIGINAL_END,Number=1,Type=Integer,Description="Unmodified upstream END before VCF span correction">\n')
                            output.write('##INFO=<ID=HITSV_REF_N_RESTORED,Number=0,Type=Flag,Description="Reference N bases encoded as A by HitSV restored from FASTA">\n')
                            if repair_padding:
                                output.write('##INFO=<ID=HITSV_PADDING_REPAIRED,Number=0,Type=Flag,Description="Common allele padding corrected using FASTA and verified contig alignment">\n')
                                output.write('##INFO=<ID=HITSV_ORIGINAL_REF,Number=1,Type=String,Description="Upstream REF before padding repair">\n')
                            output.write("\t".join(fields[:9] + [args.sample]) + "\n")
                            header = fields[:9]
                        elif fields[:9] != header:
                            raise ValueError("Inconsistent VCF column header")
                        saw_header = True
                        continue
                    if not saw_header:
                        raise ValueError("Missing VCF header: " + path)
                    fields = line.rstrip().split("\t")
                    if len(fields) != 10 or not fields[1].isdigit():
                        raise ValueError("Malformed VCF record: " + path)
                    counts["raw_records"] += 1
                    info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
                    svtype = info.get("SVTYPE", "")
                    if svtype not in {"DEL", "INS", "DUP", "INV", "BND", "TRA"}:
                        raise ValueError("Unexpected SVTYPE: " + svtype)
                    if fields[4] in {".", ""} or fields[3] in {".", ""}:
                        raise ValueError("Missing REF/ALT")
                    if svtype in {"DEL", "INS", "DUP", "INV"}:
                        if "SVLEN" not in info:
                            raise ValueError("Missing SVLEN")
                        # SRS can emit lt50bp records even with -m 50; enforce the requested cutoff.
                        if abs(int(info["SVLEN"])) < args.min_size:
                            counts["filtered_below_min_size"] += 1
                            continue
                        if "END" in info and int(info["END"]) < int(fields[1]):
                            raise ValueError("END before POS")
                    formats = fields[8].split(":")
                    if "GT" not in formats:
                        raise ValueError("Missing genotype")
                    gt = fields[9].split(":")[formats.index("GT")]
                    if not re.fullmatch(r"[01.][/|][01.]|[01.]", gt):
                        raise ValueError("Invalid genotype: " + gt)
                    counts["eligible_records"] += 1
                    counts["eligible_PASS"] += fields[6] == "PASS"
                    try:
                        fields, changes = _sequence_alleles(fields, info, reference, repair_padding)
                    except AlleleValidationError as exc:
                        if policy == "error":
                            raise AlleleValidationError(str(exc) + ": " + fields[0] + ":" + fields[1]) from exc
                        reason = str(exc)
                        counts["quarantined_records"] += 1
                        counts["quarantined_PASS"] += fields[6] == "PASS"
                        reasons[reason] += 1
                        rejected.write(json.dumps({"source": str(path), "chrom": fields[0], "pos": int(fields[1]),
                            "id": fields[2], "filter": fields[6], "svtype": svtype, "reason": reason,
                            "raw_record": line.rstrip("\n")}) + "\n")
                        continue
                    counts.update(changes)
                    counts["records"] += 1
                    counts[svtype] += 1
                    counts["FILTER_" + fields[6]] += 1
                    counts["sequence_resolved"] += bool(re.fullmatch("[ACGTNacgtn]+", fields[4]))
                    # Avoid reserved VCF INFO names with incompatible upstream meanings.
                    fields[7] = ";".join(renamed_info.get(item.split("=", 1)[0], item.split("=", 1)[0]) +
                                         ("=" + item.split("=", 1)[1] if "=" in item else "")
                                         for item in fields[7].split(";"))
                    # Chromosome-local HitSV IDs can repeat after scatter/gather.
                    fields[2] = "HitSV_" + fields[0] + "_" + str(counts["records"])
                    output.write("\t".join(fields) + "\n")
            if not saw_header:
                raise ValueError("Missing VCF header: " + path)
    report = dict(counts)
    report.update(repair_padding=repair_padding, invalid_allele_policy=policy,
                  rejected_path=str(rejected_path), quarantine_reasons=dict(reasons))
    Path(args.audit).write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("preflight")
    p.add_argument("--bam", action="append", required=True)
    p.add_argument("--index", action="append", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--chromosomes", required=True)
    p.add_argument("--identities", default="{}")
    p.add_argument("--output", required=True)
    p.set_defaults(func=preflight)
    p = sub.add_parser("call")
    for flag in ("binary", "manifest", "chromosome", "reference", "output"):
        p.add_argument("--" + flag, required=True)
    for flag in ("srs", "lrs", "tl", "fa-stat", "bam-stat", "error-profile"):
        p.add_argument("--" + flag, default="")
    p.add_argument("--preset", choices=["HIFI", "ONT_Q20", "ERR_PRONE", "ASM"], default="HIFI")
    p.add_argument("--alignment-preset", choices=["asm5", "asm10", "asm20"], default="asm5")
    p.add_argument("--min-size", type=int, default=50)
    p.set_defaults(func=call)
    p = sub.add_parser("normalize")
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--audit", required=True)
    p.add_argument("--min-size", type=int, default=50)
    p.add_argument("--repair-padding", action="store_true")
    p.add_argument("--invalid-allele-policy", choices=["error", "quarantine"], default="error")
    p.add_argument("--rejected", default="")
    p.set_defaults(func=normalize)
    args = parser.parse_args()
    if args.action == "preflight" and len(args.bam) != len(args.index):
        parser.error("Each BAM requires an explicit index")
    args.func(args)


if __name__ == "__main__":
    main()
