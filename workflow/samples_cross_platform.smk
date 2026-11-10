# Cross-platform (NGS + TGS) sample pairing for integrative analysis
# This module loads a user-maintained mapping table that links the SAME biological
# sample across the two sequencing platforms (e.g. an Illumina short-read sample and
# its matched ONT/PacBio long-read sample). The mapping is REQUIRED because NGS and
# TGS samples are discovered independently and their sample IDs generally differ, so
# the workflow cannot infer the pairing automatically.
#
# CSV format (header required):
#   pair_id,ngs_sample,tgs_sample
#   HG002,HG002_NA24385_son,HG002_ONT
#
# - pair_id:    identifier for the cross-platform result (used in output paths)
# - ngs_sample: sample ID as discovered by the NGS side (matches get_all_sample_ids())
# - tgs_sample: sample ID as discovered by the TGS side (matches get_ont_samples())

import os
import csv

CROSS_PLATFORM_PAIRS_CSV = resolve_input_file_path(
    config.get("samples_pairs_csv", "config/samples_pairs.csv")
)


def _load_cross_platform_pairs(csv_path):
    """Load cross-platform sample pairs from CSV.

    Returns a list of dicts: [{"pair_id", "ngs_sample", "tgs_sample"}, ...].
    Rows with missing fields are skipped with a warning. Duplicate pair_ids are
    rejected because they would map to the same output path.
    """
    pairs = []
    seen_ids = set()
    if not csv_path or not os.path.exists(csv_path):
        return pairs

    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        required = {"pair_id", "ngs_sample", "tgs_sample"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"Cross-platform pairs CSV '{csv_path}' must have header "
                f"'pair_id,ngs_sample,tgs_sample' (found: {reader.fieldnames})"
            )
        for row in reader:
            pair_id = (row.get("pair_id") or "").strip()
            ngs_sample = (row.get("ngs_sample") or "").strip()
            tgs_sample = (row.get("tgs_sample") or "").strip()
            if not pair_id or not ngs_sample or not tgs_sample:
                print(f"Warning: skipping incomplete cross-platform pair row: {row}")
                continue
            if pair_id in seen_ids:
                raise ValueError(
                    f"Duplicate pair_id '{pair_id}' in {csv_path}; pair_id must be unique."
                )
            seen_ids.add(pair_id)
            pairs.append({
                "pair_id": pair_id,
                "ngs_sample": ngs_sample,
                "tgs_sample": tgs_sample,
            })
    return pairs


CROSS_PLATFORM_PAIRS = _load_cross_platform_pairs(CROSS_PLATFORM_PAIRS_CSV)


def _validate_cross_platform_pairs():
    """Return only the pairs whose member samples were discovered on each platform.

    Validation is intentionally soft (warnings, not errors) so that a mapping table
    listing more pairs than the current run provides data for does not abort the whole
    workflow. Only fully-resolvable pairs are exposed to the rules.
    """
    valid = []
    try:
        ngs_ids = set(get_all_sample_ids())
    except NameError:
        ngs_ids = None
    try:
        tgs_ids = set(get_ont_samples())
    except NameError:
        tgs_ids = None

    for pair in CROSS_PLATFORM_PAIRS:
        ok = True
        if ngs_ids is not None and pair["ngs_sample"] not in ngs_ids:
            print(f"Warning: cross-platform pair '{pair['pair_id']}': NGS sample "
                  f"'{pair['ngs_sample']}' not found among discovered NGS samples.")
            ok = False
        if tgs_ids is not None and pair["tgs_sample"] not in tgs_ids:
            print(f"Warning: cross-platform pair '{pair['pair_id']}': TGS sample "
                  f"'{pair['tgs_sample']}' not found among discovered TGS samples.")
            ok = False
        if ok:
            valid.append(pair)
    return valid


def get_cross_platform_pairs():
    """Return the list of validated cross-platform sample pairs."""
    return _validate_cross_platform_pairs()


def get_cross_platform_pair_ids():
    """Return the list of pair_ids for validated cross-platform pairs."""
    return [p["pair_id"] for p in get_cross_platform_pairs()]


def get_cross_platform_pair(pair_id):
    """Return the pair dict for a given pair_id (raises KeyError if unknown)."""
    for pair in CROSS_PLATFORM_PAIRS:
        if pair["pair_id"] == pair_id:
            return pair
    raise KeyError(f"Unknown cross-platform pair_id: {pair_id}")


# Write a small discovery log for debugging (only when a mapping table is present).
if CROSS_PLATFORM_PAIRS:
    try:
        _cp_log = f"{LOGS_DIR}/cross_platform_pairs.csv"
        os.makedirs(os.path.dirname(_cp_log), exist_ok=True)
        with open(_cp_log, "w", newline="") as _fh:
            _w = csv.writer(_fh)
            _w.writerow(["pair_id", "ngs_sample", "tgs_sample"])
            for _p in CROSS_PLATFORM_PAIRS:
                _w.writerow([_p["pair_id"], _p["ngs_sample"], _p["tgs_sample"]])
    except Exception as _e:  # logging must never break workflow parsing
        print(f"Warning: could not write cross-platform pairs log: {_e}")
