"""A4-filovirus independent recompute, step 1: recover the column -> treatment map
from deposited GEO metadata only.

Frozen criteria (written before any count value was inspected)
-------------------------------------------------------------
This module reads **only** the GEO sample metadata deposited with GSE309699
(``raw/GSE309699_samples_brief.txt`` plus the matrix header of
``raw/GSE309699_raw_counts_All_samples.csv.gz``). It does not read any figure value
and does not read the original implementation's ``column_map.json``.

C1  Every ``^SAMPLE`` block must yield a GSM accession, a
    ``!Sample_characteristics_ch1 = treatment: ...`` line, a
    ``!Sample_characteristics_ch1 = batch: ...`` line, and a
    ``!Sample_description = Column name in "raw_counts_All_samples" - <ID>`` line.
    A block missing any of these is reported, not silently dropped.
C2  The set of column IDs recovered from metadata must equal the set of sample column
    names in the matrix header (order-free equality). Any mismatch is reported.
C3  An arm is a (treatment, batch) pair as recorded in metadata. Arm sizes are read
    off the result and compared with the sizes asserted in the source publication
    (8 alpha / 8 beta / 5 gamma / 8 lambda / 15 untreated).
C4  Nothing in this script infers a treatment from expression values. The map is
    therefore circularity-free with respect to the paper's inhibition ordering.

Writes ``geo_column_map.tsv`` and ``geo_map_summary.json`` into this directory.
"""
from __future__ import annotations

import gzip
import json
import os
import re
from collections import OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"G:\本迪布焦研究"
RAW = os.path.join(ROOT, "analysis", "a2_filovirus_anchor_20260926", "raw")
BRIEF = os.path.join(RAW, "GSE309699_samples_brief.txt")
COUNTS = os.path.join(RAW, "GSE309699_raw_counts_All_samples.csv.gz")

COL_RE = re.compile(r'Column name in "raw_counts_All_samples"\s*-\s*(\S+)')
LIB_RE = re.compile(r"Library name:\s*(.+?)\s*$")
TREAT_RE = re.compile(r"treatment:\s*(.+?)\s*$")
BATCH_RE = re.compile(r"batch:\s*(.+?)\s*$")


def parse_brief():
    """Return one record per ^SAMPLE block, keeping every field we care about."""
    records = []
    cur = None
    with open(BRIEF, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if cur is not None:
                    records.append(cur)
                cur = {"gsm": line.split("=")[-1].strip(), "treatment": None,
                       "batch": None, "library": None, "column": None,
                       "title": None}
                continue
            if cur is None:
                continue
            if line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[1].strip()
                continue
            m = TREAT_RE.search(line)
            if m and line.startswith("!Sample_characteristics_ch1"):
                cur["treatment"] = m.group(1)
                continue
            m = BATCH_RE.search(line)
            if m and line.startswith("!Sample_characteristics_ch1"):
                cur["batch"] = m.group(1)
                continue
            if line.startswith("!Sample_description"):
                m = LIB_RE.search(line)
                if m:
                    cur["library"] = m.group(1)
                    continue
                m = COL_RE.search(line)
                if m:
                    cur["column"] = m.group(1)
                    continue
    if cur is not None:
        records.append(cur)
    return records


def matrix_columns():
    with gzip.open(COUNTS, "rt", encoding="utf-8-sig", errors="replace") as fh:
        header = next(fh).rstrip("\n").split(",")
    return header[2:]


def main():
    recs = parse_brief()
    bad = [r["gsm"] for r in recs
           if not (r["treatment"] and r["batch"] and r["column"])]
    print(f"C1 sample blocks parsed: {len(recs)}; incomplete blocks: {len(bad)} {bad}")

    cols_meta = [r["column"] for r in recs]
    cols_matrix = matrix_columns()
    print(f"C2 metadata columns {len(cols_meta)}, matrix columns {len(cols_matrix)}")
    print(f"C2 set equality: {set(cols_meta) == set(cols_matrix)}")
    only_meta = sorted(set(cols_meta) - set(cols_matrix))
    only_matrix = sorted(set(cols_matrix) - set(cols_meta))
    print(f"C2 only in metadata: {only_meta}; only in matrix: {only_matrix}")

    # arms keyed by treatment string exactly as deposited
    arms = defaultdict(list)
    for r in recs:
        arms[(r["treatment"], r["batch"])].append(r["column"])

    print("\nC3/(arm -> columns) recovered purely from GEO metadata:")
    for (treatment, batch), cols in sorted(arms.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        print(f"  {treatment!r:32s} batch {batch}  n={len(cols):2d}  "
              f"{','.join(sorted(cols))}")

    # letter-prefix grouping (what the count matrix uses)
    letters = defaultdict(list)
    for r in recs:
        letters[r["column"][0]].append(r["column"])
    print("\nletter-prefix group sizes (count-matrix labels):")
    for L in sorted(letters):
        treatments = sorted({r["treatment"] for r in recs if r["column"][0] == L})
        print(f"  {L}: n={len(letters[L])} treatments={treatments}")

    out_rows = sorted(recs, key=lambda r: (r["column"][0], int(r["column"][1:])))
    with open(os.path.join(HERE, "geo_column_map.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        fh.write("column\tletter\tdigit\tgsm\ttreatment\tbatch\tlibrary_name\ttitle\n")
        for r in out_rows:
            fh.write("\t".join([
                r["column"], r["column"][0], r["column"][1:], r["gsm"],
                r["treatment"], r["batch"], r["library"] or "", r["title"] or ""]) + "\n")

    summary = {
        "n_sample_blocks": len(recs),
        "n_incomplete_blocks": len(bad),
        "incomplete_gsms": bad,
        "n_metadata_columns": len(cols_meta),
        "n_matrix_columns": len(cols_matrix),
        "column_set_equality": set(cols_meta) == set(cols_matrix),
        "only_in_metadata": only_meta,
        "only_in_matrix": only_matrix,
        "arms": {f"{t} | batch {b}": sorted(c) for (t, b), c in arms.items()},
        "letter_prefix_sizes": {L: len(v) for L, v in sorted(letters.items())},
        "matrix_column_order": cols_matrix,
    }
    with open(os.path.join(HERE, "geo_map_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print("\nwrote geo_column_map.tsv, geo_map_summary.json")


if __name__ == "__main__":
    main()
