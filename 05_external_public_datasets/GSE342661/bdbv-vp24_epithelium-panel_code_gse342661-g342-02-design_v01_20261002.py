"""Step 2: bind every TPM-matrix column to its GEO sample condition.

Frozen criterion (written before running):
  * a column is bound only if its name appears verbatim as one of the
    `!Sample_description` values of the series matrix that end in
    `.genes.results`, and the corresponding `!Sample_title` parses into
    (time, treatment, biol rep, route);
  * the run is declared usable only if 180/180 columns bind and the resulting
    grid is complete (6 treatments x 3 times x 2 routes x 5 replicates).

Inputs : raw/GSE342661_counts_tpm.matrix.gz (header only), the series matrix
         already archived by the status check.
Outputs: design.tsv, g342_design.json  (this directory only).
"""

import gzip
import json
import pathlib
import re

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
MATRIX = HERE / "raw" / "GSE342661_counts_tpm.matrix.gz"
SERIES = (
    HERE.parent / "protein_gse342661_20260926" / "raw" / "GSE342661_series_matrix.txt.gz"
)

TITLE_RE = re.compile(
    r"(?P<time>\d+h)\s*,\s*(?P<trt>[^,]+?)\s*,\s*biol rep\s*(?P<rep>\d+)\s*,\s*(?P<route>[ab]+)\s*$"
)


def read_series_field(lines, prefix):
    for line in lines:
        if line.startswith(prefix):
            return [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
    raise SystemExit(f"field {prefix!r} not found in series matrix")


def main() -> None:
    with gzip.open(SERIES, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    gsm = read_series_field(lines, "!Sample_geo_accession")
    titles = read_series_field(lines, "!Sample_title")
    desc_rows = [
        [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
        for line in lines
        if line.startswith("!Sample_description\t")
    ]
    file_rows = [r for r in desc_rows if all(c.endswith(".genes.results") for c in r)]
    if len(file_rows) != 1:
        raise SystemExit(f"expected exactly one .genes.results description row, got {len(file_rows)}")
    matrix_names = file_rows[0]
    assert len(gsm) == len(titles) == len(matrix_names) == 180, (
        len(gsm), len(titles), len(matrix_names)
    )

    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        header = [c.strip().strip('"') for c in fh.readline().rstrip("\n").split("\t")][1:]
    assert len(header) == 180, len(header)

    by_name = {nm: (g, t) for nm, g, t in zip(matrix_names, gsm, titles)}
    rows = []
    unmatched = []
    for col in header:
        hit = by_name.get(col)
        if hit is None:
            unmatched.append(col)
            continue
        g, t = hit
        m = TITLE_RE.match(t)
        if not m:
            unmatched.append(col)
            continue
        rows.append(
            {
                "column": col,
                "gsm": g,
                "title": t,
                "time": m.group("time"),
                "treatment": m.group("trt"),
                "rep": int(m.group("rep")),
                "route": "apical-only" if m.group("route") == "a" else "apical-basal",
            }
        )

    design = pd.DataFrame(rows)
    design.to_csv(HERE / "design.tsv", sep="\t", index=False)

    grid = (
        design.groupby(["treatment", "time", "route"]).size().rename("n").reset_index()
    )
    treatments = sorted(design["treatment"].unique())
    counts = grid["n"].value_counts().to_dict()
    summary = {
        "n_matrix_columns": len(header),
        "n_bound": len(design),
        "n_unmatched": len(unmatched),
        "unmatched": unmatched[:10],
        "treatments": treatments,
        "times": sorted(design["time"].unique()),
        "routes": sorted(design["route"].unique()),
        "grid_cells": int(len(grid)),
        "n_per_grid_cell_counts": {str(k): int(v) for k, v in counts.items()},
        "complete_grid": bool(len(grid) == 6 * 3 * 2 and set(counts) == {5}),
        "duplicate_gsm": int(design["gsm"].duplicated().sum()),
        "duplicate_column": int(design["column"].duplicated().sum()),
    }
    (HERE / "g342_design.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    pivot = (
        design.assign(one=1)
        .pivot_table(index=["time", "route"], columns="treatment", values="one", aggfunc="sum")
        .fillna(0)
        .astype(int)
    )
    print(pivot.to_string())


if __name__ == "__main__":
    main()
