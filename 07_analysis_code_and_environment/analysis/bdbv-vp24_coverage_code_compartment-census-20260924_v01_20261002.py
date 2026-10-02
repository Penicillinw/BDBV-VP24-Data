#!/usr/bin/env python
"""Census of every annotated cell type in the 19 local single-cell datasets.

Why: the 56-compartment route table is a curated subset that contains **no
lymphoid compartment**.  The research plan section 7.3 requires the model to
return a negative for lymphocytes, and external test set D4 (Kotliar 2020,
PBMC) is a myeloid-vs-lymphoid contrast, so both checks need lymphoid rows.
This census lists what is actually available before any compartment is added.

Reads only ``obs`` (cell_type, donor_id) from each h5ad - no expression matrix -
so it is cheap even for the 2 GB files.

Outputs (analysis/compartment_census_20260924/)
    celltype_census.tsv        dataset x cell_type x n_cells x n_donors
    dataset_summary.tsv        dataset-level totals and label vocabulary size
    lymphoid_availability.tsv  per dataset: which lymphoid labels exist

Usage
    python tools/compartment_census_20260924.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import h5py
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:  # noqa: BLE001
    pass

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "scRNAseq_raw"
OUT = ROOT / "analysis" / "compartment_census_20260924"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PIPELINE = load_module("entry_route_pipeline", ROOT / "tools" / "entry_route_pipeline_20260919.py")
anchoring = PIPELINE.anchoring

LYMPHOID_HINTS = ("t cell", "b cell", "natural killer", "nk cell", "plasma cell",
                  "lymphocyte", "plasmablast", "mast cell", "tissue-resident memory")


def main() -> int:
    rows: list[dict] = []
    summary: list[dict] = []
    for key, (filename, _) in PIPELINE.DATASETS.items():
        path = RAW / filename
        if not path.exists():
            print(f"  {key}: {filename} missing, skipped")
            continue
        with h5py.File(path, "r") as handle:
            cell_type = anchoring.decode_column(handle["obs"], "cell_type")
            try:
                donor = anchoring.decode_column(handle["obs"], "donor_id")
            except Exception:  # noqa: BLE001
                donor = None
        frame = pd.DataFrame({"cell_type": cell_type})
        if donor is not None:
            frame["donor_id"] = donor[: len(frame)]
        grouped = (
            frame.groupby("cell_type", dropna=False)
            .agg(n_cells=("cell_type", "size"),
                 n_donors=("donor_id", "nunique") if "donor_id" in frame else ("cell_type", "size"))
            .reset_index()
            .sort_values("n_cells", ascending=False)
        )
        for record in grouped.to_dict("records"):
            rows.append({"dataset": key, "file": filename, **record})
        summary.append({
            "dataset": key,
            "file": filename,
            "n_cells": int(len(frame)),
            "n_cell_types": int(grouped.shape[0]),
            "n_donors": int(frame["donor_id"].nunique()) if "donor_id" in frame else None,
        })
        print(f"  {key}: {len(frame)} cells, {grouped.shape[0]} cell types")

    census = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    census.to_csv(OUT / "celltype_census.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(OUT / "dataset_summary.tsv", sep="\t", index=False)

    census["is_lymphoid_hint"] = census["cell_type"].astype(str).str.lower().apply(
        lambda value: any(hint in value for hint in LYMPHOID_HINTS)
    )
    lymphoid = census[census["is_lymphoid_hint"]].copy()
    lymphoid.to_csv(OUT / "lymphoid_availability.tsv", sep="\t", index=False)
    print(f"\nlymphoid-labelled cell types found: {lymphoid.shape[0]} rows "
          f"across {lymphoid['dataset'].nunique()} datasets")
    print(lymphoid[["dataset", "cell_type", "n_cells"]].head(40).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
