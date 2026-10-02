"""B5 step 1 - inventory of cell-type labels and gene coverage in the local atlases.

Every cached h5ad in data/scRNAseq_raw/ was written with the CELLxGENE Discover
schema, so `obs/cell_type` and `obs/cell_type_ontology_term_id` are available
for all 20 files. This step reads ONLY the observation metadata and the gene
names, so it stays cheap even for the multi-GB files, and writes the frozen
inventory that the mapping step is built on.

Frozen decisions (made before any score/agreement was computed):
  * harmonisation handle = `cell_type_ontology_term_id` (Cell Ontology, CL);
    the CL label carried in `cell_type` is used as the human-readable name and
    is never used for matching by itself.
  * no cell-type filtering beyond a minimum count that is declared here:
    MIN_CELLS_PSEUDOBULK = 50 cells (same floor the previous g9 probe used).
  * genes used for the score (frozen in the manuscript): IFNAR1, IFNAR2,
    IFNLR1, IL10RB, ISG15, MX1; plus the KPNA probe set for completeness.

Outputs (analysis/b5_harmonised_atlas_20260926/raw/):
  atlas_files.tsv              one row per h5ad: sizes, genes found, total cells
  atlas_celltype_inventory.tsv one row per (atlas, CL id): label, n_cells, tissue
  atlas_cl_union.tsv           one row per CL id: label(s), n_cells, n_atlases
"""

from __future__ import annotations

import csv
import glob
import os

import h5py
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SC = os.path.join(ROOT, "data", "scRNAseq_raw")
OUT = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
RAW = os.path.join(OUT, "raw")

MIN_CELLS_PSEUDOBULK = 50

SCORE_GENES = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]
EXTRA_GENES = ["KPNA1", "KPNA5", "KPNA6", "MX2", "OAS1", "BST2", "IFITM1", "IFITM2"]
WANTED = SCORE_GENES + EXTRA_GENES

OBS_COLS = [
    "cell_type",
    "cell_type_ontology_term_id",
    "tissue",
    "tissue_ontology_term_id",
    "donor_id",
]


def read_obs_column(handle: h5py.File, name: str) -> np.ndarray | None:
    """Return an obs column as an array of str, decoding h5ad categoricals."""
    if name not in handle["obs"]:
        return None
    node = handle["obs"][name]
    if isinstance(node, h5py.Group):  # categorical
        cats = node["categories"][()]
        codes = node["codes"][()]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        out = np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], "")
        return out.astype(object)
    data = node[()]
    if isinstance(data, bytes):
        return np.array([data.decode()])
    return np.array([d.decode() if isinstance(d, bytes) else str(d) for d in data])


def decode_node(node: h5py.Dataset | h5py.Group) -> np.ndarray:
    """Decode an h5ad dataset, or a categorical group, into an array of str."""
    if isinstance(node, h5py.Group):
        cats = node["categories"][()]
        codes = node["codes"][()]
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in cats]
        return np.array([cats[c] if 0 <= c < len(cats) else "" for c in codes], dtype=object)
    data = node[()]
    if isinstance(data, bytes):
        data = [data]
    return np.array([d.decode() if isinstance(d, bytes) else str(d) for d in data], dtype=object)


def read_var_names(handle: h5py.File) -> np.ndarray:
    for key in ("feature_name", "_index"):
        if key in handle["var"]:
            return decode_node(handle["var"][key])
    raise KeyError("no feature_name/_index in var")


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    files = sorted(glob.glob(os.path.join(SC, "*.h5ad")))
    file_rows: list[dict[str, object]] = []
    ct_rows: list[dict[str, object]] = []

    for path in files:
        base = os.path.basename(path)
        with h5py.File(path, "r") as handle:
            x = handle["X"]
            shape = x.shape if hasattr(x, "shape") else tuple(x.attrs["shape"])
            n_obs, n_vars = int(shape[0]), int(shape[1])
            genes = read_var_names(handle)
            present = [g for g in WANTED if g in set(genes.tolist())]
            missing_score = [g for g in SCORE_GENES if g not in present]

            ct = read_obs_column(handle, "cell_type")
            cl = read_obs_column(handle, "cell_type_ontology_term_id")
            tissue = read_obs_column(handle, "tissue")
            if ct is None or cl is None:
                file_rows.append(
                    {
                        "atlas": base,
                        "n_cells": n_obs,
                        "n_genes_total": n_vars,
                        "score_genes_found": "",
                        "score_genes_missing": ",".join(SCORE_GENES),
                        "has_cell_type_ontology": "no",
                        "n_labels": 0,
                    }
                )
                continue

            file_rows.append(
                {
                    "atlas": base,
                    "n_cells": n_obs,
                    "n_genes_total": n_vars,
                    "score_genes_found": ",".join(present),
                    "score_genes_missing": ",".join(missing_score),
                    "has_cell_type_ontology": "yes",
                    "n_labels": int(len({(c, l) for c, l in zip(ct, cl)})),
                }
            )

            # counts per (CL id, CL label), with the dominant tissue for context
            key = np.array([f"{l}\x1f{c}" for c, l in zip(ct, cl)])
            tissues = tissue if tissue is not None else np.array([""] * len(ct))
            agg: dict[str, dict[str, object]] = {}
            for k, tis in zip(key, tissues):
                rec = agg.setdefault(k, {"n": 0, "tissues": {}})
                rec["n"] = int(rec["n"]) + 1
                tmap = rec["tissues"]
                tmap[tis] = tmap.get(tis, 0) + 1
            for k, rec in agg.items():
                cl_id, cl_label = k.split("\x1f")
                tmap = rec["tissues"]
                dom = max(tmap.items(), key=lambda kv: kv[1])[0]
                ct_rows.append(
                    {
                        "atlas": base,
                        "cl_id": cl_id,
                        "cl_label": cl_label,
                        "n_cells": rec["n"],
                        "dominant_tissue": dom,
                        "n_tissues": len(tmap),
                        "passes_min_cells": "yes" if int(rec["n"]) >= MIN_CELLS_PSEUDOBULK else "no",
                    }
                )
        print(f"{base}: {n_obs} cells, {len(present)}/{len(WANTED)} panel genes", flush=True)

    with open(os.path.join(RAW, "atlas_files.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(file_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(file_rows)

    ct_rows.sort(key=lambda r: (r["cl_id"], r["atlas"]))
    with open(os.path.join(RAW, "atlas_celltype_inventory.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ct_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(ct_rows)

    # union across atlases
    union: dict[str, dict[str, object]] = {}
    for r in ct_rows:
        rec = union.setdefault(r["cl_id"], {"labels": set(), "n": 0, "atlases": set()})
        rec["labels"].add(r["cl_label"])
        rec["n"] = int(rec["n"]) + int(r["n_cells"])
        rec["atlases"].add(r["atlas"])
    union_rows = [
        {
            "cl_id": cl_id,
            "cl_labels": " | ".join(sorted(rec["labels"])),
            "n_labels_in_atlases": len(rec["labels"]),
            "total_cells": rec["n"],
            "n_atlases": len(rec["atlases"]),
            "atlases": ",".join(sorted(a.replace(".h5ad", "") for a in rec["atlases"])),
            "passes_min_cells_somewhere": "yes"
            if any(
                r["passes_min_cells"] == "yes" and r["cl_id"] == cl_id for r in ct_rows
            )
            else "no",
        }
        for cl_id, rec in sorted(union.items())
    ]
    union_rows.sort(key=lambda r: -int(r["total_cells"]))
    with open(os.path.join(RAW, "atlas_cl_union.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(union_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(union_rows)

    print(f"\n{len(files)} atlases, {len(union_rows)} distinct CL ids")
    print("wrote:", RAW)


if __name__ == "__main__":
    main()
