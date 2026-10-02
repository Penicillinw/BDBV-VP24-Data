"""B5 step 4 - pseudobulk the six score genes per (atlas, CL cell type).

Two count conventions are computed side by side, because the earlier probe
(tools/t4b_g9_delta_robustness_20260926.py) used the second one:

  route "raw"      numerator and denominator both from the raw integer counts
                   (`raw/X` when the file has it, else `X` when `X` is itself
                   integral).  CPM = 1e6 * sum(counts) / sum(library size).
                   This is the primary route of B5.
  route "lognorm"  exactly the earlier convention: numerator = sum of the
                   values stored in `X` (which in every one of these files is
                   log1p-normalised, NOT counts) and denominator = the library
                   size.  Kept only as a diagnostic of the earlier probe.

Library size is computed from the same matrices that are scanned, so no atlas
depends on an obs column name.

Outputs: raw/atlas_pseudobulk.tsv (both routes), raw/pseudobulk_qc.json
"""

from __future__ import annotations

import csv
import glob
import json
import math
import os
from collections import defaultdict

import h5py
import numpy as np

from b5_common import MODULE, OUT, RAW, SC, SIX, decode_h5_node

EXTRA = ["KPNA1", "KPNA5", "KPNA6"]
PANEL = SIX + EXTRA
BLOCK = 20000


def var_names(group: h5py.Group) -> list[str]:
    for key in ("feature_name", "_index"):
        if key in group:
            return [str(x) for x in decode_h5_node(group[key])]
    raise KeyError("no feature_name/_index")


def first_positions(genes: list[str], panel: list[str]) -> tuple[dict[int, int], list[str]]:
    """Map var position -> panel column, using the FIRST occurrence of each symbol."""
    pos: dict[int, int] = {}
    duplicated = []
    for k, g in enumerate(panel):
        hits = [i for i, name in enumerate(genes) if name == g]
        if hits:
            pos[hits[0]] = k
            if len(hits) > 1:
                duplicated.append(g)
    return pos, duplicated


def scan_matrix(node, panel_pos: dict[int, int], n_genes: int, want_rowsum: bool):
    """Return (values[n_cells, len(panel)], rowsum[n_cells] or None, integral_fraction)."""
    shape = tuple(node.attrs["shape"]) if isinstance(node, h5py.Group) else node.shape
    n_cells = int(shape[0])
    vals = np.zeros((n_cells, len(panel_pos)), dtype=np.float64)
    rowsum = np.zeros(n_cells, dtype=np.float64) if want_rowsum else None
    if not isinstance(node, h5py.Group):  # dense
        cols = sorted(panel_pos)
        block = node[:, cols]
        for j, c in enumerate(cols):
            vals[:, panel_pos[c]] = block[:, j]
        if want_rowsum:
            rowsum = np.asarray(node[:], dtype=np.float64).sum(axis=1)
        ints = np.mean(np.abs(block - np.round(block)) < 1e-9)
        return vals, rowsum, float(ints)

    data, indices, indptr = node["data"], node["indices"], node["indptr"]
    lut = np.full(n_genes, -1, dtype=np.int16)
    for gidx, col in panel_pos.items():
        lut[gidx] = col
    n_int = 0
    n_nnz = 0
    for i0 in range(0, n_cells, BLOCK):
        i1 = min(n_cells, i0 + BLOCK)
        ip = indptr[i0 : i1 + 1]
        d0, d1 = int(ip[0]), int(ip[-1])
        if d1 <= d0:
            continue
        d = data[d0:d1].astype(np.float64)
        ix = indices[d0:d1]
        n_nnz += d.size
        n_int += int(np.count_nonzero(np.abs(d - np.round(d)) < 1e-9))
        rel = np.asarray(ip - d0, dtype=np.int64)
        if want_rowsum:
            row_ids = np.repeat(np.arange(i1 - i0), np.diff(rel))
            rowsum[i0:i1] = np.bincount(row_ids, weights=d, minlength=i1 - i0)
        col_of = lut[ix]
        m = col_of >= 0
        if m.any():
            rows = np.repeat(np.arange(i1 - i0), np.diff(rel))[m]
            vals[i0 + rows, col_of[m]] = d[m]
    return vals, rowsum, float(n_int / n_nnz) if n_nnz else float("nan")


def aggregate(vals: np.ndarray, rowsum: np.ndarray | None, labels: np.ndarray, panel: list[str]):
    """Sum per (CL term, gene); also sum library size per CL term."""
    uniq, inv = np.unique(labels, return_inverse=True)
    sums = np.zeros((len(uniq), len(panel)), dtype=np.float64)
    counts = np.bincount(inv, minlength=len(uniq))
    for j in range(len(panel)):
        sums[:, j] = np.bincount(inv, weights=vals[:, j], minlength=len(uniq))
    lib = np.bincount(inv, weights=rowsum, minlength=len(uniq)) if rowsum is not None else None
    return {u: {"sums": sums[i], "n": int(counts[i]), "lib": float(lib[i]) if lib is not None else float("nan")} for i, u in enumerate(uniq)}


def components(cpm: dict[str, float]) -> tuple[float, float, float]:
    out = []
    for genes in MODULE.values():
        out.append(sum(math.log10(cpm.get(g, 0.0) + 1.0) for g in genes) / len(genes))
    return tuple(out)  # type: ignore[return-value]


def main() -> None:
    rows = []
    qc = {}
    files = sorted(glob.glob(os.path.join(SC, "*.h5ad")))
    for path in files:
        atlas = os.path.basename(path)
        with h5py.File(path, "r") as handle:
            obs = handle["obs"]
            labels = decode_h5_node(obs["cell_type_ontology_term_id"])

            # count source
            if "raw" in handle and "X" in handle["raw"]:
                raw_node = handle["raw"]["X"]
                raw_genes = var_names(handle["raw"]["var"])
                raw_source = "raw/X"
            else:
                raw_node = handle["X"]
                raw_genes = var_names(handle["var"])
                raw_source = "X"
            x_node = handle["X"]
            x_genes = var_names(handle["var"])

            raw_pos, raw_dup = first_positions(raw_genes, PANEL)
            x_pos, x_dup = first_positions(x_genes, PANEL)

            raw_vals, lib, raw_int = scan_matrix(raw_node, raw_pos, len(raw_genes), want_rowsum=True)
            x_vals, _, x_int = scan_matrix(x_node, x_pos, len(x_genes), want_rowsum=False)

        raw_agg = aggregate(raw_vals, lib, labels, PANEL)
        x_agg = aggregate(x_vals, None, labels, PANEL)

        for cl in sorted(set(labels)):
            r = raw_agg[cl]
            if r["n"] < 1:
                continue
            if r["lib"] > 0:
                cpm = {g: 1e6 * r["sums"][k] / r["lib"] for k, g in enumerate(PANEL)}
            else:
                cpm = {g: 0.0 for g in PANEL}
            c1, c2, c3 = components(cpm)
            n_cells = r["n"]
            lcpm = {g: 1e6 * x_agg[cl]["sums"][k] / r["lib"] if r["lib"] > 0 else 0.0 for k, g in enumerate(PANEL)}
            l1, l2, l3 = components(lcpm)
            rows.append(
                {
                    "atlas": atlas,
                    "cl_id": cl,
                    "n_cells": n_cells,
                    "count_source": raw_source,
                    "lib_sum": int(round(r["lib"])),
                    **{f"cpm_{g}": round(cpm[g], 4) for g in PANEL},
                    "comp_IFN_I": round(c1, 5),
                    "comp_IFN_III": round(c2, 5),
                    "comp_ISG": round(c3, 5),
                    "delta": round((c1 + c2 + c3) / 3.0, 5),
                    **{f"lognorm_{g}": round(lcpm[g], 4) for g in PANEL},
                    "delta_lognorm_route": round((l1 + l2 + l3) / 3.0, 5),
                }
            )
        qc[atlas] = {
            "count_source": raw_source,
            "duplicated_panel_symbols_raw": raw_dup,
            "duplicated_panel_symbols_x": x_dup,
            "n_cells": int(len(labels)),
            "raw_source_integral_fraction": round(raw_int, 6),
            "x_source_integral_fraction": round(x_int, 6),
            "median_library_size": float(np.median(lib)),
            "n_labels": int(len(raw_agg)),
            "median_counts_per_gene_INFAR1": float(np.median(raw_vals[:, PANEL.index("IFNAR1")])),
            "median_counts_per_gene_ISG15": float(np.median(raw_vals[:, PANEL.index("ISG15")])),
        }
        print(f"{atlas}: {len(raw_agg)} labels, source={raw_source}, integral frac raw={raw_int:.3f} x={x_int:.3f}", flush=True)

    os.makedirs(RAW, exist_ok=True)
    with open(os.path.join(RAW, "atlas_pseudobulk.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(RAW, "pseudobulk_qc.json"), "w", encoding="utf-8") as fh:
        json.dump(qc, fh, indent=1)
    print(f"\n{len(rows)} (atlas, CL) pseudobulks over {len(files)} atlases")


if __name__ == "__main__":
    main()
