"""H2 line: pseudobulk the six Δ-panel genes (plus the KPNA genes) for every
Cell Ontology label in the newly downloaded low-end atlases, using exactly the
conventions of `analysis/b5_harmonised_atlas_20260926/b5_04_pseudobulk.py`
(route "raw": numerator and denominator both from integer counts; CPM = 1e6 *
sum(counts) / sum(library size)).

The output table has the same columns as
`analysis/b5_harmonised_atlas_20260926/raw/atlas_pseudobulk.tsv` so the two can
be concatenated without further translation.

Reads:  raw/atlas_h5ad/*.h5ad   (downloaded by h2_download_targets.py)
Writes: raw/new_atlas_pseudobulk.tsv, raw/new_atlas_pseudobulk_qc.json,
        raw/new_atlas_labels.tsv
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
SC = RAW / "atlas_h5ad"

MODULE = {
    "IFN_I_capacity": ["IFNAR1", "IFNAR2"],
    "IFN_III_capacity": ["IFNLR1", "IL10RB"],
    "ISG_priming": ["ISG15", "MX1"],
}
EXTRA = ["KPNA1", "KPNA5", "KPNA6"]
PANEL = [g for v in MODULE.values() for g in v] + EXTRA
BLOCK = 4000


def decode_h5_node(node) -> np.ndarray:
    """obs/var dataset or categorical group -> object array of str."""
    if isinstance(node, h5py.Group):
        cats = node["categories"][()]
        codes = node["codes"][()]
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in cats]
        return np.array([cats[c] if 0 <= c < len(cats) else "" for c in codes], dtype=object)
    data = node[()]
    if isinstance(data, bytes):
        data = [data]
    return np.array([d.decode() if isinstance(d, bytes) else str(d) for d in data], dtype=object)


def var_names(group: h5py.Group) -> list[str]:
    for key in ("feature_name", "_index", "index"):
        if key in group:
            return [str(x) for x in decode_h5_node(group[key])]
    raise KeyError("no var name column")


def first_positions(genes: list[str], panel: list[str]):
    pos: dict[int, int] = {}
    dup = []
    for k, g in enumerate(panel):
        hits = [i for i, name in enumerate(genes) if name == g]
        if hits:
            pos[hits[0]] = k
            if len(hits) > 1:
                dup.append(g)
    return pos, dup


def matrix_shape(node) -> tuple[int, int]:
    if isinstance(node, h5py.Group):
        return tuple(int(x) for x in node.attrs["shape"])  # type: ignore[return-value]
    return tuple(int(x) for x in node.shape)  # type: ignore[return-value]


def scan_matrix(node, panel_pos: dict[int, int], n_genes: int, want_rowsum: bool):
    n_cells = matrix_shape(node)[0]
    vals = np.zeros((n_cells, len(panel_pos)), dtype=np.float64)
    rowsum = np.zeros(n_cells, dtype=np.float64) if want_rowsum else None

    if not isinstance(node, h5py.Group):  # dense
        cols = sorted(panel_pos)
        block = np.asarray(node[:, cols], dtype=np.float64)
        for j, c in enumerate(cols):
            vals[:, panel_pos[c]] = block[:, j]
        if want_rowsum:
            rowsum = np.asarray(node[:], dtype=np.float64).sum(axis=1)
        ints = float(np.mean(np.abs(block - np.round(block)) < 1e-9))
        return vals, rowsum, ints

    data, indices, indptr = node["data"], node["indices"], node["indptr"]
    lut = np.full(n_genes, -1, dtype=np.int32)
    for gidx, col in panel_pos.items():
        lut[gidx] = col
    n_int = n_nnz = 0
    for i0 in range(0, n_cells, BLOCK):
        i1 = min(n_cells, i0 + BLOCK)
        ip = indptr[i0:i1 + 1]
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


def components(cpm: dict[str, float]) -> tuple[float, float, float]:
    return tuple(  # type: ignore[return-value]
        sum(math.log10(cpm.get(g, 0.0) + 1.0) for g in genes) / len(genes)
        for genes in MODULE.values()
    )


def main() -> None:
    files = sorted(SC.glob("*.h5ad"))
    if not files:
        raise SystemExit("no h5ad files found in %s" % SC)

    rows, qc, label_rows = [], {}, []
    for path in files:
        atlas = path.name
        with h5py.File(path, "r") as h:
            obs = h["obs"]
            cl = decode_h5_node(obs["cell_type_ontology_term_id"])
            name = decode_h5_node(obs["cell_type"]) if "cell_type" in obs else cl

            if "raw" in h and "X" in h["raw"]:
                raw_node = h["raw"]["X"]
                raw_genes = var_names(h["raw"]["var"])
                src = "raw/X"
            else:
                raw_node = h["X"]
                raw_genes = var_names(h["var"])
                src = "X"
            pos, dup = first_positions(raw_genes, PANEL)
            missing = [g for g in PANEL if g not in raw_genes]
            vals, lib, intfrac = scan_matrix(raw_node, pos, len(raw_genes), want_rowsum=True)

        uniq, inv = np.unique(cl, return_inverse=True)
        counts = np.bincount(inv, minlength=len(uniq))
        sums = np.zeros((len(uniq), len(PANEL)), dtype=np.float64)
        for j in range(len(PANEL)):
            sums[:, j] = np.bincount(inv, weights=vals[:, j].astype(np.float64), minlength=len(uniq))
        libsum = np.bincount(inv, weights=lib, minlength=len(uniq))

        label_rows.append({"atlas": atlas, "cl_id": "", "cl_label": "", "n_cells": len(cl),
                           "source": src, "integral_fraction": round(intfrac, 6)})
        for i, u in enumerate(uniq):
            if counts[i] < 1:
                continue
            cpm = {g: (1e6 * sums[i, k] / libsum[i]) if libsum[i] > 0 else 0.0
                   for k, g in enumerate(PANEL)}
            c1, c2, c3 = components(cpm)
            rows.append({
                "atlas": atlas,
                "cl_id": u,
                "n_cells": int(counts[i]),
                "count_source": src,
                "lib_sum": int(round(libsum[i])),
                **{f"cpm_{g}": round(cpm[g], 4) for g in PANEL},
                "comp_IFN_I": round(c1, 5),
                "comp_IFN_III": round(c2, 5),
                "comp_ISG": round(c3, 5),
                "delta": round((c1 + c2 + c3) / 3.0, 5),
            })
            # keep the human-readable label alongside the CL id for the report
            mask = cl == u
            lbl = np.unique(name[mask])
            label_rows.append({
                "atlas": atlas, "cl_id": u, "cl_label": " | ".join(map(str, lbl)),
                "n_cells": int(counts[i]), "source": src,
                "integral_fraction": round(intfrac, 6),
            })

        qc[atlas] = {
            "count_source": src,
            "n_cells": int(len(cl)),
            "n_labels": int(len(uniq)),
            "integral_fraction": round(intfrac, 6),
            "median_library_size": float(np.median(lib)),
            "panel_genes_missing": missing,
            "duplicated_panel_symbols": dup,
        }
        print(f"{atlas}: {len(uniq)} CL labels, source={src}, integral={intfrac:.3f}", flush=True)

    RAW.mkdir(parents=True, exist_ok=True)
    with (RAW / "new_atlas_pseudobulk.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with (RAW / "new_atlas_labels.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(label_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(label_rows)
    (RAW / "new_atlas_pseudobulk_qc.json").write_text(json.dumps(qc, indent=1), encoding="utf-8")
    print(f"\n{len(rows)} (atlas, CL) pseudobulks over {len(files)} new atlases")


if __name__ == "__main__":
    main()
