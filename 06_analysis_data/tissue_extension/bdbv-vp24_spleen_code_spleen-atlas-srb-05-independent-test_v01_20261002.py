"""SRB step 5 - independent spleen test (criteria R14-R16).

The spleen panel was estimated from HCL.  The independent test resource is
Tabula Sapiens - Spleen (CELLxGENE cee11228-9f0b-4e57-afe2-cfe15ee56312, 70,448
cells, adult donors TSP2/TSP7/TSP14/TSP25, 10x 3' v3 / 10x 5' v2 / Smart-seq2):
a different study, platform and donor set from HCL.

Frozen criteria:
  R14  the test is run on the same CL-identifier key space, floor 50 cells, same
       pseudobulk convention (raw counts), same frozen transform; labels are
       matched by exact Cell Ontology identifier only;
  R15  the statistic is Spearman rho between the two resources over the matched
       splenic labels, with a label-permutation null (5,000 draws) and the
       minimum detectable rho at the observed n reported alongside;
  R16  the test is declared *informative* only if the minimum detectable rho at
       the observed n is <= 0.8; otherwise the outcome is reported as a power
       limit and not as evidence for or against the panel;
  R17  as a platform check, the same transfer offset as in step 3 is computed for
       TS against the HPA values of the matched labels (both resources are then
       compared on the same six cell types).

Writes raw/ts_spleen_pseudobulk.tsv, raw/srb_independent_test.json,
raw/srb_platform_offset.tsv.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys

import h5py
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926")
RAW = os.path.join(OUT, "raw")
T4 = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
TS = os.path.join(RAW, "TS_spleen.h5ad")

FLOOR = 50
BLOCK = 20000
N_PERM = 5000
SEED = 20260926
SIX = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]
MODS = [("IFNAR1", "IFNAR2"), ("IFNLR1", "IL10RB"), ("ISG15", "MX1")]


def components(cpm, stats=None):
    out = []
    for genes in MODS:
        vals = []
        for g in genes:
            lv = math.log10(cpm.get(g, 0.0) + 1.0)
            if stats is not None:
                lv = (lv - stats[g]["mean_log10_nCPM"]) / stats[g]["sd_log10_nCPM"]
            vals.append(lv)
        out.append(sum(vals) / len(vals))
    return out


def spearman(a, b):
    def rank(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    da = math.sqrt(sum((v - ma) ** 2 for v in ra))
    db = math.sqrt(sum((v - mb) ** 2 for v in rb))
    return num / (da * db) if da and db else float("nan")


def min_detectable_rho(n, power=0.8, alpha=0.05):
    if n < 4:
        return float("nan")
    from scipy import stats
    z = (stats.norm.isf(alpha / 2) + stats.norm.isf(1 - power)) / math.sqrt(n - 3)
    return float(math.tanh(z))


def main() -> None:
    sys.path.insert(0, B5)
    from b5_common import decode_h5_node

    stats = json.load(open(os.path.join(RAW, "hpa_gene_stats_verified.json"), encoding="utf-8"))["gene_stats"]
    hpa = {r["cell_type"]: r for r in
           csv.DictReader(open(os.path.join(T4, "restriction_gradient.tsv"), encoding="utf-8"), delimiter="\t")}
    cl2hpa = {}
    ROUTE_RANK = {"exact_label": 0, "synonym": 1, "curated": 2}
    cand: dict[str, list] = {}
    for r in csv.DictReader(open(os.path.join(B5, "label_mapping.tsv"), encoding="utf-8"), delimiter="\t"):
        if r["cl_id"]:
            cand.setdefault(r["cl_id"], []).append(r)
    for cl, rows_ in cand.items():
        rows_.sort(key=lambda r: (ROUTE_RANK.get(r["mapping_route"], 9), r["hpa_cell_type"]))
        cl2hpa[cl] = rows_[0]["hpa_cell_type"]

    # --- scan the independent file
    with h5py.File(TS, "r") as f:
        node = f["raw/X"] if ("raw" in f and "X" in f["raw"]) else f["X"]
        varg = f["raw/var"] if node.name.startswith("/raw") else f["var"]
        genes = [str(x) for x in decode_h5_node(varg["feature_name"])]
        pos = {}
        for k, g in enumerate(SIX):
            hit = [i for i, nm in enumerate(genes) if nm == g]
            if hit:
                pos[hit[0]] = k
        missing = [g for g in SIX if g not in genes]
        lab = decode_h5_node(f["obs"]["cell_type_ontology_term_id"]).astype(str)
        nm = decode_h5_node(f["obs"]["cell_type"]).astype(str)
        n_cells = int(node.attrs["shape"][0]) if isinstance(node, h5py.Group) else node.shape[0]
        n_genes = int(node.attrs["shape"][1]) if isinstance(node, h5py.Group) else node.shape[1]
        vals = np.zeros((n_cells, len(SIX)), dtype=np.float64)
        lib = np.zeros(n_cells, dtype=np.float64)
        data, indices, indptr = node["data"], node["indices"], node["indptr"]
        lut = np.full(n_genes, -1, dtype=np.int16)
        for gidx, col in pos.items():
            lut[gidx] = col
        ints = nnz = 0
        for i0 in range(0, n_cells, BLOCK):
            i1 = min(n_cells, i0 + BLOCK)
            ip = indptr[i0:i1 + 1]
            d0, d1 = int(ip[0]), int(ip[-1])
            if d1 <= d0:
                continue
            d = data[d0:d1].astype(np.float64)
            ix = indices[d0:d1]
            nnz += d.size
            ints += int(np.count_nonzero(np.abs(d - np.round(d)) < 1e-9))
            rel = np.asarray(ip - d0, dtype=np.int64)
            rows = np.repeat(np.arange(i1 - i0), np.diff(rel))
            lib[i0:i1] = np.bincount(rows, weights=d, minlength=i1 - i0)
            col_of = lut[ix]
            m = col_of >= 0
            if m.any():
                vals[i0 + rows[m], col_of[m]] = d[m]
    intfrac = float(ints / nnz) if nnz else float("nan")
    print(f"TS spleen: {n_cells} cells, raw integral fraction {intfrac:.3f}, missing genes {missing}")

    name_of = {}
    for cid, n in zip(lab, nm):
        name_of.setdefault(cid, n)
    rows = []
    for u in sorted(set(lab.tolist())):
        m = lab == u
        if m.sum() < FLOOR:
            continue
        lsum = float(lib[m].sum())
        sums = vals[m].sum(axis=0)
        cpm = {g: 1e6 * sums[k] / lsum for k, g in enumerate(SIX)}
        rows.append({"label": u, "label_name": name_of[u], "n_cells": int(m.sum()),
                     "delta_raw": round(sum(components(cpm)) / 3.0, 5),
                     "delta_transfer": round(sum(components(cpm, stats)) / 3.0, 5)})
    with open(os.path.join(RAW, "ts_spleen_pseudobulk.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"TS labels at or above the {FLOOR}-cell floor: {len(rows)}")

    # --- matched splenic labels
    hcl = {r["label"]: r for r in
           csv.DictReader(open(os.path.join(RAW, "hcl_spleen_pseudobulk.tsv"), encoding="utf-8"), delimiter="\t")
           if r["route"] == "cl"}
    matched = sorted(set(hcl) & {r["label"] for r in rows})
    tsmap = {r["label"]: r for r in rows}
    x = [float(hcl[m]["delta_transfer"]) for m in matched]
    y = [float(tsmap[m]["delta_transfer"]) for m in matched]
    rho = spearman(x, y)
    rng = np.random.default_rng(SEED)
    null = []
    for _ in range(N_PERM):
        p = rng.permutation(y)
        null.append(spearman(x, list(p)))
    pval = float(np.mean([abs(v) >= abs(rho) for v in null]))
    mdr = min_detectable_rho(len(matched))
    print(f"R14/R15: {len(matched)} matched splenic labels; Spearman = {rho:.3f}; "
          f"permutation p = {pval:.4f} ({N_PERM} draws); min detectable rho at n = {len(matched)} is {mdr:.3f}")
    informative = (not math.isnan(mdr)) and mdr <= 0.8
    print(f"R16: test {'INFORMATIVE' if informative else 'NOT INFORMATIVE (power limit)'}")
    for m in matched:
        print(f"    {m:12s} {hcl[m]['label_name']:26s} HCL {float(hcl[m]['delta_transfer']):+.3f}  "
              f"TS {float(tsmap[m]['delta_transfer']):+.3f}")

    # --- R17 platform offset on the six labels common to HCL, TS and HPA
    plat = []
    for m in matched:
        if m not in cl2hpa:
            continue
        ht = cl2hpa[m]
        if ht not in hpa:
            continue
        plat.append({"cl_id": m, "hpa_cell_type": ht, "hpa_delta": float(hpa[ht]["delta_restriction"]),
                     "hcl_transfer": float(hcl[m]["delta_transfer"]),
                     "ts_transfer": float(tsmap[m]["delta_transfer"]),
                     "hcl_offset": float(hcl[m]["delta_transfer"]) - float(hpa[ht]["delta_restriction"]),
                     "ts_offset": float(tsmap[m]["delta_transfer"]) - float(hpa[ht]["delta_restriction"])})
    with open(os.path.join(RAW, "srb_platform_offset.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(plat[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(plat)
    if plat:
        hcl_off = sorted(p["hcl_offset"] for p in plat)
        ts_off = sorted(p["ts_offset"] for p in plat)
        print(f"R17 platform offset on {len(plat)} shared labels: HCL mean "
              f"{sum(hcl_off)/len(hcl_off):+.3f} (range {hcl_off[0]:+.3f}..{hcl_off[-1]:+.3f}); TS mean "
              f"{sum(ts_off)/len(ts_off):+.3f} (range {ts_off[0]:+.3f}..{ts_off[-1]:+.3f})")

    json.dump({
        "resource": {"path": os.path.relpath(TS, ROOT), "cells": n_cells, "raw_integral_fraction": intfrac,
                     "missing_panel_genes": missing, "n_labels_ge_floor": len(rows),
                     "dataset_id": "cee11228-9f0b-4e57-afe2-cfe15ee56312", "title": "Tabula Sapiens - Spleen"},
        "R15": {"n_matched": len(matched), "matched": matched, "spearman": rho, "permutation_p": pval,
                "n_perm": N_PERM, "min_detectable_rho": mdr, "informative": bool(informative)},
        "R17": {"n_shared_with_hpa": len(plat),
                "hcl_offset_mean": float(sum(p["hcl_offset"] for p in plat) / len(plat)) if plat else None,
                "ts_offset_mean": float(sum(p["ts_offset"] for p in plat) / len(plat)) if plat else None},
        "ts_finest_labels": [{"label": r["label"], "name": r["label_name"], "n_cells": r["n_cells"]} for r in rows],
    }, open(os.path.join(RAW, "srb_independent_test.json"), "w", encoding="utf-8"), indent=1)
    print("wrote raw/ts_spleen_pseudobulk.tsv, raw/srb_independent_test.json, raw/srb_platform_offset.tsv")


if __name__ == "__main__":
    main()
