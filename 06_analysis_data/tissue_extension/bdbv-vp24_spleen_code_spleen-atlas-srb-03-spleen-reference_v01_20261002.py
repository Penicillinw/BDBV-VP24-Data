"""SRB step 3 - rebuild the spleen side of the reference score (criteria R5-R9).

Design (frozen before the numbers were read)
-------------------------------------------
The frozen score is a property of the HPA 154-type reference panel.  The spleen
has no member of that panel, so the compartment cannot be placed.  "Rebuilding
the reference score" is taken here to mean: estimate the same score for
spleen-resident cell types from a *spleen-tissue* reference and place it on the
frozen scale.  The 154 frozen values are not touched (R5) - this is an extension
panel, not a re-derivation.

  R6  spleen reference resource = the local Human Cell Landscape file
      (analysis/h2_lowend_atlas_hunt_20260926/raw/atlas_h5ad/2adb1f8a_*.h5ad),
      restricted to the cells whose tissue is `spleen` (rows 191895-207700,
      one contiguous block, 15,806 cells).  The primary label route is
      `cell_type_ontology_term_id` - the same key the frozen pipeline uses, so
      that a CL identifier maps straight onto the frozen label_mapping.tsv - and
      the secondary route is `author_cell_type` (finer, author-supplied).
  R7  pseudobulk convention identical to B5: numerator and denominator from the
      raw integer counts (`raw/X`), CPM = 1e6 * sum(gene counts) / sum(library
      size); labels with fewer than 50 cells are dropped (frozen QC floor).
  R8  two values are reported per label:
        delta_raw       same formula, in the resource's own units
        delta_transfer  the frozen HPA transform applied verbatim, i.e.
                        z_g = (log10(CPM_g + 1) - mu_g) / sd_g with mu/sd from the
                        154 HPA cell types (no fitted parameters)
  R9  anchor gate: the same two values are computed for HCL labels *outside* the
      spleen restriction (whole-atlas pseudobulk, taken from the frozen H2
      product) and compared with the HPA delta of the same CL term.  The transfer
      is declared supported only if Spearman >= 0.5, |slope - 1| <= 0.5 and
      |intercept| <= 0.5 for HPA delta = a + b * delta_transfer.

Writes raw/hcl_spleen_pseudobulk.tsv, raw/srb_anchor_check.tsv,
raw/srb_reference_build.json.
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
H2RAW = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926", "raw")
HCL = os.path.join(H2RAW, "atlas_h5ad", "2adb1f8a_Construction_of_a_human_cell_landscape_at_single.h5ad")

sys.path.insert(0, B5)
from b5_common import MODULE, SIX, decode_h5_node  # noqa: E402

FLOOR = 50
SPLEEN_LOW, SPLEEN_HIGH = 191895, 207700  # inclusive row range, verified one block
BLOCK = 20000


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


def ols(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((v - mx) ** 2 for v in x)
    sxy = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    b = sxy / sxx if sxx else float("nan")
    return b, my - b * mx


def components(cpm, stats=None):
    out = []
    for genes in MODULE.values():
        vals = []
        for g in genes:
            lv = math.log10(cpm.get(g, 0.0) + 1.0)
            if stats is not None:
                lv = (lv - stats[g]["mean_log10_nCPM"]) / stats[g]["sd_log10_nCPM"]
            vals.append(lv)
        out.append(sum(vals) / len(vals))
    return out


def scan_spleen(path):
    """Per-cell 6-gene values and library size for the spleen row block."""
    with h5py.File(path, "r") as f:
        node = f["raw/X"] if ("raw" in f and "X" in f["raw"]) else f["X"]
        varg = f["raw/var"] if node.name.startswith("/raw") else f["var"]
        genes = [str(x) for x in decode_h5_node(varg["feature_name"])]
        pos = {}
        for k, g in enumerate(SIX):
            hit = [i for i, nm in enumerate(genes) if nm == g]
            if hit:
                pos[hit[0]] = k
        missing = [g for g in SIX if g not in genes]
        n_cells = int(node.attrs["shape"][0]) if isinstance(node, h5py.Group) else node.shape[0]
        lo, hi = SPLEEN_LOW, SPLEEN_HIGH + 1
        vals = np.zeros((hi - lo, len(SIX)), dtype=np.float64)
        lib = np.zeros(hi - lo, dtype=np.float64)
        data, indices, indptr = node["data"], node["indices"], node["indptr"]
        lut = np.full(int(node.attrs["shape"][1]), -1, dtype=np.int16)
        for gidx, col in pos.items():
            lut[gidx] = col
        ints, nnz = 0, 0
        for i0 in range(lo, hi, BLOCK):
            i1 = min(hi, i0 + BLOCK)
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
            lib[i0 - lo:i1 - lo] = np.bincount(rows, weights=d, minlength=i1 - i0)
            col_of = lut[ix]
            m = col_of >= 0
            if m.any():
                vals[i0 - lo + rows[m], col_of[m]] = d[m]
    return vals, lib, genes, missing, (n_cells, float(ints / nnz) if nnz else float("nan"))


def main() -> None:
    stats = json.load(open(os.path.join(RAW, "hpa_gene_stats_verified.json"), encoding="utf-8"))["gene_stats"]
    hpa = {r["cell_type"]: r for r in
           csv.DictReader(open(os.path.join(T4, "restriction_gradient.tsv"), encoding="utf-8"), delimiter="\t")}

    with h5py.File(HCL, "r") as f:
        obs = f["obs"]
        tissue = decode_h5_node(obs["tissue"]).astype(str)
        lab_cl = decode_h5_node(obs["cell_type_ontology_term_id"]).astype(str)
        lab_name = decode_h5_node(obs["cell_type"]).astype(str)
        lab_author = decode_h5_node(obs["author_cell_type"]).astype(str)
    sl = np.zeros(len(tissue), dtype=bool)
    sl[SPLEEN_LOW:SPLEEN_HIGH + 1] = True
    assert (tissue[SPLEEN_LOW:SPLEEN_HIGH + 1] == "spleen").all(), "R6 FAIL: row block is not all spleen"
    print(f"R6 PASS: spleen block {SPLEEN_LOW}-{SPLEEN_HIGH} = {sl.sum()} cells, all tissue == spleen")

    vals, lib, genes, missing, (n_cells, intfrac) = scan_spleen(HCL)
    print(f"scan: {vals.shape[0]} cells, raw integral fraction = {intfrac:.3f}, missing genes = {missing}")

    # R6b  invert the frozen label mapping.  Five CL ids are shared by two HPA
    #      types; the tie is broken deterministically by (i) mapping route
    #      exact_label > synonym > curated and (ii) token overlap with the CL
    #      label.  The choice is reported for every anchor.
    ROUTE_RANK = {"exact_label": 0, "synonym": 1, "curated": 2}
    candidates: dict[str, list[dict]] = {}
    for r in csv.DictReader(open(os.path.join(B5, "label_mapping.tsv"), encoding="utf-8"), delimiter="\t"):
        if r["cl_id"]:
            candidates.setdefault(r["cl_id"], []).append(r)

    def overlap(a: str, b: str) -> int:
        stop = {"cell", "cells", "of", "the", "human"}
        ta = {t for t in a.lower().replace("-", " ").split() if t not in stop}
        tb = {t for t in b.lower().replace("-", " ").split() if t not in stop}
        return len(ta & tb)

    cl2hpa, dup = {}, []
    for cl, rows_ in candidates.items():
        rows_.sort(key=lambda r: (ROUTE_RANK.get(r["mapping_route"], 9),
                                  -overlap(r["hpa_cell_type"], r["cl_label"]),
                                  r["hpa_cell_type"]))
        cl2hpa[cl] = rows_[0]["hpa_cell_type"]
        if len(rows_) > 1:
            dup.append({"cl_id": cl, "chosen": rows_[0]["hpa_cell_type"],
                        "rejected": [r["hpa_cell_type"] for r in rows_[1:]],
                        "route": rows_[0]["mapping_route"]})
    print(f"R6b: CL->HPA mapping from label_mapping.tsv: {len(cl2hpa)} ids ({len(dup)} ambiguous, resolved)")

    rows = []
    name_of_cl = {}
    for cid, nm in zip(lab_cl, lab_name):
        name_of_cl.setdefault(cid, nm)

    for route, labels in (("cl", lab_cl), ("author", lab_author)):
        sub = labels[SPLEEN_LOW:SPLEEN_HIGH + 1]
        for lab in sorted(set(sub)):
            m = sub == lab
            n = int(m.sum())
            if n < FLOOR:
                continue
            sums = vals[m].sum(axis=0)
            lsum = float(lib[m].sum())
            if lsum <= 0:
                continue
            cpm = {g: 1e6 * sums[k] / lsum for k, g in enumerate(SIX)}
            d_raw = sum(components(cpm)) / 3.0
            d_tr = sum(components(cpm, stats)) / 3.0
            hpa_type = cl2hpa.get(lab, "")
            rows.append({
                "route": route, "label": lab, "label_name": name_of_cl.get(lab, lab) if route == "cl" else lab,
                "n_cells": n,
                "hpa_member_of_same_cl": hpa_type,
                "hpa_delta_for_same_cl": hpa[hpa_type]["delta_restriction"] if hpa_type in hpa else "",
                **{f"cpm_{g}": round(cpm[g], 3) for g in SIX},
                "delta_raw": round(d_raw, 5), "delta_transfer": round(d_tr, 5),
            })
    with open(os.path.join(RAW, "hcl_spleen_pseudobulk.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"R7/R8: {len(rows)} spleen labels at or above the {FLOOR}-cell floor (both routes)")

    # R9 anchor gate: whole-atlas HCL pseudobulk (frozen H2 product) vs HPA
    anchors = []
    for r in csv.DictReader(open(os.path.join(H2RAW, "new_atlas_pseudobulk.tsv"), encoding="utf-8"), delimiter="\t"):
        if not r["atlas"].startswith("2adb1f8a"):
            continue
        hpa_type = cl2hpa.get(r["cl_id"])
        if hpa_type not in hpa:
            continue
        cpm = {g: float(r[f"cpm_{g}"]) for g in SIX}
        anchors.append({"cl_id": r["cl_id"], "hpa_cell_type": hpa_type, "n_cells_whole_atlas": int(r["n_cells"]),
                        "delta_raw_whole": float(r["delta"]),
                        "delta_transfer_whole": sum(components(cpm, stats)) / 3.0,
                        "hpa_delta": float(hpa[hpa_type]["delta_restriction"])})
    anchors.sort(key=lambda a: a["hpa_cell_type"])
    xr = [a["delta_raw_whole"] for a in anchors]
    xt = [a["delta_transfer_whole"] for a in anchors]
    y = [a["hpa_delta"] for a in anchors]
    rho_raw = spearman(xr, y)
    rho_tr = spearman(xt, y)
    slope, intercept = ols(xt, y)
    gate = (rho_tr >= 0.5) and (abs(slope - 1.0) <= 0.5) and (abs(intercept) <= 0.5)
    print(f"R9 anchors: n = {len(anchors)}; Spearman(raw) = {rho_raw:.3f}; Spearman(transfer) = {rho_tr:.3f}; "
          f"slope = {slope:.3f}, intercept = {intercept:+.3f}; gate = {'SUPPORTED' if gate else 'NOT SUPPORTED'}")
    for a in anchors:
        a["residual_transfer"] = a["delta_transfer_whole"] - a["hpa_delta"]
    with open(os.path.join(RAW, "srb_anchor_check.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(anchors[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(anchors)

    med = sorted(r["delta_transfer"] for r in rows if r["route"] == "cl")
    med_raw = sorted(r["delta_raw"] for r in rows if r["route"] == "cl")
    print(f"spleen (cl route, n = {len(med)}): delta_transfer {med[0]:+.3f} .. {med[-1]:+.3f}; "
          f"spleen median = {med[len(med)//2]:+.3f} (raw-units median {med_raw[len(med_raw)//2]:+.3f})")

    json.dump({
        "R5": "154 frozen HPA deltas untouched; this is an extension panel",
        "R6": {"resource": os.path.relpath(HCL, ROOT), "cells": int(sl.sum()),
               "row_block": [SPLEEN_LOW, SPLEEN_HIGH], "all_spleen": True},
        "R7": {"floor": FLOOR, "count_source": "raw/X", "raw_integral_fraction": intfrac,
               "missing_panel_genes": missing},
        "R8": {"transform": "frozen HPA per-gene z (no fitted parameters)",
               "labels_kept": len(rows)},
        "R9": {"n_anchors": len(anchors), "spearman_raw": rho_raw, "spearman_transfer": rho_tr,
               "slope": slope, "intercept": intercept, "gate": "SUPPORTED" if gate else "NOT SUPPORTED",
               "gate_rule": "rho >= 0.5 and |slope-1| <= 0.5 and |intercept| <= 0.5"},
        "cl_mapping": {"n_ids": len(cl2hpa), "collisions": dup},
    }, open(os.path.join(RAW, "srb_reference_build.json"), "w", encoding="utf-8"), indent=1)
    print("wrote raw/hcl_spleen_pseudobulk.tsv, raw/srb_anchor_check.tsv, raw/srb_reference_build.json")


if __name__ == "__main__":
    main()
