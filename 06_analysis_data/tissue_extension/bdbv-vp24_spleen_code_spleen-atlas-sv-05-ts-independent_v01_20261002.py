"""SRB-VERIFY step 5: independent recompute of the Tabula Sapiens spleen test.

Task item 2.5 asked for the independent resource to be re-derived rather than
trusted.  This script re-reads analysis/spleen_reference_rebuild_20260926/raw/
TS_spleen.h5ad (read-only) with its own accumulator and re-derives:
  * cell count, number of labels at the 50-cell floor
  * per-label pseudobulk CPM and both delta scales
  * the exact-CL intersection with the HCL spleen labels
  * Spearman, label-permutation p (5,000 draws) and the minimum detectable rho

Writes raw/sv_ts_pseudobulk.tsv, raw/sv_ts_check.json
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections import defaultdict

import h5py
import numpy as np

from sv_03_verify import (SIX, min_detectable_rho, read_categorical, spearman,
                          upper_median)
from sv_04_null_frozen_scale import frozen_delta, native_delta

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUT = os.path.join(ROOT, "analysis", "srb_verify_20260926")
RAW = os.path.join(OUT, "raw")
SRB_RAW = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926", "raw")
TS = os.path.join(SRB_RAW, "TS_spleen.h5ad")
FLOOR = 50
BLOCK = 20000
N_PERM = 5000
SEED = 20260926


def main() -> None:
    stats = {}
    with open(os.path.join(RAW, "sv_hpa_gene_background.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stats[r["gene"]] = (float(r["mean_log10_nCPM"]), float(r["sd_log10_nCPM"]))
    mus = [stats[g][0] for g in SIX]
    sds = [stats[g][1] for g in SIX]

    with h5py.File(TS, "r") as f:
        node = f["raw/X"] if ("raw" in f and "X" in f["raw"]) else f["X"]
        vg = f["raw/var"] if "raw" in f else f["var"]
        genes = [str(x) for x in read_categorical(vg["feature_name"])]
        cl = read_categorical(f["obs"]["cell_type_ontology_term_id"])
        n_cells = len(cl)
        if hasattr(node, "attrs") and "shape" in node.attrs:
            n_genes = int(node.attrs["shape"][1])
        else:
            n_genes = int(node.shape[1])
        pos = {}
        for i, g in enumerate(genes):
            pos.setdefault(g, i)
        want = {pos[g]: k for k, g in enumerate(SIX) if g in pos}
        missing = [g for g in SIX if g not in pos]
        labels, inv = np.unique(cl, return_inverse=True)
        n_lab = len(labels)
        acc = np.zeros((n_lab, len(SIX)), dtype=np.float64)
        lib = np.zeros(n_lab, dtype=np.float64)
        data, indices, indptr = node["data"], node["indices"], node["indptr"]
        intfrac_sum = 0.0
        intfrac_n = 0
        for i0 in range(0, n_cells, BLOCK):
            i1 = min(n_cells, i0 + BLOCK)
            ip = indptr[i0:i1 + 1]
            d0, d1 = int(ip[0]), int(ip[-1])
            if d1 <= d0:
                continue
            d = data[d0:d1].astype(np.float64)
            ix = indices[d0:d1]
            rows = np.repeat(np.arange(i1 - i0, dtype=np.int64), np.diff(ip - d0))
            lib += np.bincount(inv[i0:i1][rows], weights=d, minlength=n_lab)
            if intfrac_n < 2_000_000:
                take = min(200_000, len(d))
                intfrac_sum += float(np.sum(d[:take] == np.round(d[:take])))
                intfrac_n += take
            for gidx, c in want.items():
                m = ix == gidx
                if m.any():
                    acc[:, c] += np.bincount(inv[i0:i1][rows[m]], weights=d[m], minlength=n_lab)
    counts = np.bincount(inv, minlength=n_lab)

    rows_out = []
    for i, lab in enumerate(labels):
        n = int(counts[i])
        if n < FLOOR:
            continue
        cpm = [1e6 * acc[i, k] / lib[i] if lib[i] > 0 else 0.0 for k in range(len(SIX))]
        row = {"label": str(lab), "n_cells": n,
               "delta_raw": native_delta(cpm),
               "delta_transfer": frozen_delta(cpm, mus, sds)}
        for k, g in enumerate(SIX):
            row["cpm_" + g] = cpm[k]
        rows_out.append(row)
    with open(os.path.join(RAW, "sv_ts_pseudobulk.tsv"), "w", encoding="utf-8", newline="") as fh:
        fh.write("label\tn_cells\t" + "\t".join("cpm_" + g for g in SIX)
                 + "\tdelta_raw\tdelta_transfer\n")
        for r in rows_out:
            fh.write(f"{r['label']}\t{r['n_cells']}\t"
                     + "\t".join(f"{r['cpm_' + g]:.4f}" for g in SIX)
                     + f"\t{r['delta_raw']:.5f}\t{r['delta_transfer']:.5f}\n")

    srb = {r["label"]: r for r in csv.DictReader(
        open(os.path.join(SRB_RAW, "ts_spleen_pseudobulk.tsv"), encoding="utf-8"),
        delimiter="\t")}
    mine = {r["label"]: r for r in rows_out}
    worst_n = worst_dr = worst_df = 0.0
    for u, r in srb.items():
        m = mine.get(u)
        if m is None:
            worst_n = float("inf")
            continue
        worst_n = max(worst_n, abs(m["n_cells"] - int(r["n_cells"])))
        worst_dr = max(worst_dr, abs(m["delta_raw"] - float(r["delta_raw"])))
        worst_df = max(worst_df, abs(m["delta_transfer"] - float(r["delta_transfer"])))

    hcl_labels = set()
    with open(os.path.join(RAW, "sv_hcl_spleen_pseudobulk.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["n_cells"]) >= FLOOR:
                hcl_labels.add(r["label"])
    shared = sorted(set(mine) & hcl_labels)
    xs = [mine[u]["delta_transfer"] for u in shared]
    with open(os.path.join(RAW, "sv_hcl_spleen_pseudobulk.tsv"), encoding="utf-8") as fh:
        hcl = {r["label"]: r for r in csv.DictReader(fh, delimiter="\t")}
    ys = [frozen_delta([float(hcl[u][g]) for g in SIX], mus, sds) for u in shared]
    rho = spearman(ys, xs)
    rng = np.random.default_rng(SEED)
    hits = 0
    base = list(ys)
    for _ in range(N_PERM):
        if abs(spearman(list(rng.permutation(base)), xs)) >= abs(rho):
            hits += 1
    p_perm = hits / N_PERM
    mdr = min_detectable_rho(len(shared))

    out = {
        "resource": {"path": os.path.relpath(TS, ROOT), "cells": n_cells,
                     "n_genes": n_genes, "labels_at_floor": len(rows_out),
                     "raw_integral_fraction_sampled": intfrac_sum / intfrac_n if intfrac_n else None,
                     "six_gene_missing": missing,
                     "bytes": os.path.getsize(TS)},
        "vs_SRB": {"max_abs_dN": worst_n, "max_abs_d_delta_raw": worst_dr,
                   "max_abs_d_delta_transfer": worst_df,
                   "verdict": "CONFIRMED" if (worst_n == 0 and worst_dr < 1e-4
                                              and worst_df < 1e-4) else "DIFFERS"},
        "R15": {"n_matched": len(shared), "matched": shared, "spearman": rho,
                "permutation_p": p_perm, "n_perm": N_PERM,
                "min_detectable_rho": mdr,
                "informative": bool(rho >= mdr)},
    }
    with open(os.path.join(RAW, "sv_ts_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
