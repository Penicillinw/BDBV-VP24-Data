"""SRB step 6 - nulls for the anchor statistic (criterion R18, evidence discipline).

The claim that needs a null is "the spleen panel's ordering agrees with the HPA
panel".  The bridge to HPA is the 20-label anchor regression, whose Spearman is
0.615.  Two nulls are computed:

  (a) label permutation - permute the HPA values across the 20 anchors (5,000
      draws); this asks whether the anchor agreement is an artefact of the anchor
      label set.
  (b) expression-matched random modules - 200 random six-gene modules in which
      each gene is drawn from HPA genes whose mean log10(nCPM+1) across the 154
      reference cell types is within +/-0.10 of its cognate score gene.  Each
      module is pushed through the identical HCL whole-atlas pipeline and its
      anchor Spearman recorded.  This asks whether any six genes with the right
      expression level would have agreed as well.

Criterion R18 (frozen): the anchor agreement is called *specific to the score
genes* only if the observed Spearman exceeds the 95th percentile of the
expression-matched null.

Writes raw/hpa_gene_background_matched.tsv, raw/srb_null.json.
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
import sys
from collections import defaultdict

import h5py
import numpy as np
from scipy import sparse

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926")
RAW = os.path.join(OUT, "raw")
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
H2RAW = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926", "raw")
HCL = os.path.join(H2RAW, "atlas_h5ad", "2adb1f8a_Construction_of_a_human_cell_landscape_at_single.h5ad")
HPA = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")

SIX = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]
MODS = [(0, 1), (2, 3), (4, 5)]
N_MODULES = 200
N_PERM = 5000
SEED = 20260926
BLOCK = 20000
WINDOW = 0.10


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


def gene_background():
    sums: dict[str, list[float]] = defaultdict(list)
    with open(HPA, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            try:
                sums[r["Gene name"]].append(math.log10(float(r["nCPM"]) + 1.0))
            except (ValueError, KeyError):
                continue
    return {g: sum(v) / len(v) for g, v in sums.items() if len(v) == 154}


def delta_of(cpm):
    out = []
    for a, b in MODS:
        out.append((math.log10(cpm[a] + 1.0) + math.log10(cpm[b] + 1.0)) / 2.0)
    return sum(out) / 3.0


def main() -> None:
    sys.path.insert(0, B5)
    from b5_common import decode_h5_node

    bg = gene_background()
    print(f"HPA gene background with complete 154-type coverage: {len(bg)} genes")
    rng = random.Random(SEED)
    targets = [bg[g] for g in SIX]
    pool = {g: bg[g] for g in bg if g not in SIX}
    by_level = sorted(pool.items(), key=lambda kv: kv[1])
    levels = [v for _, v in by_level]
    import bisect

    modules: list[list[str]] = []
    used = set()
    tries = 0
    while len(modules) < N_MODULES and tries < N_MODULES * 200:
        tries += 1
        m = []
        for t in targets:
            lo = bisect.bisect_left(levels, t - WINDOW)
            hi = bisect.bisect_right(levels, t + WINDOW)
            if hi <= lo:
                continue
            for _ in range(50):
                g = by_level[rng.randrange(lo, hi)][0]
                if g not in m and g not in used:
                    m.append(g)
                    break
        if len(m) == 6:
            modules.append(m)
            used.update(m)
    print(f"R18: built {N_MODULES} expression-matched six-gene modules (matched means: "
          f"{', '.join(f'{t:.2f}' for t in targets)})")

    panel = SIX + [g for m in modules for g in m]
    with h5py.File(HCL, "r") as f:
        genes = [str(x) for x in decode_h5_node(f["raw/var"]["feature_name"])]
        lab = decode_h5_node(f["obs"]["cell_type_ontology_term_id"]).astype(str)
        node = f["raw/X"]
        n_cells = int(node.attrs["shape"][0])
        n_genes = int(node.attrs["shape"][1])
        pos = np.full(len(panel), -1, dtype=np.int64)
        first = {}
        for i, g in enumerate(genes):
            first.setdefault(g, i)
        for k, g in enumerate(panel):
            if g in first:
                pos[k] = first[g]
        missing = [panel[k] for k in range(len(panel)) if pos[k] < 0]
        keep = [k for k in range(len(panel)) if pos[k] >= 0]
        cols_idx = pos[keep]
        labels, inv = np.unique(lab, return_inverse=True)
        n_lab = len(labels)
        onehot = sparse.csr_matrix((np.ones(n_cells), (np.arange(n_cells), inv)), shape=(n_cells, n_lab))
        acc = np.zeros((n_lab, len(keep)), dtype=np.float64)
        libsum = np.zeros(n_lab, dtype=np.float64)
        data, indices, indptr = node["data"], node["indices"], node["indptr"]
        lut = np.full(n_genes, -1, dtype=np.int64)
        for j, c in enumerate(cols_idx):
            lut[c] = j
        for i0 in range(0, n_cells, BLOCK):
            i1 = min(n_cells, i0 + BLOCK)
            ip = indptr[i0:i1 + 1]
            d0, d1 = int(ip[0]), int(ip[-1])
            if d1 <= d0:
                continue
            d = data[d0:d1].astype(np.float64)
            ix = indices[d0:d1]
            rel = np.asarray(ip - d0, dtype=np.int64)
            rows = np.repeat(np.arange(i1 - i0), np.diff(rel))
            libsum += np.bincount(inv[i0:i1][rows], weights=d, minlength=n_lab)
            col_of = lut[ix]
            m = col_of >= 0
            if m.any():
                sub = sparse.csr_matrix((d[m], (rows[m], col_of[m])), shape=(i1 - i0, len(keep)))
                acc += np.asarray((onehot[i0:i1].T @ sub).todense())
    print(f"scanned {n_cells} cells x {len(keep)}/{len(panel)} panel genes "
          f"({len(missing)} symbols absent from HCL), {n_lab} atlas labels")

    anchors = list(csv.DictReader(open(os.path.join(RAW, "srb_anchor_check.tsv"), encoding="utf-8"), delimiter="\t"))
    a_hpa = [float(a["hpa_delta"]) for a in anchors]
    lab_index = {l: i for i, l in enumerate(labels.tolist())}
    anchor_rows = [lab_index[a["cl_id"]] for a in anchors]

    col_of_gene = {panel[k]: j for j, k in enumerate(keep)}
    obs_cols = [col_of_gene[g] for g in SIX]
    mod_cols = [[col_of_gene[g] for g in m] if all(g in col_of_gene for g in m) else None for m in modules]

    def anchor_rho(panel_cols):
        xs = []
        for r in anchor_rows:
            cpm = [1e6 * acc[r, j] / libsum[r] for j in panel_cols]
            xs.append(delta_of(cpm))
        return spearman(xs, a_hpa)

    obs = anchor_rho(obs_cols)
    null = [anchor_rho(c) for c in mod_cols if c is not None]
    null_s = sorted(null)
    p95 = null_s[int(0.95 * len(null_s))]
    frac = sum(1 for v in null if v >= obs) / len(null)
    print(f"R18 anchor Spearman: observed {obs:.3f}; expression-matched null median "
          f"{null_s[len(null_s)//2]:.3f}, 95th pct {p95:.3f}, max {null_s[-1]:.3f}; "
          f"fraction >= observed = {frac:.3f}")

    prng = np.random.default_rng(SEED)
    perm = []
    for _ in range(N_PERM):
        y = list(prng.permutation(a_hpa))
        perm.append(spearman([delta_of([1e6 * acc[r, j] / libsum[r] for j in range(6)]) for r in anchor_rows], y))
    perm_p = float(np.mean([abs(v) >= abs(obs) for v in perm]))
    print(f"label-permutation null: p = {perm_p:.4f} ({N_PERM} draws)")

    with open(os.path.join(RAW, "hpa_gene_background_matched.tsv"), "w", encoding="utf-8", newline="") as fh:
        fh.write("score_gene\thpa_mean_log10_nCPM\tmodule_index\tmatched_gene\n")
        for i, m in enumerate(modules):
            for k, g in enumerate(m):
                fh.write(f"{SIX[k]}\t{targets[k]:.6f}\t{i}\t{g}\n")
    json.dump({"R18": {"observed_anchor_spearman": obs, "n_modules": N_MODULES,
                       "window": WINDOW, "null_median": null_s[len(null_s) // 2], "null_p95": p95,
                       "null_max": null_s[-1], "fraction_ge_observed": frac,
                       "specific": bool(obs > p95)},
               "label_permutation": {"n": N_PERM, "p": perm_p},
               "panel": {"genes_requested": len(panel), "genes_used": len(keep), "missing": missing}},
              open(os.path.join(RAW, "srb_null.json"), "w", encoding="utf-8"), indent=1)
    print("wrote raw/hpa_gene_background_matched.tsv, raw/srb_null.json")


if __name__ == "__main__":
    main()
