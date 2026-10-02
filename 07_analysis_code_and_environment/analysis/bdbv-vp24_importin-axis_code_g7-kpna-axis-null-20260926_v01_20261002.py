"""G7 post-hoc robustness: is rho(Delta, KPNA_axis) = -0.270 a real axis separation
or an artefact of per-cell-type overall expression depth?

POST-HOC: NOT part of the G7 pre-registration.  It is a confound check and is
reported as such; it cannot change the frozen quadrant verdict.

Two checks against the FULL local HPA single-cell table
(data/hpa/rna_single_cell_type.tsv, 3,087,080 rows):
  1. partial Spearman rho(Delta, KPNA_axis | depth), depth = per-cell-type median
     log10(nCPM+1) over all genes in the table.
  2. null distribution of |rho| for 2000 random 3-gene axes (same construction),
     so the observed |rho| can be placed on a percentile scale.

Output: analysis/g7_kpna_axis_20260926/g7_null_check.json  (+ stdout report)
"""

import csv
import json
import math
import os
import random
import statistics
from array import array

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")
LAND = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")
OUT = os.path.join(ROOT, "analysis", "g7_kpna_axis_20260926", "g7_null_check.json")

N_NULL = 2000
SEED = 20260926


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def spearman(xs, ys):
    return pearson(ranks(xs), ranks(ys))


def partial_spearman(x, y, z):
    """Spearman rho(x, y | z) via linear residualisation of the ranks."""
    rx, ry, rz = ranks(x), ranks(y), ranks(z)
    n = len(rx)
    mz = sum(rz) / n
    vzz = sum((a - mz) ** 2 for a in rz)
    if vzz == 0:
        return spearman(x, y), 0.0
    def resid(rr):
        mr = sum(rr) / n
        cov = sum((a - mz) * (b - mr) for a, b in zip(rz, rr))
        beta = cov / vzz
        return [b - mr - beta * (a - mz) for a, b in zip(rz, rr)]
    return pearson(resid(rx), resid(ry)), pearson(rz, ranks(x))


def zscores(vals):
    mu = sum(vals) / len(vals)
    sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals))
    if sd == 0:
        return [0.0] * len(vals)
    return [(v - mu) / sd for v in vals]


def main():
    # --- pass 1: cell types -------------------------------------------------
    cell_types, genes = [], set()
    with open(RAW, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        gi = header.index("Gene name")
        ci = header.index("Cell type")
        vi = header.index("nCPM")
        cts = set()
        for row in rd:
            cts.add(row[ci])
            genes.add(row[gi])
    cell_types = sorted(cts)
    ct_idx = {c: i for i, c in enumerate(cell_types)}
    n_ct = len(cell_types)
    print("cell types: %d   distinct genes: %d" % (n_ct, len(genes)))

    # --- pass 2: matrix (memory-lean: 4 bytes per cell) ---------------------
    mat = {g: array("f", bytes(4 * n_ct)) for g in genes}
    depth = [array("f") for _ in range(n_ct)]
    with open(RAW, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for row in rd:
            try:
                v = float(row[vi])
            except ValueError:
                continue
            i = ct_idx[row[ci]]
            mat[row[gi]][i] = v
            depth[i].append(math.log10(v + 1.0))
    del genes

    depth_med = [statistics.median(list(d)) for d in depth]

    grad = {r["cell_type"]: float(r["delta_restriction"]) for r in csv.DictReader(
        open(os.path.join(LAND, "restriction_gradient.tsv"), encoding="utf-8"),
        delimiter="\t")}
    delta = [grad[c] for c in cell_types]

    def axis_of(gene_list):
        cols = []
        for g in gene_list:
            cols.append([math.log10(mat[g][i] + 1.0) for i in range(n_ct)])
        zs = [zscores(col) for col in cols]
        return [sum(zs[k][i] for k in range(len(cols))) / len(cols) for i in range(n_ct)]

    obs_axis = axis_of(["KPNA1", "KPNA5", "KPNA6"])
    rho_obs = spearman(delta, obs_axis)
    rho_partial, rho_depth_delta = partial_spearman(delta, obs_axis, depth_med)
    rho_depth_axis = spearman(depth_med, obs_axis)

    # --- null: random 3-gene axes ------------------------------------------
    rng = random.Random(SEED)
    usable = [g for g in mat if max(mat[g]) > 0 and min(mat[g]) < max(mat[g])]
    print("usable genes for the null: %d" % len(usable))
    null_abs = []
    for _ in range(N_NULL):
        tri = rng.sample(usable, 3)
        null_abs.append(abs(spearman(delta, axis_of(tri))))
    null_abs.sort()
    pct = 100.0 * sum(1 for v in null_abs if v < abs(rho_obs)) / len(null_abs)

    res = {
        "status": "POST-HOC (not pre-registered); cannot change the frozen quadrant verdict",
        "source": "data/hpa/rna_single_cell_type.tsv",
        "n_cell_types": n_ct,
        "n_genes_in_table": len(mat),
        "rho_obs_delta_kpna": rho_obs,
        "rho_depth_vs_delta": rho_depth_delta,
        "rho_depth_vs_kpna_axis": rho_depth_axis,
        "partial_rho_delta_kpna_given_depth": rho_partial,
        "null_random_3gene_axes": {
            "n": len(null_abs),
            "seed": SEED,
            "median_abs_rho": statistics.median(null_abs),
            "p95_abs_rho": null_abs[int(0.95 * (len(null_abs) - 1))],
            "p99_abs_rho": null_abs[int(0.99 * (len(null_abs) - 1))],
            "percentile_of_observed_abs_rho": pct,
        },
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)

    print("")
    print("rho(Delta, KPNA_axis)                 = %+.3f" % rho_obs)
    print("rho(depth, Delta)                      = %+.3f" % rho_depth_delta)
    print("rho(depth, KPNA_axis)                  = %+.3f" % rho_depth_axis)
    print("partial rho(Delta, KPNA_axis | depth)  = %+.3f" % rho_partial)
    print("null |rho| median=%.3f p95=%.3f p99=%.3f -> observed sits at %.1f percentile"
          % (res["null_random_3gene_axes"]["median_abs_rho"],
             res["null_random_3gene_axes"]["p95_abs_rho"],
             res["null_random_3gene_axes"]["p99_abs_rho"], pct))
    print("wrote %s" % os.path.relpath(OUT, ROOT))


if __name__ == "__main__":
    main()
