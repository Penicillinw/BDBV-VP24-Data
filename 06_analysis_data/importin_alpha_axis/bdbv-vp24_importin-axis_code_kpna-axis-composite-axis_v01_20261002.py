"""Composite-axis review (2026-09-30).

Question: does a composite restriction parameter that folds the importin-alpha axis
(KPNA1/KPNA5/KPNA6) into the manuscript's Delta beat Delta alone?

Everything is recomputed here from frozen inputs; nothing under manuscript/, figures/,
analysis/t4_ifn_landscape_20260925/ or analysis/g7_kpna_axis_20260926/ is modified.

Blocks
  A  functional phenotype panel (GSE46599; HIV-1 restriction label)
  B  rank drift of composites over the 154 cell-type reference atlas
  C  effect on the one directional external test (GSE114905, Huh7)
  D  weight identifiability on the 154 atlas (collinearity, PCA, depth, null)

Run:  python analysis/composite_axis_review_20260930/composite_axis_review.py
Seed: 20260926 (fixed; used only for the random 3-gene-axis null in block D)
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import random
import statistics
import sys
import time
from array import array
from itertools import combinations

import numpy as np
from sklearn.linear_model import LogisticRegression

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "composite_axis_review_20260930")
RAW46599 = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926", "raw")
LAND = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")
KPNA = os.path.join(ROOT, "analysis", "g7_kpna_axis_20260926")
HPA_SC = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")

SEED = 20260926
N_NULL = 2000

# --- GSE46599 functional groups (series !Sample_characteristics_ch1) -----------
RESISTANT = ["primary-macrophages", "PMA-THP-1", "THP-1", "U87-MG", "primary-CD4+-T-cells"]
PARTIAL = ["HT1080", "PMA-U937"]
PERMISSIVE = ["CEM", "CEM-SS", "Jurkat", "U937"]

# Delta composition (manuscript) and the importin-alpha axis
COMPONENTS = {
    "typeI": ["IFNAR1", "IFNAR2"],
    "typeIII": ["IFNLR1", "IL10RB"],
    "ISG_prime": ["ISG15", "MX1"],
    "KPNA_genes": ["KPNA1", "KPNA5", "KPNA6"],
}


# =============================================================================
# helpers
# =============================================================================
def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    sd = v.std()
    if sd == 0:
        return np.zeros_like(v)
    return (v - v.mean()) / sd


def auc(scores, labels) -> float:
    """Rank AUC, higher score = positive class (labels True/1)."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float("nan")
    return sum(1.0 if p > q else 0.5 if p == q else 0.0
               for p in pos for q in neg) / (len(pos) * len(neg))


def loo_auc(X: np.ndarray, y: np.ndarray, C: float = 1.0) -> float:
    """Leave-one-out AUC of a logistic fit on X (columns = features)."""
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    y = np.asarray(y, int)
    n = len(y)
    scores = np.zeros(n)
    for i in range(n):
        tr = np.ones(n, bool)
        tr[i] = False
        clf = LogisticRegression(C=C, max_iter=5000).fit(X[tr], y[tr])
        scores[i] = float(clf.decision_function(X[i:i + 1])[0])
    return auc(scores, y)


def exact_perm_null(feature_cols, n_pos, max_iter=1000):
    """Exact label-permutation null of the LOO logistic AUC.

    Each permutation keeps the feature matrix and reassigns which n_pos units are
    positive, so the null is the distribution of the *procedure*'s LOO AUC under
    exchangeable labels.  Returns one array per candidate, aligned by permutation.
    """
    Xs = [np.column_stack(cols) for cols in feature_cols]
    n = Xs[0].shape[0]
    out = [[] for _ in Xs]
    for combo in combinations(range(n), n_pos):
        y = np.zeros(n, int)
        y[list(combo)] = 1
        for k, X in enumerate(Xs):
            scores = np.zeros(n)
            for i in range(n):
                tr = np.ones(n, bool)
                tr[i] = False
                clf = LogisticRegression(C=1.0, max_iter=max_iter).fit(X[tr], y[tr])
                scores[i] = float(clf.decision_function(X[i:i + 1])[0])
            out[k].append(auc(scores, y))
    return [np.array(v) for v in out]


def spearman(xs, ys) -> float:
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

    rx, ry = ranks(list(xs)), ranks(list(ys))
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float("nan")


def pearson(xs, ys) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def partial_spearman(x, y, z):
    """Spearman rho(x,y | z) by residualising ranks on the rank of z."""
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

    rx, ry, rz = ranks(list(x)), ranks(list(y)), ranks(list(z))
    n = len(rx)
    mz = sum(rz) / n
    vzz = sum((a - mz) ** 2 for a in rz)
    if vzz == 0:
        return spearman(x, y), float("nan")

    def resid(rr):
        mr = sum(rr) / n
        beta = sum((a - mz) * (b - mr) for a, b in zip(rz, rr)) / vzz
        return [b - mr - beta * (a - mz) for a, b in zip(rz, rr)]

    return pearson(resid(rx), resid(ry)), spearman(z, x)


# =============================================================================
# A. GSE46599 functional panel
# =============================================================================
def load_gse46599():
    """probe -> symbol (first non-empty symbol in annotation col 3); probe -> vector."""
    probe2sym: dict[str, str] = {}
    with gzip.open(os.path.join(RAW46599, "GPL10558.annot.gz"), "rt",
                   encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                sym_col = next((i for i, h in enumerate(header)
                                if h.strip() in ("Symbol", "GENE_SYMBOL", "ILMN_Gene")), 2)
                continue
            if header is None or line.startswith(("!", "#")) or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > sym_col:
                s = p[sym_col].strip()
                if s and s != "---":
                    probe2sym[p[0].strip()] = s.split("///")[0].strip()

    titles: list[str] = []
    probe_vec: dict[str, list[float]] = {}
    with gzip.open(os.path.join(RAW46599, "GSE46599_series_matrix.txt.gz"), "rt",
                   encoding="utf-8", errors="replace") as fh:
        in_t = False
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_t = True
                next(fh)
            elif in_t:
                if line.startswith("!series_matrix_table_end"):
                    break
                p = line.rstrip("\n").split("\t")
                if len(p) > 1:
                    try:
                        probe_vec[p[0].strip('"')] = [float(x) for x in p[1:]]
                    except ValueError:
                        pass

    gene_vec: dict[str, list[float]] = {}
    gene_probe: dict[str, str] = {}
    for probe, vec in probe_vec.items():
        s = probe2sym.get(probe)
        if s and s not in gene_vec:
            gene_vec[s] = vec
            gene_probe[s] = probe
    return titles, gene_vec, gene_probe


def block_a():
    t0 = time.time()
    titles, genes, gene_probe = load_gse46599()
    cells = []
    for t in titles:
        c = t.rsplit("_", 2)[0]
        if c not in cells:
            cells.append(c)

    def val(syms, cell, cond):
        idx = [i for i, t in enumerate(titles)
               if t.startswith(cell + "_") and f"_{cond}_" in t]
        vals = []
        for g in syms:
            v = genes.get(g)
            if v is not None:
                vals.append(float(np.mean([v[i] for i in idx])))
        return float(np.mean(vals)) if vals else float("nan")

    # raw (un-z) module values per cell type, untreated condition
    raw = {}
    for c in cells:
        rec = {}
        for name, syms in COMPONENTS.items():
            rec[name] = val(syms, c, "None")
        rec["ISG_post"] = val(["ISG15", "MX1"], c, "IFN")
        rec["ISG_induction"] = rec["ISG_post"] - rec["ISG_prime"]
        raw[c] = rec

    probe_map = {}
    for name, syms in COMPONENTS.items():
        probe_map[name] = {g: gene_probe.get(g, None) for g in syms}

    def build(defn_name, groups):
        """groups = (positive, negative) cell lists; z within the union."""
        pos, neg = groups
        sel = pos + neg
        z = {}
        for name in COMPONENTS:
            z[name] = dict(zip(sel, zscore([raw[c][name] for c in sel])))
        z_ISG_post = dict(zip(sel, zscore([raw[c]["ISG_post"] for c in sel])))

        delta = {c: (z["typeI"][c] + z["typeIII"][c] + z["ISG_prime"][c]) / 3.0 for c in sel}
        # KPNA axis = mean of the three per-gene z-scores (the G7 atlas definition)
        kpna_gene_z = {g: dict(zip(sel, zscore([val([g], c, "None") for c in sel])))
                       for g in COMPONENTS["KPNA_genes"]}
        kpna = {c: sum(kpna_gene_z[g][c] for g in COMPONENTS["KPNA_genes"]) / 3.0
                for c in sel}

        d = np.array([delta[c] for c in sel], float)
        k = np.array([kpna[c] for c in sel], float)
        zd, zk = zscore(d), zscore(k)
        y = np.array([c in pos for c in sel], int)

        cand = {}
        for label, score in [
            ("C0_delta", d),
            ("C1_kpna", k),
            ("C2_equal_weight", zd + zk),
            ("C3_weakest_link", np.minimum(zd, zk)),
            ("C4_product", zd * zk),
        ]:
            cand[label] = {
                "features": [score.tolist()],
                "auc_insample": round(auc(score, y), 4),
                "loo_auc": round(loo_auc(score, y), 4),
                "loo_auc_rank_parameter_free": round(auc(score, y), 4),
            }
        # C5: logistic on [z(delta), z(kpna)]
        X5 = np.column_stack([zd, zk])
        clf5 = LogisticRegression(C=1.0, max_iter=5000).fit(X5, y)
        cand["C5_fitted_two_axis"] = {
            "features": X5.T.tolist(),
            "coef": [round(float(c), 4) for c in clf5.coef_[0]],
            "intercept": round(float(clf5.intercept_[0]), 4),
            "auc_insample": round(auc(clf5.decision_function(X5), y), 4),
            "loo_auc": round(loo_auc(X5, y), 4),
        }
        # C6: logistic on the three Delta components
        X6 = np.column_stack([[z["typeI"][c] for c in sel],
                              [z["typeIII"][c] for c in sel],
                              [z["ISG_prime"][c] for c in sel]])
        clf6 = LogisticRegression(C=1.0, max_iter=5000).fit(X6, y)
        cand["C6_fitted_delta_components"] = {
            "features": X6.T.tolist(),
            "coef": [round(float(c), 4) for c in clf6.coef_[0]],
            "intercept": round(float(clf6.intercept_[0]), 4),
            "auc_insample": round(auc(clf6.decision_function(X6), y), 4),
            "loo_auc": round(loo_auc(X6, y), 4),
        }

        # exact label-permutation null of the LOO-AUC of every candidate
        order = ["C0_delta", "C1_kpna", "C2_equal_weight", "C3_weakest_link",
                 "C4_product", "C5_fitted_two_axis", "C6_fitted_delta_components"]
        nulls = exact_perm_null([cand[o]["features"] for o in order], int(y.sum()))
        for o, null in zip(order, nulls):
            obs = cand[o]["loo_auc"]
            cand[o]["loo_auc_perm_null"] = {
                "n_permutations_exact": int(len(null)),
                "median": round(float(np.median(null)), 4),
                "p95": round(float(np.quantile(null, 0.95)), 4),
                "observed": obs,
                "percentile_of_observed": round(
                    100.0 * float((null < obs).sum()) / len(null), 1),
                "exact_p_one_sided_ge_observed": round(
                    float((null >= obs).sum()) / len(null), 4),
            }
        # paired difference vs C0 (same permutation index)
        i0 = order.index("C0_delta")
        for o, null in zip(order, nulls):
            if o == "C0_delta":
                continue
            diff = null - nulls[i0]
            obs_diff = cand[o]["loo_auc"] - cand["C0_delta"]["loo_auc"]
            cand[o]["loo_auc_gain_vs_C0"] = {
                "observed_gain": round(obs_diff, 4),
                "perm_p_one_sided": round(float((diff >= obs_diff).sum()) / len(diff), 4),
            }
        cand["C0_delta"]["loo_auc_gain_vs_C0"] = {
            "observed_gain": 0.0, "perm_p_one_sided": 1.0}
        return {
            "definition": defn_name,
            "n": len(sel),
            "n_positive": int(y.sum()),
            "n_negative": int(len(sel) - y.sum()),
            "units": sel,
            "raw_module_values": {c: {k2: round(v2, 4) for k2, v2 in raw[c].items()}
                                 for c in sel},
            "z_delta_series": {c: round(delta[c], 4) for c in sel},
            "z_kpna_series": {c: round(kpna[c], 4) for c in sel},
            "candidates": cand,
        }

    res = {
        "source": "GSE46599 (GPL10558), untreated baseline modules; phenotype from "
                  "!Sample_characteristics_ch1 'resistance to hiv-1 following ifn treatment'",
        "n_cell_types_in_series": len(cells),
        "cell_types": cells,
        "gene_probe_resolution": probe_map,
        "definition_1_strict": build(
            "1: positive = 5 resistant, negative = 4 permissive (n=9)", (RESISTANT, PERMISSIVE)),
        "definition_2_partial_neg": build(
            "2: positive = 5 resistant, negative = 4 permissive + 2 partial (n=11)",
            (RESISTANT, PERMISSIVE + PARTIAL)),
        "runtime_s": round(time.time() - t0, 1),
    }
    return res


# =============================================================================
# B / C. atlas rank drift
# =============================================================================
def load_atlas():
    grad = {}
    with open(os.path.join(LAND, "restriction_gradient.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            grad[r["cell_type"]] = float(r["delta_restriction"])
    kpna = {}
    with open(os.path.join(KPNA, "kpna_axis.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kpna[r["cell_type"]] = {
                "delta": float(r["delta_restriction"]),
                "delta_rank": int(r["delta_rank"]),
                "KPNA_axis": float(r["KPNA_axis"]),
                "quadrant": r["quadrant"],
            }
    return grad, kpna


def rank_desc(values: dict) -> dict:
    order = sorted(values, key=lambda c: -values[c])
    return {c: i + 1 for i, c in enumerate(order)}


def block_b(weight_sets=None):
    grad, kpna = load_atlas()
    cells = sorted(grad)
    assert set(cells) == set(kpna), "atlas cell-type mismatch"
    d = np.array([grad[c] for c in cells], float)
    k = np.array([kpna[c]["KPNA_axis"] for c in cells], float)
    zd, zk = zscore(d), zscore(k)

    # frozen Delta rank (from the frozen table, not recomputed)
    old_rank = {c: kpna[c]["delta_rank"] for c in cells}
    # sanity: recomputed rank must match the frozen delta_rank
    recomp = rank_desc(dict(zip(cells, d)))
    delta_rank_matches_frozen = all(recomp[c] == old_rank[c] for c in cells)

    comps = {
        "C2_equal_weight": zd + zk,
        "C3_weakest_link": np.minimum(zd, zk),
    }
    for name, (alpha2, beta2) in (weight_sets or {}).items():
        comps[name] = alpha2 * zd + beta2 * zk

    out = {
        "n_cell_types": len(cells),
        "delta_rank_matches_frozen_table": delta_rank_matches_frozen,
        "delta_axis_source": "analysis/t4_ifn_landscape_20260925/restriction_gradient.tsv",
        "kpna_axis_source": "analysis/g7_kpna_axis_20260926/kpna_axis.tsv",
        "composites": {},
    }
    rows = []
    for name, score in comps.items():
        new_rank = rank_desc(dict(zip(cells, score)))
        rho = spearman(d, score)
        top20_old = {c for c in cells if old_rank[c] <= 20}
        top20_new = {c for c in cells if new_rank[c] <= 20}
        top15_old = {c for c in cells if old_rank[c] <= 15}
        top15_new = {c for c in cells if new_rank[c] <= 15}
        shifts = sorted(cells, key=lambda c: -abs(new_rank[c] - old_rank[c]))[:8]
        out["composites"][name] = {
            "spearman_rho_vs_frozen_delta": round(rho, 4),
            "top20_overlap": len(top20_old & top20_new),
            "top15_overlap": len(top15_old & top15_new),
            "top20_gained": sorted(top20_new - top20_old),
            "top20_lost": sorted(top20_old - top20_new),
            "largest_8_shifts": [
                {"cell_type": c, "old_rank": old_rank[c], "new_rank": new_rank[c],
                 "shift": new_rank[c] - old_rank[c],
                 "delta": round(grad[c], 4),
                 "KPNA_axis": round(kpna[c]["KPNA_axis"], 4),
                 "quadrant": kpna[c]["quadrant"]}
                for c in shifts],
        }
        for i, c in enumerate(cells):
            rows.append({
                "composite": name,
                "cell_type": c,
                "delta_restriction": grad[c],
                "delta_rank_frozen": old_rank[c],
                "KPNA_axis": kpna[c]["KPNA_axis"],
                "quadrant": kpna[c]["quadrant"],
                "composite_value": float(score[i]),
                "composite_rank": new_rank[c],
                "rank_shift": new_rank[c] - old_rank[c],
            })
    out["_rows"] = rows
    return out


def block_c(b_result):
    grad, kpna = load_atlas()
    cells = sorted(grad)
    target = "hepatocytes"
    d = np.array([grad[c] for c in cells], float)
    k = np.array([kpna[c]["KPNA_axis"] for c in cells], float)
    zd, zk = zscore(d), zscore(k)
    old_rank = {c: kpna[c]["delta_rank"] for c in cells}
    out = {
        "target_cell_type": target,
        "delta_rank_frozen": old_rank[target],
        "delta_rank_percentile": round(100.0 * old_rank[target] / len(cells), 1),
        "quadrant": kpna[target]["quadrant"],
        "delta_value": round(float(grad[target]), 4),
        "KPNA_axis_value": round(float(kpna[target]["KPNA_axis"]), 4),
        "external_dataset": {
            "accession": "GSE114905",
            "system": "Huh7 (hepatocyte-lineage line), single infections at MOI 0.1",
            "measurement": "summed viral mRNA reads (7 transcripts), 1/2/3 dpi",
            "BDBV_reads": [2344, 42497, 188428],
            "EBOV_reads": [30012, 518797, 559072],
            "BDBV_over_EBOV": [0.0781020924963348, 0.08191450605921005, 0.337037090034915],
        },
        "composites": {},
    }
    for key in b_result["composites"]:
        rec = next(r for r in b_result["_rows"]
                   if r["composite"] == key and r["cell_type"] == target)
        out["composites"][key] = {
            "new_rank": rec["composite_rank"],
            "new_rank_percentile": round(100.0 * rec["composite_rank"] / len(cells), 1),
            "rank_shift": rec["rank_shift"],
            "composite_value": rec["composite_value"],
        }
    out["within_high_delta_quarter"] = {
        "definition": "top 38 cell types by frozen Delta (upper quartile)",
        "delta_rank_cut_in_Delta_alone": old_rank[target],
        "in_top_quartile_by_delta": old_rank[target] <= len(cells) // 4,
    }
    return out


# =============================================================================
# D. identifiability
# =============================================================================
def block_d():
    t0 = time.time()
    grad, kpna = load_atlas()
    cells = sorted(grad)
    d = np.array([grad[c] for c in cells], float)
    k = np.array([kpna[c]["KPNA_axis"] for c in cells], float)
    rho_dk = spearman(d, k)
    r2 = rho_dk ** 2
    vif = 1.0 / (1.0 - r2) if r2 < 1 else float("inf")

    # PCA on [z(Delta), z(KPNA)]
    Z = np.column_stack([zscore(d), zscore(k)])
    cov = np.cov(Z, rowvar=False, bias=True)
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    pc1 = evecs[:, 0]
    var_frac = evals / evals.sum()

    # depth proxy: per-cell-type median log10(nCPM+1) over ALL genes in the HPA table
    cell_idx, genes = {}, set()
    with open(HPA_SC, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        gi = header.index("Gene name")
        ci = header.index("Cell type")
        vi = header.index("nCPM")
        cts = set()
        for row in rd:
            cts.add(row[ci])
            genes.add(row[gi])
    hpa_cells = sorted(cts)
    cell_idx = {c: i for i, c in enumerate(hpa_cells)}
    n_ct = len(hpa_cells)
    mat = {g: array("f", bytes(4 * n_ct)) for g in genes}
    depth = [array("f") for _ in range(n_ct)]
    with open(HPA_SC, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for row in rd:
            try:
                v = float(row[vi])
            except ValueError:
                continue
            i = cell_idx[row[ci]]
            mat[row[gi]][i] = v
            depth[i].append(math.log10(v + 1.0))
    depth_med = {c: statistics.median(list(depth[cell_idx[c]])) for c in cells}
    depth_vec = [depth_med[c] for c in cells]

    rho_depth_delta = spearman(depth_vec, d)
    rho_depth_kpna = spearman(depth_vec, k)
    partial_rho, _ = partial_spearman(d, k, depth_vec)

    # random 3-gene axes null
    def axis_of(gene_list):
        cols = [zscore([math.log10(mat[g][i] + 1.0) for i in range(n_ct)]) for g in gene_list]
        return np.mean(np.column_stack(cols), axis=1)

    delta_by_hpa = np.array([grad[c] for c in hpa_cells], float)
    # sorted so the null gene list (and therefore the whole null) is hash-seed independent
    usable = sorted(g for g in mat if max(mat[g]) > 0 and min(mat[g]) < max(mat[g]))
    rng = random.Random(SEED)
    null_abs = []
    for _ in range(N_NULL):
        tri = rng.sample(usable, 3)
        null_abs.append(abs(spearman(delta_by_hpa, axis_of(tri))))
    null_abs.sort()
    pct = 100.0 * sum(1 for v in null_abs if v < abs(rho_dk)) / len(null_abs)

    out = {
        "n_cell_types": len(cells),
        "rho_delta_kpna": round(rho_dk, 4),
        "vif_each_axis": round(vif, 4),
        "vif_note": "two-regressor case: VIF = 1/(1-rho^2), identical for both axes",
        "r_squared_between_axes": round(r2, 4),
        "pca_pc1_variance_fraction": round(float(var_frac[0]), 4),
        "pca_pc2_variance_fraction": round(float(var_frac[1]), 4),
        "pca_pc1_loadings": {"z_delta": round(float(pc1[0]), 4),
                             "z_kpna": round(float(pc1[1]), 4)},
        "pca_note": "sign of loadings is arbitrary; both entries share the same sign",
        "depth_proxy": "per-cell-type median log10(nCPM+1) over all 20,151 genes in "
                       "data/hpa/rna_single_cell_type.tsv",
        "rho_depth_vs_delta": round(rho_depth_delta, 4),
        "rho_depth_vs_kpna_axis": round(rho_depth_kpna, 4),
        "partial_rho_delta_kpna_given_depth": round(partial_rho, 4),
        "null_random_3gene_axes": {
            "n": len(null_abs),
            "seed": SEED,
            "median_abs_rho": round(statistics.median(null_abs), 4),
            "p95_abs_rho": round(null_abs[int(0.95 * (len(null_abs) - 1))], 4),
            "p99_abs_rho": round(null_abs[int(0.99 * (len(null_abs) - 1))], 4),
            "observed_abs_rho": round(abs(rho_dk), 4),
            "percentile_of_observed": round(pct, 1),
        },
        "hpa_table": {"path": "data/hpa/rna_single_cell_type.tsv",
                      "n_hpa_cell_types": n_ct, "n_genes_in_table": len(mat)},
        "runtime_s": round(time.time() - t0, 1),
    }
    return out


# =============================================================================
def main():
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()

    print("[A] GSE46599 functional panel ...", flush=True)
    a = block_a()

    print("[B] atlas rank drift ...", flush=True)
    # fitted weights from the n=11 (definition 2) logistic, transposed onto the atlas
    coef = a["definition_2_partial_neg"]["candidates"]["C5_fitted_two_axis"]["coef"]
    coef1 = a["definition_1_strict"]["candidates"]["C5_fitted_two_axis"]["coef"]
    b = block_b(weight_sets={
        "C5_fitted_weights": (coef[0], coef[1]),
        "C5_fitted_weights_def1_sensitivity": (coef1[0], coef1[1]),
    })

    print("[C] external direction test ...", flush=True)
    c = block_c(b)

    print("[D] identifiability ...", flush=True)
    d = block_d()

    # --- write outputs -------------------------------------------------------
    rows = b.pop("_rows")
    decision = {
        "generated": "2026-09-30",
        "seed": SEED,
        "A_functional_panel_GSE46599": a,
        "B_atlas_rank_drift": b,
        "C_external_direction_GSE114905": c,
        "D_weight_identifiability": d,
        "runtime_total_s": round(time.time() - t_start, 1),
    }
    with open(os.path.join(OUT, "composite_decision.json"), "w", encoding="utf-8") as fh:
        json.dump(decision, fh, ensure_ascii=False, indent=1)

    # AUC table
    auc_rows = []
    for defkey, label in [("definition_1_strict", "def1_strict_n9"),
                          ("definition_2_partial_neg", "def2_partial_neg_n11")]:
        blk = a[defkey]
        for cand, rec in blk["candidates"].items():
            pn = rec["loo_auc_perm_null"]
            gain = rec["loo_auc_gain_vs_C0"]
            auc_rows.append({
                "negative_class_definition": label,
                "n": blk["n"],
                "n_positive": blk["n_positive"],
                "n_negative": blk["n_negative"],
                "candidate": cand,
                "auc_insample": rec["auc_insample"],
                "loo_auc_logistic": rec["loo_auc"],
                "loo_auc_exact_perm_null_median": pn["median"],
                "loo_auc_exact_perm_null_p95": pn["p95"],
                "loo_auc_null_percentile": pn["percentile_of_observed"],
                "loo_auc_exact_p_one_sided": pn["exact_p_one_sided_ge_observed"],
                "loo_auc_gain_vs_C0": gain["observed_gain"],
                "gain_exact_perm_p": gain["perm_p_one_sided"],
                "fitted_coefficients": (json.dumps(rec["coef"]) if "coef" in rec else ""),
            })
    with open(os.path.join(OUT, "a_auc_table.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(auc_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(auc_rows)

    with open(os.path.join(OUT, "b_rank_drift.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print("\nwrote composite_decision.json, a_auc_table.tsv, b_rank_drift.tsv")
    print("total runtime %.1f s" % (time.time() - t_start))

    # console summary
    print("\n--- A: AUC ---")
    for r in auc_rows:
        print("  %-24s %-22s n=%2d  in=%.3f  loo=%.3f" %
              (r["negative_class_definition"], r["candidate"],
               r["n"], r["auc_insample"], r["loo_auc_logistic"]))
    print("\n--- B: composites ---")
    for name, rec in b["composites"].items():
        print("  %-22s rho=%+.3f top20=%d top15=%d" %
              (name, rec["spearman_rho_vs_frozen_delta"],
               rec["top20_overlap"], rec["top15_overlap"]))
    print("\n--- C: hepatocytes ---")
    print("  delta rank %d/%d, quadrant %s" %
          (c["delta_rank_frozen"], 154, c["quadrant"]))
    for name, rec in c["composites"].items():
        print("  %-22s new_rank=%d shift=%+d" % (name, rec["new_rank"], rec["rank_shift"]))
    print("\n--- D: identifiability ---")
    print("  rho(delta,kpna)=%+.3f VIF=%.3f PC1var=%.3f" %
          (d["rho_delta_kpna"], d["vif_each_axis"], d["pca_pc1_variance_fraction"]))
    print("  null |rho| med=%.3f p95=%.3f observed pct=%.1f" %
          (d["null_random_3gene_axes"]["median_abs_rho"],
           d["null_random_3gene_axes"]["p95_abs_rho"],
           d["null_random_3gene_axes"]["percentile_of_observed"]))


if __name__ == "__main__":
    main()
