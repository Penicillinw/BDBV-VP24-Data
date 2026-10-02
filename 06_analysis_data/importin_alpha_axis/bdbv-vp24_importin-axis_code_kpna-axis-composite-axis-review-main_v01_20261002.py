"""Main-agent independent review of "should Delta and the importin-alpha axis be fused
into one composite cell-restriction parameter?".

Deliberately uses a different code path from tools/a1c_components_and_kpna_20260926.py:
  * axes are rebuilt from the wide HPA table before anything else (self-check vs frozen files)
  * AUC is rank-based (Mann-Whitney with mid-ranks), not pairwise counting
  * probe aggregation is reported under two rules (first probe / median across probes)
  * LOO-AUC is used for every fitted form; fixed-weight forms have LOO-AUC == in-sample AUC

Read-only with respect to every frozen artefact. Writes only under
analysis/composite_axis_review_main_20260930/.
"""

from __future__ import annotations

import gzip
import json
import os

import numpy as np
import pandas as pd

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "composite_axis_review_main_20260930")
os.makedirs(OUT, exist_ok=True)

WIDE = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                    "hpa_ifn_landscape_wide.tsv")
FROZEN_DELTA = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                            "restriction_gradient.tsv")
FROZEN_KPNA = os.path.join(ROOT, "analysis", "g7_kpna_axis_20260926", "kpna_axis.tsv")
GSE46599_RAW = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926", "raw")

COMPONENTS = {
    "typeI": ["IFNAR1", "IFNAR2"],
    "typeIII": ["IFNLR1", "IL10RB"],
    "ISG_prime": ["ISG15", "MX1"],
}
KPNA_GENES = ["KPNA1", "KPNA5", "KPNA6"]

# HIV-1 restriction phenotype from the GSE46599 series summary (Goujon & Malim 2010;
# Goujon & Schaller 2013), transcribed unchanged from the deposited summary text.
PHENOTYPE = {
    "primary-macrophages": 1.0, "PMA-THP-1": 1.0, "THP-1": 1.0, "U87-MG": 1.0,
    "primary-CD4+-T-cells": 0.75, "HT1080": 0.5, "PMA-U937": 0.5,
    "CEM": 0.0, "CEM-SS": 0.0, "Jurkat": 0.0, "U937": 0.0,
}


# --------------------------------------------------------------------------- helpers
def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    sd = v.std(ddof=0)
    return (v - v.mean()) / (sd if sd > 1e-12 else 1.0)


def auc_midrank(scores, labels) -> float:
    """Mann-Whitney AUC with mid-ranks for ties."""
    s = np.asarray(scores, float)
    y = np.asarray(labels, bool)
    order = s.argsort(kind="mergesort")
    ranks = np.empty(len(s), float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    # mid-rank pass for ties
    uniq, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    for j, u in enumerate(uniq):
        if cnt[j] > 1:
            ranks[inv == j] = ranks[inv == j].mean()
    n1 = int(y.sum())
    n0 = int((~y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r1 = ranks[y].sum()
    return float((r1 - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def loo_logistic_auc(X: np.ndarray, y: np.ndarray) -> float:
    from sklearn.linear_model import LogisticRegression

    n = len(y)
    out = np.zeros(n)
    for i in range(n):
        keep = np.ones(n, bool)
        keep[i] = False
        if len(np.unique(y[keep])) < 2:
            out[i] = 0.0
            continue
        clf = LogisticRegression(C=1.0, max_iter=5000).fit(X[keep], y[keep])
        out[i] = float(clf.decision_function(X[i:i + 1])[0])
    return auc_midrank(out, y)


def spearman(a, b) -> float:
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)


# ------------------------------------------------------------------ rebuild the axes
def rebuild_atlas_axes() -> dict:
    wide = pd.read_csv(WIDE, sep="\t")
    ct = wide["cell_type"].tolist()
    d = {"cell_type": ct}
    zcols = {}
    for name, genes in COMPONENTS.items():
        zg = [zscore(np.log10(wide[g].to_numpy(float) + 1.0)) for g in genes]
        d[name] = np.mean(zg, axis=0)
        zcols[name] = d[name]
    d["delta_rebuilt"] = np.mean([zcols[k] for k in COMPONENTS], axis=0)
    zk = [zscore(np.log10(wide[g].to_numpy(float) + 1.0)) for g in KPNA_GENES]
    d["kpna_axis_rebuilt"] = np.mean(zk, axis=0)
    for g in KPNA_GENES:
        d["z_" + g] = zscore(np.log10(wide[g].to_numpy(float) + 1.0))

    frozen_d = pd.read_csv(FROZEN_DELTA, sep="\t")[["cell_type", "delta_restriction", "rank"]]
    frozen_k = pd.read_csv(FROZEN_KPNA, sep="\t")[["cell_type", "KPNA_axis", "quadrant"]]
    df = pd.DataFrame(d).merge(frozen_d, on="cell_type", how="left").merge(
        frozen_k, on="cell_type", how="inner")
    df["z_delta"] = zscore(df["delta_rebuilt"].to_numpy())
    df["z_kpna"] = zscore(df["kpna_axis_rebuilt"].to_numpy())
    return df


def part_b_and_c(df: pd.DataFrame) -> dict:
    df = df.copy()
    df["C0_delta"] = df["z_delta"]
    df["C2_equal"] = df["z_delta"] + df["z_kpna"]
    df["C3_min"] = np.minimum(df["z_delta"], df["z_kpna"])
    df["C4_product"] = df["z_delta"] * df["z_kpna"]

    res = {}
    for col in ["C0_delta", "C2_equal", "C3_min", "C4_product"]:
        df["rank_" + col] = df[col].rank(ascending=False, method="min").astype(int)

    top20 = set(df.nsmallest(20, "rank_C0_delta")["cell_type"])
    top15 = set(df.nsmallest(15, "rank_C0_delta")["cell_type"])
    for col in ["C2_equal", "C3_min", "C4_product"]:
        res[col] = {
            "spearman_vs_frozen_delta": round(spearman(df["delta_restriction"], df[col]), 4),
            "top20_overlap": len(top20 & set(df.nsmallest(20, "rank_" + col)["cell_type"])),
            "top15_overlap": len(top15 & set(df.nsmallest(15, "rank_" + col)["cell_type"])),
        }
        sh = df.copy()
        sh["shift"] = sh["rank_" + col] - sh["rank_C0_delta"]
        res[col]["largest_movers"] = [
            {"cell_type": r.cell_type, "delta_rank": int(r.rank_C0_delta),
             "composite_rank": int(getattr(r, "rank_" + col))}
            for r in sh.reindex(sh["shift"].abs().sort_values(ascending=False).index).head(8).itertuples()
        ]

    # Part C: the hepatocyte lineage used by the one multiplicity-matched dataset
    hep = df[df["cell_type"].str.contains("hepatocyte", case=False)]
    median_composite = float(df["C2_equal"].median())
    c = {}
    for r in hep.itertuples():
        c[r.cell_type] = {
            "delta_rank": int(r.rank_C0_delta),
            "delta": round(float(r.delta_restriction), 4),
            "kpna_axis": round(float(r.KPNA_axis), 4),
            "quadrant": r.quadrant,
            "composite_C2_rank": int(getattr(r, "rank_C2_equal")),
            "above_composite_median": bool(r.C2_equal > median_composite),
        }
    res["_partC_hepatocyte_lineage"] = c

    # named cell types used by the manuscript's falsification design
    named = ["urothelial cells", "enterocytes", "colonocytes", "hepatocytes",
             "respiratory ciliated cells", "pdcs", "monocytes", "macrophages",
             "microglia", "paneth cells", "extravillous trophoblasts"]
    named_out = {}
    for name in named:
        sub = df[df["cell_type"] == name]
        if sub.empty:
            named_out[name] = None
            continue
        r = sub.iloc[0]
        named_out[name] = {
            "delta": round(float(r.delta_restriction), 3),
            "delta_rank": int(r.rank_C0_delta),
            "kpna_axis": round(float(r.KPNA_axis), 3),
            "quadrant": r.quadrant,
            "rank_C2_equal": int(r.rank_C2_equal),
            "rank_C3_min": int(r.rank_C3_min),
            "rank_C4_product": int(r.rank_C4_product),
        }
    res["_named_cell_types"] = named_out
    res["_median_composite_C2"] = round(median_composite, 4)

    # store the drift table for the report
    df[["cell_type", "delta_restriction", "rank", "KPNA_axis", "quadrant",
        "C0_delta", "C2_equal", "C3_min", "C4_product",
        "rank_C0_delta", "rank_C2_equal", "rank_C3_min", "rank_C4_product"]].to_csv(
        os.path.join(OUT, "b_rank_drift_main.tsv"), sep="\t", index=False)
    return res


# ----------------------------------------------------------------- functional panel
def load_probe_to_symbol() -> dict:
    p2s = {}
    with gzip.open(os.path.join(GSE46599_RAW, "GPL10558.annot.gz"), "rt",
                   encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                continue
            if header is None or not line.strip() or line.startswith(("!", "#")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > 2:
                sym = p[2].strip().split("///")[0].strip()
                if sym and sym != "---":
                    p2s.setdefault(p[0].strip(), sym)
    return p2s


def load_gse46599():
    p2s = load_probe_to_symbol()
    titles, rows = [], {}
    with gzip.open(os.path.join(GSE46599_RAW, "GSE46599_series_matrix.txt.gz"), "rt",
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
                        rows[p[0].strip('"')] = np.array([float(x) for x in p[1:]])
                    except ValueError:
                        pass
    return titles, rows, p2s


def gene_matrix(titles, rows, p2s, probe_rule: str) -> dict:
    """gene -> {cell -> {cond -> value}} averaged over replicates within cell/cond."""
    cells, conds = [], {}
    for i, t in enumerate(titles):
        head = t.rsplit("_", 2)[0]
        cond = "IFN" if "_IFN_" in t else "None"
        if head not in cells:
            cells.append(head)
        conds[i] = (head, cond)

    by_sym = {}
    for probe, vec in rows.items():
        sym = p2s.get(probe)
        if sym:
            by_sym.setdefault(sym, []).append(vec)

    wanted = [g for gs in COMPONENTS.values() for g in gs] + KPNA_GENES
    out = {}
    for sym in wanted:
        vecs = by_sym.get(sym)
        if not vecs:
            out[sym] = None
            continue
        if probe_rule == "first":
            vec = vecs[0]
        else:
            vec = np.median(np.vstack(vecs), axis=0)
        per = {}
        for cell in cells:
            per[cell] = {}
            for cond in ("None", "IFN"):
                idx = [i for i, (c, k) in conds.items() if c == cell and k == cond]
                per[cell][cond] = float(np.mean(vec[idx])) if idx else float("nan")
        out[sym] = per
    return out, cells


def part_a(probe_rule: str) -> dict:
    titles, rows, p2s = load_gse46599()
    genes, cells = gene_matrix(titles, rows, p2s, probe_rule)
    cells = [c for c in cells if c in PHENOTYPE]

    def comp(names, cond):
        vals = [np.array([genes[g][c][cond] for c in cells]) for g in names if genes.get(g)]
        return np.mean(vals, axis=0)

    typeI = comp(COMPONENTS["typeI"], "None")
    typeIII = comp(COMPONENTS["typeIII"], "None")
    isg = comp(COMPONENTS["ISG_prime"], "None")
    kpna = comp(KPNA_GENES, "None")

    zI, zIII, zS, zK = zscore(typeI), zscore(typeIII), zscore(isg), zscore(kpna)
    zD = (zI + zIII + zS) / 3.0

    pheno = np.array([PHENOTYPE[c] for c in cells])
    y_all = pheno >= 0.75                       # 5 positive, 6 negative (partials counted negative)
    y_strict_pos = pheno >= 0.75
    y_strict_neg = pheno == 0.0
    keep_strict = y_strict_pos | y_strict_neg

    forms = {
        "C0_delta_manuscript": zD,
        "C1_kpna_axis": zK,
        "C1b_kpna_negated": -zK,
        "C2_equal_sum": zD + zK,
        "C3_weakest_link_min": np.minimum(zD, zK),
        "C4_product": zD * zK,
        "typeI_only": zI,
        "typeIII_only": zIII,
        "ISG_priming_only": zS,
    }
    out = {"probe_rule": probe_rule, "cells": cells, "n_cells": len(cells),
           "phenotype": {c: PHENOTYPE[c] for c in cells},
           "negative_group_definitions": {
               "all_11_partials_as_negative": "positive = phenotype>=0.75 (5); negative = the other 6",
               "strict_partials_excluded": "positive = 5; negative = the 4 phenotype==0",
           },
           "forms": {}, "fitted_two_factor": {}, "fitted_three_component": {}}

    for name, s in forms.items():
        out["forms"][name] = {
            "auc_all11": round(auc_midrank(s, y_all), 3),
            "auc_strict8": round(auc_midrank(s[keep_strict], y_all[keep_strict]), 3),
            "loo_auc_all11": round(auc_midrank(s, y_all), 3),  # fixed weights: no fitting
        }

    X2 = np.column_stack([zD, zK])
    from sklearn.linear_model import LogisticRegression
    clf2 = LogisticRegression(C=1.0, max_iter=5000).fit(X2, y_all)
    out["fitted_two_factor"] = {
        "coefficients_on_zD_and_zKPNA": [round(float(c), 3) for c in clf2.coef_[0]],
        "auc_all11_in_sample": round(auc_midrank(clf2.decision_function(X2), y_all), 3),
        "loo_auc_all11": round(loo_logistic_auc(X2, y_all), 3),
    }
    X3 = np.column_stack([zI, zIII, zS])
    clf3 = LogisticRegression(C=1.0, max_iter=5000).fit(X3, y_all)
    out["fitted_three_component"] = {
        "coefficients_on_typeI_typeIII_ISGprime": [round(float(c), 3) for c in clf3.coef_[0]],
        "auc_all11_in_sample": round(auc_midrank(clf3.decision_function(X3), y_all), 3),
        "loo_auc_all11": round(loo_logistic_auc(X3, y_all), 3),
    }
    return out


def main() -> None:
    atlas = rebuild_atlas_axes()
    self_check = {
        "n_cell_types": int(len(atlas)),
        "max_abs_diff_delta_vs_frozen": float(
            (atlas["delta_rebuilt"] - atlas["delta_restriction"]).abs().max()),
        "max_abs_diff_kpna_axis_vs_frozen": float(
            (atlas["kpna_axis_rebuilt"] - atlas["KPNA_axis"]).abs().max()),
        "spearman_delta_vs_frozen": round(
            spearman(atlas["delta_rebuilt"], atlas["delta_restriction"]), 4),
        "spearman_kpna_vs_frozen": round(
            spearman(atlas["kpna_axis_rebuilt"], atlas["KPNA_axis"]), 4),
    }
    b = part_b_and_c(atlas)

    # D: collinearity / identifiability on the 154-cell-type atlas
    zD = atlas["z_delta"].to_numpy()
    zK = atlas["z_kpna"].to_numpy()
    rho = spearman(zD, zK)
    from scipy.stats import pearsonr
    r = float(pearsonr(zD, zK).statistic)
    vif = 1.0 / max(1e-12, 1.0 - r ** 2)
    Zc = np.column_stack([zD, zK])
    Zc = Zc - Zc.mean(axis=0)
    sv = np.linalg.svd(Zc, compute_uv=False)
    pc1_ratio = float(sv[0] ** 2 / (sv ** 2).sum())
    pc1_load = sv[0] and None
    u, s, vt = np.linalg.svd(Zc, full_matrices=False)
    d = {
        "spearman_delta_vs_kpna": round(rho, 4),
        "pearson_r": round(r, 4),
        "r_squared": round(r ** 2, 4),
        "vif_each_axis": round(vif, 4),
        "pc1_variance_ratio": round(pc1_ratio, 4),
        "pc1_loadings_zD_zKPNA": [round(float(x), 4) for x in vt[0]],
        "n_cell_types": int(len(atlas)),
    }

    result = {
        "self_check_axis_rebuild": self_check,
        "partA_functional_panel_GSE46599": {},
        "partB_rank_drift_154": b,
        "partD_identifiability": d,
    }
    for rule in ("first", "median"):
        result["partA_functional_panel_GSE46599"][rule] = part_a(rule)

    with open(os.path.join(OUT, "composite_decision_main.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: result[k] for k in
                      ["self_check_axis_rebuild", "partD_identifiability"]},
                     ensure_ascii=False, indent=1))
    print(json.dumps(result["partB_rank_drift_154"], ensure_ascii=False, indent=1)[:4000])
    for rule in ("first", "median"):
        print("---- probe rule:", rule)
        print(json.dumps(result["partA_functional_panel_GSE46599"][rule]["forms"],
                         ensure_ascii=False, indent=1))
        print(json.dumps(result["partA_functional_panel_GSE46599"][rule]["fitted_two_factor"],
                         ensure_ascii=False))


if __name__ == "__main__":
    main()
