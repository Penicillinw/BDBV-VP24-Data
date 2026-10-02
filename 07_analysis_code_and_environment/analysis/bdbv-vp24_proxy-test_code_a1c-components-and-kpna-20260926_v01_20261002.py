"""A1c: test the manuscript's OWN three components, and the sign of the KPNA axis,
against the one functional phenotype available (GSE46599 HIV-1 restriction).

The manuscript's Δ = mean of (type-I receptor capacity: IFNAR1+IFNAR2), (type-III
receptor capacity: IFNLR1+IL10RB), (baseline ISG priming: ISG15+MX1), each z-scored
across cell types. This script rebuilds exactly those components from the GSE46599
matrices, adds the importin-alpha axis (KPNA1+KPNA5+KPNA6), and evaluates every
candidate against the phenotype with leave-one-out AUC.
"""

from __future__ import annotations

import gzip
import json
import os
from collections import defaultdict

import numpy as np

ROOT = r"G:\本迪布焦研究"
RAW = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926", "raw")
OUT = os.path.join(ROOT, "analysis", "a1b_functional_form_20260926")

PHENOTYPE = {
    "primary-macrophages": 1.0, "PMA-THP-1": 1.0, "THP-1": 1.0, "U87-MG": 1.0,
    "primary-CD4+-T-cells": 0.75, "HT1080": 0.5, "PMA-U937": 0.5,
    "CEM": 0.0, "CEM-SS": 0.0, "Jurkat": 0.0, "U937": 0.0,
}
GROUPS = {
    "typeI": ["IFNAR1", "IFNAR2"],
    "typeIII": ["IFNLR1", "IL10RB"],
    "ISG_prime": ["ISG15", "MX1"],
    "KPNA": ["KPNA1", "KPNA5", "KPNA6"],
    "ISG_wide": ["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3", "BST2",
                 "RSAD2", "USP18", "STAT1", "IRF7", "IFI6"],
}


def load() -> tuple[list[str], list[str], dict[str, list[float]]]:
    header = None
    p2s: dict[str, str] = {}
    for line in open(os.path.join(RAW, "GPL10558.annot.gz"), "rb"):
        break
    with gzip.open(os.path.join(RAW, "GPL10558.annot.gz"), "rt",
                   encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                continue
            if header is None or line.startswith(("!", "#")) or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > 2 and p[2].strip() not in {"", "---"}:
                p2s[p[0].strip()] = p[2].strip().split("///")[0].strip()

    titles: list[str] = []
    rows: dict[str, list[float]] = {}
    with gzip.open(os.path.join(RAW, "GSE46599_series_matrix.txt.gz"), "rt",
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
                        rows[p[0].strip('"')] = [float(x) for x in p[1:]]
                    except ValueError:
                        pass
    gene_vec: dict[str, list[float]] = {}
    for probe, vec in rows.items():
        s = p2s.get(probe)
        if s and s not in gene_vec:
            gene_vec[s] = vec
    return titles, list(rows), gene_vec


def auc(scores, labels) -> float:
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float("nan")
    return sum(1.0 if p > q else 0.5 if p == q else 0.0 for p in pos for q in neg) / (
        len(pos) * len(neg))


def loo_auc(X: np.ndarray, y: np.ndarray) -> float:
    from sklearn.linear_model import LogisticRegression

    n = len(y)
    scores = np.zeros(n)
    for i in range(n):
        tr = np.ones(n, bool)
        tr[i] = False
        clf = LogisticRegression(C=1.0, max_iter=5000).fit(X[tr], y[tr])
        scores[i] = float(clf.decision_function(X[i : i + 1])[0])
    return auc(scores, y)


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    titles, _, genes = load()
    cells = []
    for t in titles:
        c = t.rsplit("_", 2)[0]
        if c not in cells:
            cells.append(c)

    def val(sym_list, cell, cond):
        idx = [i for i, t in enumerate(titles)
               if t.startswith(cell + "_") and f"_{cond}_" in t]
        vals = []
        for g in sym_list:
            v = genes.get(g)
            if v:
                vals.append(float(np.mean([v[i] for i in idx])))
        return float(np.mean(vals)) if vals else float("nan")

    rows = []
    for c in cells:
        if c not in PHENOTYPE:
            continue
        rec = {"cell": c, "phenotype": PHENOTYPE[c]}
        for name, syms in GROUPS.items():
            rec[name] = val(syms, c, "None")
        rec["ISG_induction"] = val(["ISG15", "MX1"], c, "IFN") - rec["ISG_prime"]
        rec["ISG_postlevel"] = val(["ISG15", "MX1"], c, "IFN")
        rows.append(rec)

    def z(key):
        v = np.array([r[key] for r in rows], float)
        return (v - v.mean()) / (v.std() or 1.0)

    y = np.array([r["phenotype"] >= 0.75 for r in rows])
    forms = {
        "manuscript_delta (typeI+typeIII+ISG)/3": (z("typeI") + z("typeIII") + z("ISG_prime")) / 3,
        "typeI_receptor_only": z("typeI"),
        "typeIII_receptor_only": z("typeIII"),
        "baseline_ISG_only": z("ISG_prime"),
        "ISG_wide_only": z("ISG_wide"),
        "KPNA_axis": z("KPNA"),
        "KPNA_axis_negated": -z("KPNA"),
        "ISG_induction": z("ISG_induction"),
        "ISG_postlevel": z("ISG_postlevel"),
    }
    out = {"n_units": len(rows), "n_restricted": int(y.sum()), "units": [r["cell"] for r in rows],
           "notes": "GSE46599; phenotype from the series summary (HIV-1 restriction)",
           "forms": {}}
    for name, s in forms.items():
        out["forms"][name] = {"auc": round(auc(s, y), 3)}

    # fitted model over the manuscript's three components
    C = np.column_stack([z("typeI"), z("typeIII"), z("ISG_prime")])
    from sklearn.linear_model import LogisticRegression

    clf = LogisticRegression(C=1.0, max_iter=5000).fit(C, y)
    out["fitted_three_components"] = {
        "auc_all": round(auc(clf.decision_function(C), y), 3),
        "loo_auc": round(loo_auc(C, y), 3),
        "z_scaled_coefficients": [round(c, 3) for c in clf.coef_[0]],
        "component_order": ["typeI", "typeIII", "ISG_prime"],
    }
    out["manuscript_form_loo_auc"] = round(loo_auc(C * 0 + (
        C.mean(axis=1, keepdims=True)), y), 3)  # equal weight, no fitting

    json.dump(out, open(os.path.join(OUT, "components_and_kpna.json"), "w",
                        encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
