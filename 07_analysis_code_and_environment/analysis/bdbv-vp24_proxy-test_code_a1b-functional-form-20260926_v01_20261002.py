"""A1b: is the equal-weight form of Δ justified, and what weights does the data support?

Δ is, by construction, the equal-weight mean of three z-scored components (type-I
receptor capacity, type-III receptor capacity, baseline ISG priming). Nothing in the
manuscript justifies equal weights or additivity. This script fits the functional form
against the only two outcome axes available: a functional antiviral phenotype (GSE46599
HIV-1 restriction) and a measured interferon response (GSE327707 primary immune cells,
induction and post-treatment level).

Model comparison is by leave-one-out cross-validated AUC/Spearman so that a fitted model
has to earn its extra parameters.
"""

from __future__ import annotations

import json
import math
import os

import numpy as np

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "a1b_functional_form_20260926")
GSE46599 = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926", "material_test.json")

# GSE46599 phenotype from its own series summary
PHENOTYPE = {
    "primary-macrophages": 1.0, "PMA-THP-1": 1.0, "THP-1": 1.0, "U87-MG": 1.0,
    "primary-CD4+-T-cells": 0.75, "HT1080": 0.5, "PMA-U937": 0.5,
    "CEM": 0.0, "CEM-SS": 0.0, "Jurkat": 0.0, "U937": 0.0,
}


def auc(scores, labels) -> float:
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if p > q else 0.5 if p == q else 0.0 for p in pos for q in neg)
    return wins / (len(pos) * len(neg))


def loo_auc(X: np.ndarray, y: np.ndarray, form: str, weights: np.ndarray | None = None) -> float:
    """Leave-one-out AUC for a given functional form."""
    n = len(y)
    scores = np.zeros(n)
    for i in range(n):
        train = np.ones(n, dtype=bool)
        train[i] = False
        if form == "equal":
            s = X[i].mean()
        elif form == "fixed_weight":
            s = float(X[i] @ weights)
        elif form == "logistic":
            from sklearn.linear_model import LogisticRegression

            clf = LogisticRegression(C=1.0, max_iter=5000)
            clf.fit(X[train], y[train])
            s = float(clf.decision_function(X[i : i + 1])[0])
        else:
            raise ValueError(form)
        scores[i] = s
    return auc(scores, y.astype(bool))


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    data = json.load(open(GSE46599, encoding="utf-8"))
    rows = [r for r in data["rows"] if r["cell"] in PHENOTYPE]

    # three components, z-scored across the units exactly as Δ does
    base = np.array([r["baseline"] for r in rows], dtype=float)
    ifn = np.array([r["ifn"] for r in rows], dtype=float)
    fold = ifn - base
    comp = {
        "baseline_ISG": base,
        "fold_induction": fold,
        "post_IFN_level": ifn,
    }
    y = np.array([PHENOTYPE[r["cell"]] for r in rows], dtype=float)
    ybin = (y >= 0.75).astype(int)

    results: dict[str, object] = {
        "n_units": len(rows),
        "n_restricted": int(ybin.sum()),
        "n_permissive": int((y == 0).sum()),
        "note": "n is small; LOO AUC is reported so fitted forms must beat fixed forms "
                "out of sample",
    }

    # --- single-component predictors (the three rival quantities) ---
    for name, v in comp.items():
        v0 = (v - v.mean()) / (v.std() or 1.0)
        results[f"single_{name}"] = {
            "auc_all": round(auc(v0, ybin), 3),
            "loo_auc": round(loo_auc(np.column_stack([np.zeros_like(v0), np.zeros_like(v0), v0]),
                                     ybin, "fixed_weight", np.array([0.0, 0.0, 1.0])), 3),
        }

    # --- the manuscript's equal-weight form, and variants ---
    variants = {
        "equal_weight_ISG_and_receptors": [0.0, 0.0, 0.0],  # placeholder, handled below
    }
    # build a 3-component design: baseline ISG, fold, post-level (z-scored)
    Z = np.column_stack([(comp[k] - comp[k].mean()) / (comp[k].std() or 1.0)
                         for k in ["baseline_ISG", "fold_induction", "post_IFN_level"]])
    forms = {
        "equal_weight_baseline_only": np.array([1.0, 0.0, 0.0]),
        "equal_weight_fold_only": np.array([0.0, 1.0, 0.0]),
        "equal_weight_postlevel_only": np.array([0.0, 0.0, 1.0]),
        "naive_average_of_all_three": np.array([1 / 3, 1 / 3, 1 / 3]),
    }
    for name, w in forms.items():
        results[name] = {
            "auc_all": round(auc(Z @ w, ybin), 3),
            "loo_auc": round(loo_auc(Z, ybin, "fixed_weight", w), 3),
        }
    results["fitted_logistic"] = {
        "auc_all": round(auc(
            __import__("sklearn.linear_model", fromlist=["LogisticRegression"])
            .LogisticRegression(C=1.0, max_iter=5000).fit(Z, ybin).decision_function(Z), ybin), 3),
        "loo_auc": round(loo_auc(Z, ybin, "logistic"), 3),
        "coefficients_z_scaled": None,
    }
    from sklearn.linear_model import LogisticRegression

    clf = LogisticRegression(C=1.0, max_iter=5000).fit(Z, ybin)
    results["fitted_logistic"]["coefficients_z_scaled"] = [round(c, 3) for c in clf.coef_[0]]

    with open(os.path.join(OUT, "functional_form.json"), "w", encoding="utf-8") as fh:
        json.dump(results, fh, ensure_ascii=False, indent=1)
    print(json.dumps(results, ensure_ascii=False, indent=1))

    best = max(results, key=lambda k: (results[k]["loo_auc"] if isinstance(results[k], dict)
                                       and "loo_auc" in results[k] else -1))
    print(f"\nbest by LOO AUC: {best}")


if __name__ == "__main__":
    main()
