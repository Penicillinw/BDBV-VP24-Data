"""A1: which quantity tracks the antiviral phenotype - baseline tone, the fold induction,
or the level reached after interferon?

Dataset: GSE46599, whose own series summary states the functional phenotype of every
member ("Some cell types become resistant to HIV-1 infection following type 1 IFN
treatment (such as macrophages, THP-1, PMA-THP-1, U87-MG cells and to a lesser extent,
primary CD4+ T cells) while others either become only partially resistant (e.g., HT1080,
PMA-U937) or remain permissive (e.g., CEM, CEM-SS, Jurkat T cell lines and U937)").

This is the level-versus-fold question the manuscript's Δ turns on: if the phenotype
follows the post-interferon LEVEL, the baseline score is a usable proxy and the gradient
keeps its original sign; if it follows the FOLD induction, the sign is reversed.
"""

from __future__ import annotations

import json
import math
import os

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "a1_level_vs_fold_20260926")
SRC = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926", "material_test.json")

# phenotype from the series summary; score 1 = restricted, 0.5 = partial, 0 = permissive
PHENOTYPE = {
    "primary-macrophages": 1.0,
    "PMA-THP-1": 1.0,
    "THP-1": 1.0,
    "U87-MG": 1.0,
    "primary-CD4+-T-cells": 0.75,   # "to a lesser extent"
    "HT1080": 0.5,
    "PMA-U937": 0.5,
    "CEM": 0.0,
    "CEM-SS": 0.0,
    "Jurkat": 0.0,
    "U937": 0.0,
}


def spearman(a: list[float], b: list[float]) -> float:
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


def auc(pos: list[float], neg: list[float]) -> float:
    """P(a random positive exceeds a random negative), ties counted as half."""
    if not pos or not neg:
        return float("nan")
    wins = 0.0
    for p in pos:
        for q in neg:
            wins += 1.0 if p > q else 0.5 if p == q else 0.0
    return wins / (len(pos) * len(neg))


def exact_permutation_p(values: list[float], labels: list[bool]) -> float:
    """Exact one-sided permutation p for the AUC of a two-group contrast.

    Enumerates every assignment of the observed values to the two group sizes, so the
    p-value carries no simulation error at these small n.
    """
    from itertools import combinations

    k = sum(labels)
    n = len(values)
    obs = auc([v for v, l in zip(values, labels) if l],
              [v for v, l in zip(values, labels) if not l])
    total = 0
    hits = 0
    for combo in combinations(range(n), k):
        sel = set(combo)
        pos = [values[i] for i in range(n) if i in sel]
        neg = [values[i] for i in range(n) if i not in sel]
        total += 1
        if auc(pos, neg) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (total + 1)


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    data = json.load(open(SRC, encoding="utf-8"))
    rows = [r for r in data["rows"] if r["cell"] in PHENOTYPE]
    for r in rows:
        r["phenotype"] = PHENOTYPE[r["cell"]]
    rows.sort(key=lambda r: -r["phenotype"])

    print(f"{'cell':22s} {'kind':13s} {'pheno':>5s} {'baseline':>9s} {'induced':>8s} "
          f"{'fold':>7s}")
    for r in rows:
        print(f"{r['cell']:22s} {r['kind']:13s} {r['phenotype']:5.2f} {r['baseline']:9.3f} "
              f"{r['ifn']:8.3f} {r['induction']:7.3f}")

    pheno = [r["phenotype"] for r in rows]
    tests = {
        "baseline_level": [r["baseline"] for r in rows],
        "fold_induction": [r["induction"] for r in rows],
        "post_IFN_level": [r["ifn"] for r in rows],
    }
    results = {}
    restricted = [r for r in rows if r["phenotype"] >= 0.75]
    permissive = [r for r in rows if r["phenotype"] == 0.0]
    for name, vals in tests.items():
        labels = [r["phenotype"] >= 0.75 or r["phenotype"] == 0.0 for r in rows]
        keep = [i for i, r in enumerate(rows) if r["phenotype"] >= 0.75 or r["phenotype"] == 0.0]
        results[name] = {
            "spearman_vs_phenotype": round(spearman(vals, pheno), 3),
            "auc_restricted_vs_permissive": round(
                auc([v for v, r in zip(vals, rows) if r["phenotype"] >= 0.75],
                    [v for v, r in zip(vals, rows) if r["phenotype"] == 0.0]), 3),
            "exact_permutation_p": round(
                exact_permutation_p([vals[i] for i in keep],
                                    [rows[i]["phenotype"] >= 0.75 for i in keep]), 4),
            "mean_restricted": round(
                sum(v for v, r in zip(vals, rows) if r["phenotype"] >= 0.75)
                / max(1, len(restricted)), 3),
            "mean_permissive": round(
                sum(v for v, r in zip(vals, rows) if r["phenotype"] == 0.0)
                / max(1, len(permissive)), 3),
        }
    results["n_units"] = len(rows)
    results["n_immortalized_only"] = sum(1 for r in rows if r["kind"] == "immortalized")

    # same three tests restricted to immortalised lines only
    imm = [r for r in rows if r["kind"] == "immortalized"]
    for name in ["baseline", "induction", "ifn"]:
        results[f"immortalized_only_{name}"] = {
            "spearman_vs_phenotype": round(
                spearman([r[name] for r in imm], [r["phenotype"] for r in imm]), 3),
            "auc_restricted_vs_permissive": round(
                auc([r[name] for r in imm if r["phenotype"] >= 0.75],
                    [r[name] for r in imm if r["phenotype"] == 0.0]), 3),
        }

    with open(os.path.join(OUT, "level_vs_fold.json"), "w", encoding="utf-8") as fh:
        json.dump({"rows": rows, "results": results,
                   "phenotype_source": "GSE46599 !Series_summary (Goujon & Malim)"},
                  fh, ensure_ascii=False, indent=1)
    print()
    print(json.dumps(results, ensure_ascii=False, indent=1))

    best = max(["baseline_level", "fold_induction", "post_IFN_level"],
               key=lambda k: abs(results[k]["auc_restricted_vs_permissive"] - 0.5))
    print(f"\nbest discriminator by |AUC-0.5|: {best}")


if __name__ == "__main__":
    main()
