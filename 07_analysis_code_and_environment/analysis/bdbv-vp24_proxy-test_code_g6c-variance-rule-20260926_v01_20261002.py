"""G6c: symmetric 'informative population' sensitivity across both primary-immune
datasets.

MOTIVATION (declared after the first pass, therefore labelled EXPLORATORY).
GSE306664's donor-level baseline ISG contrast is essentially zero in B, NK and T
cells (between-donor range 0.07-0.08 log units) and large in monocytes (range
2.41). A donor-level correlation computed inside a population with no baseline
contrast is noise, and pooling three such populations dilutes the one population
that can actually be tested. GSE327707 has real contrast in all five populations.

RULE (scale-free, applied identically to both datasets):
  a population contributes to the donor-level test only if its between-donor SD of
  baseline ISG score is at least 20% of the between-population SD of baseline
  means in that dataset.

This script also reports the monocyte-only test in both datasets, because CD14+
monocytes are the manuscript's key myeloid compartment.
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
import sys

ROOT = r"G:\本迪布焦研究"
BASE = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926")

ISG = ["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3",
       "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6"]


def spearman(a, b) -> float:
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


def rank_within(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def perm_p(x, y, rho, draws=20000, seed=20260926):
    rng = random.Random(seed)
    yy = list(y)
    hits = 0
    for _ in range(draws):
        rng.shuffle(yy)
        if abs(spearman(x, yy)) >= abs(rho) - 1e-12:
            hits += 1
    return (hits + 1) / (draws + 1)


def load_dataset(name: str, path: str) -> dict[str, dict[tuple[str, str], float]]:
    """(cell_type, donor) -> {'baseline': v, 'induction': v} for the control/treatment pair."""
    expr: dict[str, dict[str, float]] = {}
    with open(path, encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)[1:]
        for row in reader:
            expr[row[0]] = {header[i]: float(row[i + 1]) for i in range(len(header))}
    return expr


def report(dataset: str, per_unit: dict[str, dict[str, tuple[float, float]]]) -> dict:
    cell_types = sorted(per_unit)
    means = {ct: sum(v[0] for v in per_unit[ct].values()) / len(per_unit[ct]) for ct in cell_types}
    grand = sum(means.values()) / len(means)
    across_sd = math.sqrt(sum((m - grand) ** 2 for m in means.values()) / (len(means) - 1))
    out = {"across_population_sd_of_baseline_means": round(across_sd, 3), "populations": {}}
    informative: list[str] = []
    for ct in cell_types:
        base = [v[0] for v in per_unit[ct].values()]
        ind = [v[1] for v in per_unit[ct].values()]
        sd = math.sqrt(sum((b - sum(base) / len(base)) ** 2 for b in base) / (len(base) - 1))
        ratio = sd / across_sd if across_sd else float("nan")
        entry = {
            "n_donors": len(base),
            "baseline_values": [round(b, 3) for b in base],
            "induction_values": [round(i, 3) for i in ind],
            "between_donor_baseline_sd": round(sd, 3),
            "sd_ratio_to_across_population_sd": round(ratio, 3),
            "informative_by_rule": bool(ratio >= 0.2),
            "spearman_within_population": round(spearman(base, ind), 3),
        }
        if entry["informative_by_rule"]:
            informative.append(ct)
        out["populations"][ct] = entry
    xs: list[float] = []
    ys: list[float] = []
    for ct in informative:
        xs += rank_within([v[0] for v in per_unit[ct].values()])
        ys += rank_within([v[1] for v in per_unit[ct].values()])
    if len(xs) >= 6:
        rho = spearman(xs, ys)
        out["pooled_informative_populations"] = {
            "populations": informative,
            "n_units": len(xs),
            "spearman": round(rho, 3),
            "p_permutation": round(perm_p(xs, ys, rho), 4),
        }
    else:
        out["pooled_informative_populations"] = {"populations": informative, "n_units": len(xs)}
    return out


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass

    # ---- GSE327707: CD14 pair from the sample score table -------------------
    rows = list(csv.DictReader(open(os.path.join(BASE, "gse327707", "sample_scores.tsv"),
                                    encoding="utf-8"), delimiter="\t"))
    by_donor: dict[str, dict[str, dict[str, float]]] = {}
    for r in rows:
        by_donor.setdefault(r["cell_type"], {}).setdefault(r["donor"], {})[r["treatment"]] = float(r["isg"])
    per_unit_1: dict[str, dict[str, tuple[float, float]]] = {}
    for ct, donors in by_donor.items():
        per_unit_1[ct] = {}
        for donor, tr in donors.items():
            if "Baseline" in tr and "IFNa" in tr:
                per_unit_1[ct][donor] = (tr["Baseline"], tr["IFNa"] - tr["Baseline"])

    # ---- GSE306664 ----------------------------------------------------------
    qc = list(csv.DictReader(open(os.path.join(BASE, "gse306664", "sample_qc.tsv"),
                                  encoding="utf-8"), delimiter="\t"))
    expr = load_dataset("GSE306664", os.path.join(BASE, "gse306664", "pseudobulk.tsv"))
    present = [g for g in ISG if g in expr]

    def score(gsm: str) -> float:
        return sum(expr[g][gsm] for g in present) / len(present)

    per_unit_2: dict[str, dict[str, tuple[float, float]]] = {}
    idx = {(r["cell_type"], r["donor"], r["treatment"]): r["gsm"] for r in qc}
    for (ct, donor, tr), gsm in idx.items():
        if tr not in {"none", "IFNa"}:
            continue
        per_unit_2.setdefault(ct, {}).setdefault(donor, {})[tr] = score(gsm)
    per_unit_2 = {
        ct: {d: (tr["none"], tr["IFNa"] - tr["none"])
             for d, tr in donors.items() if "none" in tr and "IFNa" in tr}
        for ct, donors in per_unit_2.items()
    }
    per_unit_2 = {ct: v for ct, v in per_unit_2.items() if v}

    result = {
        "rule": "population contributes only if between-donor baseline SD >= 20% of the "
                "between-population SD of baseline means (exploratory, declared 2026-09-26)",
        "GSE327707": report("GSE327707", per_unit_1),
        "GSE306664": report("GSE306664", per_unit_2),
    }
    for name in ("GSE327707", "GSE306664"):
        mono = result[name]["populations"].get(
            "CD14" if name == "GSE327707" else "Monocyte")
        result[name]["monocyte_only"] = mono

    json.dump(result, open(os.path.join(BASE, "g6c_variance_rule.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(json.dumps(result, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
