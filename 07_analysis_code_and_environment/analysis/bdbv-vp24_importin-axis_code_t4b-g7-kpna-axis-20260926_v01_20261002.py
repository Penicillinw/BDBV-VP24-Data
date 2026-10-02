"""G7: add the missing virus-side host axis to the restriction score.

The manuscript's Δ is built purely from interferon receptor capacity and baseline
ISG priming. But the BDBV VP24 defect acts by competing with PY-STAT1 for the
non-classical NLS site on karyopherin-α. A cell that is highly interferon
responsive but expresses little KPNA1/5/6 should therefore not pay the full
penalty, and Δ as currently defined cannot express that.

This script adds a second host axis — PY-STAT1 nuclear-import capacity
(KPNA1/KPNA5/KPNA6, the NPI-1 subfamily VP24 binds, plus the ISGF3 components
STAT1/STAT2/IRF9) — from the same local HPA single-cell-type table, then
characterises the joint (Δ, KPNA-capacity) space.

Outputs: analysis/t4b_g7_kpna_axis_20260926/
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")
WIDE = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925", "hpa_ifn_landscape_wide.tsv")
GRADIENT = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925", "restriction_gradient.tsv")
OUT = os.path.join(ROOT, "analysis", "t4b_g7_kpna_axis_20260926")

# genes that have to come from the source table because the frozen extraction
# did not request them
EXTRA_GENES = ["IRF9", "IFIT1", "JAK1", "JAK2"]


def read_wide() -> tuple[list[str], dict[str, dict[str, str]]]:
    rows = list(csv.DictReader(open(WIDE, encoding="utf-8"), delimiter="\t"))
    return [r["cell_type"] for r in rows], {r["cell_type"]: r for r in rows}


def read_extra(genes: set[str]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = defaultdict(dict)
    with open(SRC, encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        gi, ci, vi = idx["Gene name"], idx["Cell type"], idx["nCPM"]
        for row in reader:
            if row[gi] in genes:
                try:
                    out[row[gi]][row[ci]] = float(row[vi])
                except ValueError:
                    continue
    return out


def zmap(values: dict[str, float]) -> dict[str, float]:
    vals = list(values.values())
    mu = sum(vals) / len(vals)
    sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals)) or 1.0
    return {k: (v - mu) / sd for k, v in values.items()}


def spearman(a: list[float], b: list[float]) -> float:
    def rank(x: list[float]) -> list[float]:
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


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    cell_types, wide = read_wide()
    extra = read_extra(set(EXTRA_GENES))

    grad = {
        r["cell_type"]: float(r["delta_restriction"])
        for r in csv.DictReader(open(GRADIENT, encoding="utf-8"), delimiter="\t")
    }

    def series(gene: str) -> dict[str, float]:
        if gene in wide[cell_types[0]]:
            out = {}
            for ct in cell_types:
                raw = wide[ct].get(gene, "")
                try:
                    out[ct] = math.log10(float(raw) + 1.0)
                except (TypeError, ValueError):
                    out[ct] = 0.0
            return out
        raw = extra.get(gene, {})
        return {ct: math.log10(raw.get(ct, 0.0) + 1.0) for ct in cell_types}

    z = {g: zmap(series(g)) for g in ["KPNA1", "KPNA5", "KPNA6", "STAT1", "STAT2", "IRF9"]}

    def axis(genes: list[str]) -> dict[str, float]:
        return {ct: sum(z[g][ct] for g in genes) / len(genes) for ct in cell_types}

    kpna_axis = axis(["KPNA1", "KPNA5", "KPNA6"])
    isgf3_axis = axis(["STAT1", "STAT2", "IRF9"])
    import_axis = {ct: (kpna_axis[ct] + isgf3_axis[ct]) / 2 for ct in cell_types}

    rows = []
    for ct in cell_types:
        rows.append(
            {
                "cell_type": ct,
                "delta": grad[ct],
                "KPNA_axis": kpna_axis[ct],
                "ISGF3_axis": isgf3_axis[ct],
                "import_axis": import_axis[ct],
                "penalty_sum": grad[ct] + kpna_axis[ct],
                "KPNA1_z": z["KPNA1"][ct],
                "KPNA5_z": z["KPNA5"][ct],
                "KPNA6_z": z["KPNA6"][ct],
                "KPNA1_nCPM": float(wide[ct]["KPNA1"]),
                "KPNA5_nCPM": float(wide[ct]["KPNA5"]),
                "KPNA6_nCPM": float(wide[ct]["KPNA6"]),
            }
        )
    rows.sort(key=lambda r: -r["penalty_sum"])
    for i, r in enumerate(rows, 1):
        r["penalty_rank"] = i
    by_delta = sorted(rows, key=lambda r: -r["delta"])
    for i, r in enumerate(by_delta, 1):
        r["delta_rank"] = i
    by_kpna = sorted(rows, key=lambda r: -r["KPNA_axis"])
    for i, r in enumerate(by_kpna, 1):
        r["kpna_rank"] = i

    with open(os.path.join(OUT, "kpna_axis.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    high_delta = [r for r in rows if r["delta_rank"] <= 40]
    escape = sorted(high_delta, key=lambda r: r["kpna_rank"], reverse=True)[:12]
    top_penalty = rows[:15]
    rho_delta_kpna = spearman([r["delta"] for r in rows], [r["KPNA_axis"] for r in rows])
    rho_delta_import = spearman([r["delta"] for r in rows], [r["import_axis"] for r in rows])

    summary = {
        "n_cell_types": len(rows),
        "axes": {
            "delta": "frozen manuscript score (IFN-I capacity, IFN-III capacity, ISG priming)",
            "KPNA_axis": "mean z of KPNA1, KPNA5, KPNA6 (log10 nCPM, z across cell types)",
            "ISGF3_axis": "mean z of STAT1, STAT2, IRF9",
            "import_axis": "mean of KPNA_axis and ISGF3_axis",
        },
        "spearman_delta_vs_KPNA_axis": round(rho_delta_kpna, 3),
        "spearman_delta_vs_import_axis": round(rho_delta_import, 3),
        "top15_by_joint_penalty": [r["cell_type"] for r in top_penalty],
        "top12_high_delta_low_KPNA_escape_candidates": [
            {
                "cell_type": r["cell_type"],
                "delta_rank": r["delta_rank"],
                "kpna_rank": r["kpna_rank"],
                "KPNA1_nCPM": r["KPNA1_nCPM"],
                "KPNA5_nCPM": r["KPNA5_nCPM"],
                "KPNA6_nCPM": r["KPNA6_nCPM"],
            }
            for r in escape
        ],
        "kpna_axis_top10": [r["cell_type"] for r in by_kpna[:10]],
        "kpna_axis_bottom10": [r["cell_type"] for r in by_kpna[-10:]],
    }
    with open(os.path.join(OUT, "kpna_axis_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    print(f"cell types: {len(rows)}")
    print(f"Spearman(Delta, KPNA axis) = {rho_delta_kpna:.3f}")
    print(f"Spearman(Delta, import axis) = {rho_delta_import:.3f}")
    print("\nKPNA axis top 10:", ", ".join(r["cell_type"] for r in by_kpna[:10]))
    print("KPNA axis bottom 10:", ", ".join(r["cell_type"] for r in by_kpna[-10:]))
    print("\nhigh-Delta cells with the LOWEST KPNA capacity (escape candidates):")
    for r in escape:
        print(
            f"  {r['cell_type']:42s} delta_rank={r['delta_rank']:>3d} kpna_rank={r['kpna_rank']:>3d} "
            f"KPNA1={r['KPNA1_nCPM']:.1f} KPNA5={r['KPNA5_nCPM']:.1f} KPNA6={r['KPNA6_nCPM']:.1f}"
        )
    print("\njoint penalty top 15:", ", ".join(r["cell_type"] for r in top_penalty))


if __name__ == "__main__":
    main()
