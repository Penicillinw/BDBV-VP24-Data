"""G6c: population-level bridge to the manuscript's Delta.

The manuscript ranks CELL TYPES by a baseline composite (type-I receptor,
type-III receptor, ISG15+MX1). Here the same composite is computed from the
baseline samples of each sorted population in both primary-immune datasets and
compared with the population's measured IFN-alpha induction, i.e. the closest
available analogue of the manuscript's cross-cell-type use of Delta.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys

ROOT = r"G:\本迪布焦研究"
BASE = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926")
ISG15MX1 = ["ISG15", "MX1"]
TYPEI = ["IFNAR1", "IFNAR2"]
TYPEIII = ["IFNLR1", "IL10RB"]


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


def from_gse327707() -> list[dict]:
    rows = list(csv.DictReader(open(os.path.join(BASE, "gse327707", "sample_scores.tsv"),
                                    encoding="utf-8"), delimiter="\t"))
    by_ct: dict[str, dict[str, dict[str, float]]] = {}
    for r in rows:
        by_ct.setdefault(r["cell_type"], {})[r["treatment"]] = {
            k: float(r[k]) for k in ("isg", "isg_2gene", "typeI", "typeIII", "delta_like")
        }
    out = []
    for ct, tr in by_ct.items():
        if "Baseline" not in tr or "IFNa" not in tr:
            continue
        out.append({
            "population": ct,
            "mean_baseline_delta_like": round(tr["Baseline"]["delta_like"], 3),
            "mean_induction_isg": round(tr["IFNa"]["isg"] - tr["Baseline"]["isg"], 3),
        })
    return out


def from_gse306664() -> list[dict]:
    expr: dict[str, dict[str, float]] = {}
    with open(os.path.join(BASE, "gse306664", "pseudobulk.tsv"), encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)[1:]
        for row in reader:
            expr[row[0]] = {header[i]: float(row[i + 1]) for i in range(len(header))}
    qc = list(csv.DictReader(open(os.path.join(BASE, "gse306664", "sample_qc.tsv"),
                                  encoding="utf-8"), delimiter="\t"))
    by_ct: dict[str, dict[str, list[dict]]] = {}
    for r in qc:
        by_ct.setdefault(r["cell_type"], {}).setdefault(r["treatment"], []).append(r)

    def module(genes, gsm):
        present = [g for g in genes if g in expr]
        return sum(expr[g][gsm] for g in present) / len(present)

    out = []
    for ct, tr in by_ct.items():
        if "none" not in tr or "IFNa" not in tr:
            continue
        base = tr["none"]
        ifn = tr["IFNa"]
        # sample-matched delta-like and induction averages (same donors by construction)
        deltas, inds = [], []
        for b in base:
            gsm = b["gsm"]
            deltas.append((module(TYPEI, gsm) + module(TYPEIII, gsm) + module(ISG15MX1, gsm)) / 3)
            match = next((x for x in ifn if x["donor"] == b["donor"]), None)
            if match:
                inds.append(module(["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3",
                                    "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6"],
                                   match["gsm"])
                            - module(["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3",
                                      "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6"],
                                     gsm))
        out.append({
            "population": ct,
            "mean_baseline_delta_like": round(sum(deltas) / len(deltas), 3),
            "mean_induction_isg": round(sum(inds) / len(inds), 3),
        })
    return out


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    result = {}
    for name, rows in (("GSE327707", from_gse327707()), ("GSE306664", from_gse306664())):
        rho = spearman([r["mean_baseline_delta_like"] for r in rows],
                       [r["mean_induction_isg"] for r in rows]) if len(rows) >= 3 else None
        result[name] = {"rows": rows, "spearman_delta_like_vs_induction": round(rho, 3) if rho is not None else None}
    json.dump(result, open(os.path.join(BASE, "g6c_population_bridge.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(json.dumps(result, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
