"""G6c: monocyte-only, per-cytokine donor-level test in GSE306664 (rank statistic,
exact permutation p over the 120 orderings), because monocytes are the only
population in that cohort with donor-level baseline contrast.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import os
import sys

BASE = os.path.join(r"G:\本迪布焦研究", "analysis", "g6c_primary_immune_20260926")
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


def exact_p(base: list[float], ind: list[float]) -> float:
    obs = abs(spearman(base, ind))
    hits = total = 0
    for perm in itertools.permutations(range(len(ind))):
        total += 1
        if abs(spearman(base, [ind[k] for k in perm])) >= obs - 1e-12:
            hits += 1
    return hits / total


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    expr: dict[str, dict[str, float]] = {}
    with open(os.path.join(BASE, "gse306664", "pseudobulk.tsv"), encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)[1:]
        for row in reader:
            expr[row[0]] = {header[i]: float(row[i + 1]) for i in range(len(header))}
    qc = list(csv.DictReader(open(os.path.join(BASE, "gse306664", "sample_qc.tsv"),
                                  encoding="utf-8"), delimiter="\t"))

    def score(gsm: str) -> float:
        return sum(expr[g][gsm] for g in ISG) / len(ISG)

    idx: dict[tuple[str, str], dict[str, str]] = {}
    for r in qc:
        idx.setdefault((r["cell_type"], r["treatment"]), {})[r["donor"]] = r["gsm"]

    out = {}
    for ct in ("Monocyte", "Bcell", "NK", "Tcell"):
        base = idx.get((ct, "none"), {})
        entry = {}
        for tr in ("IFNa", "IFNb", "IFNg", "IFN-L1"):
            trt = idx.get((ct, tr), {})
            donors = sorted(set(base) & set(trt))
            b = [score(base[d]) for d in donors]
            i = [score(trt[d]) - score(base[d]) for d in donors]
            rho = spearman(b, i)
            entry[tr] = {
                "n_donors": len(donors),
                "spearman_rank": round(rho, 3),
                "exact_two_sided_p": round(exact_p(b, i), 4) if len(donors) <= 7 else None,
            }
        out[ct] = entry
    path = os.path.join(BASE, "gse306664", "monocyte_by_cytokine.json")
    json.dump(out, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
