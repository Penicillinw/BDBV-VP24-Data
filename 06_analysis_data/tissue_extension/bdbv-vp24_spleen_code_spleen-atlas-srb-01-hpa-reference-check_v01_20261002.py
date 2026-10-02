"""SRB step 1 - freeze and re-verify the reference transform DELTA (criteria R1/R2).

The frozen score in the manuscript is

    comp_M  = mean over genes g in M of z_g,   z_g = (log10(nCPM_g + 1) - mu_g) / sd_g
    DELTA   = (comp_IFN_I + comp_IFN_III + comp_ISG_priming) / 3

with mu_g, sd_g computed **across the 154 HPA reference cell types** (population
sd, ddof = 0), and the modules

    IFN_I_capacity   = IFNAR1, IFNAR2
    IFN_III_capacity = IFNLR1, IL10RB
    ISG_priming      = ISG15,  MX1

This script re-derives mu/sd from the raw HPA wide table and checks that it
reproduces (a) the gene_stats recorded in restriction_gradient_summary.json and
(b) every one of the 154 frozen delta_restriction values.  It exists so that the
spleen extension can reuse a *verified* transform rather than an assumed one.

Criteria (frozen before the numbers were read):
  R1a  all six genes present in the wide table for all 154 cell types;
  R1b  |mu_recomputed - mu_recorded| < 1e-9 and |sd_recomputed - sd_recorded| < 1e-9;
  R1c  max |DELTA_recomputed - delta_restriction| < 1e-9 and ranks identical.

Writes raw/hpa_gene_stats_verified.json and raw/hpa_delta_recheck.tsv.
"""

from __future__ import annotations

import csv
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926")
RAW = os.path.join(OUT, "raw")
T4 = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")

MODULE = {
    "IFN_I_capacity": ["IFNAR1", "IFNAR2"],
    "IFN_III_capacity": ["IFNLR1", "IL10RB"],
    "ISG_priming": ["ISG15", "MX1"],
}
SIX = [g for v in MODULE.values() for g in v]


def spearman(a, b):
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


def main() -> None:
    wide = list(csv.DictReader(open(os.path.join(T4, "hpa_ifn_landscape_wide.tsv"), encoding="utf-8"), delimiter="\t"))
    frozen = {r["cell_type"]: r for r in
              csv.DictReader(open(os.path.join(T4, "restriction_gradient.tsv"), encoding="utf-8"), delimiter="\t")}
    summary = json.load(open(os.path.join(T4, "restriction_gradient_summary.json"), encoding="utf-8"))

    # R1a
    missing = [(r["cell_type"], g) for r in wide for g in SIX if g not in r or r[g] in ("", None)]
    assert not missing, f"R1a FAIL: {missing[:5]}"
    print(f"R1a PASS: {len(wide)} cell types x {len(SIX)} genes complete")

    logc = {g: {r["cell_type"]: math.log10(float(r[g]) + 1.0) for r in wide} for g in SIX}
    types = [r["cell_type"] for r in wide]
    assert len(types) == 154 and len(set(types)) == 154, f"expected 154 unique cell types, got {len(types)}"

    # R1b
    stats = {}
    worst = 0.0
    for g in SIX:
        xs = list(logc[g].values())
        mu = sum(xs) / len(xs)
        sd = math.sqrt(sum((v - mu) ** 2 for v in xs) / len(xs))
        rec = summary["gene_stats"][g]
        d_mu = abs(mu - rec["mean_log10_nCPM"])
        d_sd = abs(sd - rec["sd_log10_nCPM"])
        worst = max(worst, d_mu, d_sd)
        stats[g] = {"mean_log10_nCPM": mu, "sd_log10_nCPM": sd,
                    "recorded_mean": rec["mean_log10_nCPM"], "recorded_sd": rec["sd_log10_nCPM"]}
    assert worst < 1e-9, f"R1b FAIL: worst |delta| = {worst:g}"
    print(f"R1b PASS: recomputed mu/sd match the recorded gene_stats (worst |delta| = {worst:.2e})")

    # R1c
    def delta_of(logmap):
        comps = []
        for genes in MODULE.values():
            comps.append(sum((logmap[g] - stats[g]["mean_log10_nCPM"]) / stats[g]["sd_log10_nCPM"]
                             for g in genes) / len(genes))
        return comps

    rows, worst_d, mismatched = [], 0.0, 0
    for t in types:
        comps = delta_of({g: logc[g][t] for g in SIX})
        d = sum(comps) / 3.0
        ref = float(frozen[t]["delta_restriction"])
        worst_d = max(worst_d, abs(d - ref))
        if abs(d - ref) > 1e-9:
            mismatched += 1
        rows.append({"cell_type": t, "delta_recomputed": d, "delta_frozen": ref, "abs_diff": abs(d - ref),
                     "IFN_I": comps[0], "IFN_III": comps[1], "ISG_priming": comps[2]})
    assert worst_d < 1e-9 and mismatched == 0, f"R1c FAIL: {mismatched} mismatches, worst {worst_d:g}"
    rho = spearman([r["delta_recomputed"] for r in rows], [r["delta_frozen"] for r in rows])
    print(f"R1c PASS: 154/154 delta identical (worst |delta| = {worst_d:.2e}, Spearman = {rho:.6f})")

    os.makedirs(RAW, exist_ok=True)
    with open(os.path.join(RAW, "hpa_delta_recheck.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(RAW, "hpa_gene_stats_verified.json"), "w", encoding="utf-8") as fh:
        json.dump({"modules": MODULE, "transform": "log10(nCPM+1), z-scored across the 154 HPA cell types",
                   "gene_stats": stats, "reference_panel": {"n_cell_types": len(types),
                                                            "source": "analysis/t4_ifn_landscape_20260925/hpa_ifn_landscape_wide.tsv"},
                   "checks": {"R1a": "PASS", "R1b": f"PASS (worst {worst:.2e})",
                              "R1c": f"PASS ({mismatched} mismatches, worst {worst_d:.2e}, spearman {rho:.6f})"}},
                  fh, indent=1)
    print("wrote raw/hpa_gene_stats_verified.json and raw/hpa_delta_recheck.tsv")


if __name__ == "__main__":
    main()
