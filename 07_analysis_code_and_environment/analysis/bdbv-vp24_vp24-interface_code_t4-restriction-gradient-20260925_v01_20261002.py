"""Compute the pre-registered restriction-gradient input score Δ(c) per cell type.

Formula (frozen, see report/T4_预注册_...):
    z(IFN-I capacity)   = mean z(IFNAR1, IFNAR2)
    z(IFN-III capacity) = mean z(IFNLR1, IL10RB)
    z(ISG priming)      = mean z(ISG15, MX1)      # IFITM3 excluded: strong lineage
                                                  # skew (decidual/epididymal epithelium)
    Delta(c) ~ equal-weight mean of the three

All z-scores are computed across the 154 HPA cell types (value + 1, log10
transformed first, to tame the ~4-order-of-magnitude nCPM range).
"""

import csv
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WIDE = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925", "hpa_ifn_landscape_wide.tsv")
OUT = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")

GROUPS = {
    "IFN_I_capacity": ["IFNAR1", "IFNAR2"],
    "IFN_III_capacity": ["IFNLR1", "IL10RB"],
    "ISG_priming": ["ISG15", "MX1"],
}


def zscores(values):
    vals = [v for v in values.values() if v is not None]
    mu = sum(vals) / len(vals)
    sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals)) or 1.0
    return {k: ((v - mu) / sd if v is not None else None) for k, v in values.items()}, mu, sd


def main():
    rows = list(csv.DictReader(open(WIDE, encoding="utf-8"), delimiter="\t"))
    cell_types = [r["cell_type"] for r in rows]

    def series(gene):
        out = {}
        for r in rows:
            raw = r.get(gene, "")
            try:
                out[r["cell_type"]] = math.log10(float(raw) + 1.0)
            except ValueError:
                out[r["cell_type"]] = None
        return out

    z = {}
    meta = {}
    for group, genes in GROUPS.items():
        per_gene = {}
        for g in genes:
            sg, mu, sd = zscores(series(g))
            per_gene[g] = sg
            meta[g] = {"mean_log10_nCPM": mu, "sd_log10_nCPM": sd}
        z[group] = {
            ct: sum(per_gene[g][ct] for g in genes if per_gene[g][ct] is not None)
            / max(1, sum(1 for g in genes if per_gene[g][ct] is not None))
            for ct in cell_types
        }

    scored = []
    for ct in cell_types:
        parts = [z[g][ct] for g in GROUPS if z[g][ct] is not None]
        delta = sum(parts) / len(parts) if parts else None
        scored.append(
            {
                "cell_type": ct,
                "delta_restriction": delta,
                **{g: z[g][ct] for g in GROUPS},
                "IFNLR1_nCPM": float(next(r["IFNLR1"] for r in rows if r["cell_type"] == ct) or 0),
            }
        )
    scored.sort(key=lambda d: -(d["delta_restriction"] or -99))
    for rank, rec in enumerate(scored, 1):
        rec["rank"] = rank

    with open(os.path.join(OUT, "restriction_gradient.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(scored[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(scored)

    quantiles = {}
    for q in (0.1, 0.25, 0.5, 0.75, 0.9):
        idx = int(q * (len(scored) - 1))
        quantiles[str(q)] = scored[idx]["delta_restriction"]

    json.dump(
        {
            "formula": {g: GROUPS[g] for g in GROUPS},
            "transform": "log10(nCPM + 1), z-scored across 154 HPA cell types",
            "gene_stats": meta,
            "n_cell_types": len(scored),
            "delta_quantiles": quantiles,
            "top20": scored[:20],
            "bottom20": scored[-20:],
        },
        open(os.path.join(OUT, "restriction_gradient_summary.json"), "w", encoding="utf-8"),
        ensure_ascii=False,
        indent=1,
    )

    print(f"{'rank':>4s} {'Delta':>7s} {'IFN-I':>7s} {'IFN-III':>8s} {'ISG':>7s}  cell_type")
    for rec in scored[:25]:
        print(f"{rec['rank']:>4d} {rec['delta_restriction']:>7.2f} {rec['IFN_I_capacity']:>7.2f} "
              f"{rec['IFN_III_capacity']:>8.2f} {rec['ISG_priming']:>7.2f}  {rec['cell_type']}")
    print("...")
    for rec in scored[-12:]:
        print(f"{rec['rank']:>4d} {rec['delta_restriction']:>7.2f} {rec['IFN_I_capacity']:>7.2f} "
              f"{rec['IFN_III_capacity']:>8.2f} {rec['ISG_priming']:>7.2f}  {rec['cell_type']}")


if __name__ == "__main__":
    main()
