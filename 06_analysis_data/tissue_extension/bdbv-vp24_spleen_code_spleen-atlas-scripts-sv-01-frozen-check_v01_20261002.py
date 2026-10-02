"""SRB-VERIFY step 1: independent recompute of the frozen transform (criterion V1).

Task: analysis/_audit_tmp/_TASK_srb_verify.md  (item 2.1)
Independence rule: this module never imports anything from
analysis/spleen_reference_rebuild_20260926/ and re-derives every number from the
read-only HPA inputs.

Frozen transform (as re-derived here, not imported):
    x_g   = log10(nCPM_g + 1)
    z_g   = (x_g - mu_g) / sd_g,  mu/sd taken over the 154 HPA single-cell types
    comp1 = mean(z_IFNAR1, z_IFNAR2)        (type-I capacity)
    comp2 = mean(z_IFNLR1, z_IL10RB)        (type-III capacity)
    comp3 = mean(z_ISG15,  z_MX1)           (ISG priming)
    delta = mean(comp1, comp2, comp3)

Pre-frozen verdict rule
    V1a  all six score genes present for all 154 types            -> PASS/FAIL
    V1b  recomputed mu/sigma == restriction_gradient_summary.json -> |diff| < 1e-9
    V1c  recomputed delta     == restriction_gradient.tsv         -> |diff| < 1e-9

Also dumps the HPA whole-genome gene background (mean log10(nCPM+1) over the 154
types) that the frozen-scale null in sv_04 needs.

Writes:
    raw/sv_frozen_stats.json
    raw/sv_frozen_recheck.tsv
    raw/sv_hpa_gene_background.tsv
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUT = os.path.join(ROOT, "analysis", "srb_verify_20260926")
RAW = os.path.join(OUT, "raw")
os.makedirs(RAW, exist_ok=True)

T4 = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")
HPA_SC = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")

# component definition frozen by the main text (Pairs of genes, equally weighted)
COMPONENTS = {
    "IFN_I_capacity": ("IFNAR1", "IFNAR2"),
    "IFN_III_capacity": ("IFNLR1", "IL10RB"),
    "ISG_priming": ("ISG15", "MX1"),
}
SCORE_GENES = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]


def mean(xs):
    return sum(xs) / len(xs)


def sd(xs, mu):
    return math.sqrt(sum((v - mu) ** 2 for v in xs) / len(xs))


def main() -> None:
    wide_path = os.path.join(T4, "hpa_ifn_landscape_wide.tsv")
    with open(wide_path, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames
        rows = list(rd)
    cell_types = [r["cell_type"] for r in rows]
    print(f"V1: wide table {len(rows)} cell types, {len(cols) - 1} genes")

    missing_cells = {g: [r["cell_type"] for r in rows if r.get(g, "") == ""] for g in SCORE_GENES}
    v1a = all(not v for v in missing_cells.values())

    stats = {}
    for g in SCORE_GENES:
        xs = [math.log10(float(r[g]) + 1.0) for r in rows]
        mu, s = mean(xs), sd(xs, mean(xs))
        stats[g] = {"mean_log10_nCPM": mu, "sd_log10_nCPM": s, "n": len(xs)}

    with open(os.path.join(T4, "restriction_gradient_summary.json"), encoding="utf-8") as fh:
        ref = json.load(fh)["gene_stats"]
    stat_diffs = {g: {"d_mean": stats[g]["mean_log10_nCPM"] - ref[g]["mean_log10_nCPM"],
                      "d_sd": stats[g]["sd_log10_nCPM"] - ref[g]["sd_log10_nCPM"]}
                  for g in SCORE_GENES}
    v1b = max(max(abs(v["d_mean"]), abs(v["d_sd"])) for v in stat_diffs.values())

    def frozen_delta(vals):
        comps = {}
        for name, (a, b) in COMPONENTS.items():
            za = (math.log10(vals[a] + 1.0) - stats[a]["mean_log10_nCPM"]) / stats[a]["sd_log10_nCPM"]
            zb = (math.log10(vals[b] + 1.0) - stats[b]["mean_log10_nCPM"]) / stats[b]["sd_log10_nCPM"]
            comps[name] = (za + zb) / 2.0
        return mean(list(comps.values())), comps

    with open(os.path.join(T4, "restriction_gradient.tsv"), encoding="utf-8") as fh:
        ref_rows = {r["cell_type"]: r for r in csv.DictReader(fh, delimiter="\t")}

    out = []
    worst = 0.0
    for r in rows:
        vals = {g: float(r[g]) for g in SCORE_GENES}
        d, comps = frozen_delta(vals)
        rr = ref_rows.get(r["cell_type"])
        if rr is None:
            out.append({"cell_type": r["cell_type"], "delta": d, "ref_delta": None, "abs_diff": None,
                        "ref_rank": None})
            continue
        ref_delta = float(rr["delta_restriction"])
        ad = abs(d - ref_delta)
        worst = max(worst, ad)
        out.append({"cell_type": r["cell_type"], "delta": d, "ref_delta": ref_delta, "abs_diff": ad,
                    "ref_rank": rr["rank"],
                    "d_comp1": comps["IFN_I_capacity"] - float(rr["IFN_I_capacity"]),
                    "d_comp2": comps["IFN_III_capacity"] - float(rr["IFN_III_capacity"]),
                    "d_comp3": comps["ISG_priming"] - float(rr["ISG_priming"])})
    v1c = worst

    order = sorted(out, key=lambda x: -x["delta"])
    rank_mismatch = [r["cell_type"] for i, r in enumerate(order, 1)
                     if r["ref_rank"] is not None and int(r["ref_rank"]) != i]

    # --- whole-genome HPA background (for the frozen-scale null in sv_04) ---
    # Rule (see sv_01b): a gene enters the background iff it has exactly one row
    # per frozen reference type, i.e. 154 rows covering exactly the frozen panel.
    frozen_set = set(cell_types)
    per_gene = defaultdict(list)
    with open(HPA_SC, encoding="utf-8") as fh:
        rd2 = csv.DictReader(fh, delimiter="\t")
        hdr = rd2.fieldnames
        for r in rd2:
            try:
                per_gene[r["Gene name"]].append(
                    (r["Cell type"], math.log10(float(r["nCPM"]) + 1.0)))
            except (ValueError, KeyError, TypeError):
                continue
    widths = {len(v) for v in per_gene.values()}
    bg = {}
    for g, v in per_gene.items():
        if len(v) != len(cell_types) or {ct for ct, _ in v} != frozen_set:
            continue
        vals = [val for _, val in v]
        mu = mean(vals)
        bg[g] = {"mean_log10_nCPM": mu, "sd_log10_nCPM": sd(vals, mu)}
    n_types_full = len(cell_types)
    n_genes_covering = sum(1 for v in per_gene.values() if frozen_set <= {ct for ct, _ in v})
    with open(os.path.join(RAW, "sv_hpa_gene_background.tsv"), "w", encoding="utf-8", newline="") as fh:
        fh.write("gene\tmean_log10_nCPM\tsd_log10_nCPM\tn_observations\n")
        for g in sorted(bg, key=lambda k: bg[k]["mean_log10_nCPM"]):
            fh.write(f"{g}\t{bg[g]['mean_log10_nCPM']:.6f}\t{bg[g]['sd_log10_nCPM']:.6f}"
                     f"\t{n_types_full}\n")

    with open(os.path.join(RAW, "sv_frozen_recheck.tsv"), "w", encoding="utf-8", newline="") as fh:
        fh.write("cell_type\tdelta_recomputed\tdelta_reference\tabs_diff\tref_rank\n")
        for r in out:
            ref_txt = "" if r["ref_delta"] is None else f"{r['ref_delta']:.12f}"
            diff_txt = "" if r["abs_diff"] is None else f"{r['abs_diff']:.3e}"
            rank_txt = "" if r["ref_rank"] is None else str(r["ref_rank"])
            fh.write(f"{r['cell_type']}\t{r['delta']:.12f}\t{ref_txt}\t{diff_txt}\t{rank_txt}\n")

    summary = {
        "V1a": {"verdict": "PASS" if v1a else "FAIL",
                "n_cell_types": len(rows), "missing": {k: v for k, v in missing_cells.items() if v}},
        "V1b": {"verdict": "PASS" if v1b < 1e-9 else "FAIL", "max_abs_diff": v1b,
                "recomputed": stats, "reference": ref, "diffs": stat_diffs},
        "V1c": {"verdict": "PASS" if v1c < 1e-9 else "FAIL", "max_abs_diff": v1c,
                "n_compared": sum(1 for r in out if r["abs_diff"] is not None),
                "n_cell_types_wide": len(rows), "n_cell_types_reference": len(ref_rows),
                "rank_mismatches": rank_mismatch[:20], "n_rank_mismatches": len(rank_mismatch),
                "hpa_gene_background": {"n_types_matched": n_types_full, "n_genes": len(bg),
                                        "n_genes_covering_all_frozen_types": n_genes_covering,
                                        "rule": "exactly one row per frozen reference type",
                                        "width_values": sorted(widths)[:8]},
                "hpa_single_cell_header": hdr},
    }
    with open(os.path.join(RAW, "sv_frozen_stats.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    print(f"V1a all six genes x {len(rows)} types present: {'PASS' if v1a else 'FAIL'}")
    print(f"V1b max |mu/sigma diff| = {v1b:.3e} -> {'PASS' if v1b < 1e-9 else 'FAIL'}")
    print(f"V1c max |delta diff| over {summary['V1c']['n_compared']} types = {v1c:.3e} "
          f"-> {'PASS' if v1c < 1e-9 else 'FAIL'}")
    print(f"    rank mismatches: {len(rank_mismatch)}")
    print(f"    HPA background genes with complete {n_types_full}-type coverage: {len(bg)}")


if __name__ == "__main__":
    main()
