"""A1 calibration (v2): is the "level vs fold" question *decidable* in GSE46599,
and does the manuscript's module-sensitivity margin survive a fixed group definition?

WHY A v2 EXISTS (revision record, kept explicit)
------------------------------------------------
v1 of this script pre-registered five rules (R1-R5, below) and computed the
baseline -> induced-level continuum for six modules using
analysis/g6e_gse46599_20260926/g6e_results.json (a single, self-consistent probe
resolution). It returned three facts that contradict the manuscript's A1 sentence:

  (i)  R3 failed: in 3 of 6 modules the induced level did NOT beat both rivals;
  (ii) the panel is degenerate for most modules (rank correlation between baseline
       and induced level 0.91-0.98; sd(fold)/sd(baseline) 0.20-0.50), so "level" and
       "fold" are not independent contrasts here;
  (iii) the manuscript's sentence "0.90 against 0.70 for a two-gene module to 0.83
       against 0.67 for a seventeen-gene module" could not be reproduced from any one
       artifact. Tracing it showed the two ends come from two different scripts that
       use different definitions of the NEGATIVE group:
         - analysis/a1_level_vs_fold_20260926/level_vs_fold.json
             negatives = the 4 permissive units only  -> induced 0.90, baseline 0.75
         - analysis/a1b_functional_form_20260926/functional_form.json
             negatives = 4 permissive + 2 "partially resistant" units -> induced 0.833,
             rivals 0.667 / 0.767
       So the quoted margin mixes a change of MODULE with a change of GROUP DEFINITION.

v2 therefore adds the missing axis and nothing else: the same six modules are now
evaluated under three fixed group rules. The pre-registered rules R1-R5 are unchanged;
the group rule becomes an explicit factor instead of an implicit one. This is a
pre-registration *revision*, recorded here and in FINDINGS.md, not a silent change.

DATA (all on disk; no downloads)
--------------------------------
analysis/g6e_gse46599_20260926/g6e_results.json
  per_cell_type: 11 units x {module}_baseline / {module}_IFN (log2 mean expression,
  one probe resolution, one annotation pass)
  hiv1_restriction_label_from_series: {cell: [group, ...]} from the series text

------------------------------------------------------------------------------
PRE-REGISTERED DECISION RULES (fixed before the numbers were read; revision above)

R1 DEGENERACY. Baseline tone and induced level count as one axis in this panel if,
   over the 11 units, Spearman(baseline, induced) >= 0.80 AND the induced-on-baseline
   OLS slope >= 0.80 with R^2 >= 0.60; equivalently sd(fold)/sd(baseline) < 0.50.

R2 FLAT PROFILE. For a given module x group rule, the continuum is uninformative if
   max_theta AUC - min_theta AUC <= 0.10.

R3 RANKING STABILITY. "Post-interferon level separates the groups better than baseline
   tone and than fold induction" counts as supported ONLY if
   AUC(induced) >= max(AUC(baseline), AUC(fold)) holds for EVERY module AND EVERY
   group rule tested. Every violating cell is reported by name.

R4 FAMILY-WISE SIGNIFICANCE. A contrast is called significant only if its exact
   permutation p on the max-statistic over the whole theta grid (i.e. corrected for
   having searched the continuum) is < 0.05. Pointwise p is reported alongside.

R5 POWER. The exact permutation p can only take the values k/N for N = C(n, 5) label
   assignments, so N, the smallest achievable p and the largest head-count that still
   gives p < 0.05 are reported for every group rule.
------------------------------------------------------------------------------

Group rules (fixed):
  A partials_excluded     : positive = 5 resistant, negative = 4 permissive  (n = 9)
  B partials_as_permissive: positive = 5 resistant, negative = 4 + 2 partial (n = 11)
  C ordinal_all_units     : Spearman of the score against the ordinal grade
                            (1.0 / 0.75 / 0.5 / 0.0) over all 11 units

Outputs (this directory only): contrast_grid_v2.tsv, contrast_profile_v2.json.
"""

from __future__ import annotations

import itertools
import json
import math
import os
from statistics import mean, stdev

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926", "g6e_results.json")
OUT = os.path.join(ROOT, "analysis", "a1_level_calibration_20260926")

MODULES = {
    "ISG_17gene": ("ISG_baseline", "ISG_IFN"),
    "priming_ISG15_MX1": ("priming_ISG15_MX1_baseline", "priming_ISG15_MX1_IFN"),
    "typeI_receptor": ("typeI_IFNAR1_IFNAR2_baseline", "typeI_IFNAR1_IFNAR2_IFN"),
    "typeIII_receptor": ("typeIII_IFNLR1_IL10RB_baseline", "typeIII_IFNLR1_IL10RB_IFN"),
    "receptor_4gene": ("receptor_baseline", "receptor_IFN"),
    "housekeeping_negcontrol": ("housekeeping_baseline", "housekeeping_IFN"),
}
GRADE = {"resistant": 1.0, "partially resistant": 0.5, "permissive": 0.0}
THETAS = list(range(0, 91, 1))  # degrees; 0 = baseline only, 90 = induced only


def _rank(vals: list[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    r = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(xs: list[float], ys: list[float]) -> float:
    rx, ry = _rank(xs), _rank(ys)
    mx, my = mean(rx), mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def ols(y: list[float], x: list[float]) -> tuple[float, float]:
    mx, my = mean(x), mean(y)
    sxx = sum((v - mx) ** 2 for v in x)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    slope = sxy / sxx
    intercept = my - slope * mx
    ss_tot = sum((b - my) ** 2 for b in y)
    ss_res = sum((b - (intercept + slope * a)) ** 2 for a, b in zip(x, y))
    return slope, (1.0 - ss_res / ss_tot if ss_tot else float("nan"))


def zscore(vals: list[float]) -> list[float]:
    s = stdev(vals)
    m = mean(vals)
    return [(v - m) / s for v in vals] if s else [0.0] * len(vals)


def auc(pos: list[float], neg: list[float]) -> float:
    gt = eq = 0
    for a in pos:
        for b in neg:
            if a > b:
                gt += 1
            elif a == b:
                eq += 1
    return (gt + 0.5 * eq) / (len(pos) * len(neg))


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    data = json.load(open(SRC, encoding="utf-8"))
    units = data["per_cell_type"]
    labels = {k: v[0] for k, v in data["hiv1_restriction_label_from_series"].items()}
    cell_index = {u["cell"]: k for k, u in enumerate(units)}

    groups = {
        g: [u["cell"] for u in units if labels[u["cell"]] == g]
        for g in ("resistant", "partially resistant", "permissive")
    }

    rules = {
        "A_partials_excluded": {
            "pos": [cell_index[c] for c in groups["resistant"]],
            "neg": [cell_index[c] for c in groups["permissive"]],
        },
        "B_partials_as_permissive": {
            "pos": [cell_index[c] for c in groups["resistant"]],
            "neg": [cell_index[c] for c in groups["permissive"] + groups["partially resistant"]],
        },
    }

    summary: dict = {
        "source": os.path.relpath(SRC, ROOT).replace("\\", "/"),
        "units": [u["cell"] for u in units],
        "groups": groups,
        "rules": {
            "R1_degeneracy": "spearman>=0.80 and slope>=0.80 and R2>=0.60; sd(fold)/sd(base)<0.50",
            "R2_flat_profile": "max AUC - min AUC <= 0.10 over theta grid",
            "R3_ranking_stability": "AUC(induced) >= max(AUC(baseline), AUC(fold)) in every module x rule",
            "R4_family_wise": "max-statistic exact permutation p over the theta grid < 0.05",
            "R5_power": "exact p grid = k/N for N = C(n, 5)",
        },
        "permutation_space": {},
        "degeneracy": {},
        "grid": {},
        "ordinal_all_units": {},
    }
    for name, r in rules.items():
        n = len(r["pos"]) + len(r["neg"])
        N = math.comb(n, len(r["pos"]))
        summary["permutation_space"][name] = {
            "n_pos": len(r["pos"]),
            "n_neg": len(r["neg"]),
            "n_assignments": N,
            "min_achievable_p": round(1.0 / N, 5),
            "max_headcount_for_p_below_0.05": max(k for k in range(1, N + 1) if k / N < 0.05),
        }

    rows_out: list[dict] = []

    for mod, (kb, ki) in MODULES.items():
        base_all = [u[kb] for u in units]
        ind_all = [u[ki] for u in units]
        fold_all = [i - b for b, i in zip(base_all, ind_all)]

        rho = spearman(base_all, ind_all)
        slope, r2 = ols(ind_all, base_all)
        sd_ratio = stdev(fold_all) / stdev(base_all)
        summary["degeneracy"][mod] = {
            "spearman_baseline_vs_induced_n11": round(rho, 4),
            "ols_slope_induced_on_baseline_n11": round(slope, 4),
            "ols_r2_n11": round(r2, 4),
            "sd_fold_over_sd_baseline_n11": round(sd_ratio, 4),
            "fold_range_log2_n11": [round(min(fold_all), 3), round(max(fold_all), 3)],
            "R1_degenerate": bool(rho >= 0.80 and slope >= 0.80 and r2 >= 0.60),
        }

        zb, zi = zscore(base_all), zscore(ind_all)
        grade_vec = [GRADE[labels[u["cell"]]] for u in units]
        summary["ordinal_all_units"][mod] = {
            "spearman_baseline_vs_grade": round(spearman(base_all, grade_vec), 4),
            "spearman_fold_vs_grade": round(spearman(fold_all, grade_vec), 4),
            "spearman_induced_vs_grade": round(spearman(ind_all, grade_vec), 4),
        }

        summary["grid"][mod] = {}
        for rule_name, r in rules.items():
            pos_idx, neg_idx = r["pos"], r["neg"]

            def profile(pos=pos_idx, neg=neg_idx, zbb=zb, zii=zi) -> list[float]:
                out = []
                for th in THETAS:
                    c, s = math.cos(math.radians(th)), math.sin(math.radians(th))
                    sc = [c * zbb[k] + s * zii[k] for k in range(len(units))]
                    out.append(auc([sc[k] for k in pos], [sc[k] for k in neg]))
                return out

            obs = profile()
            N = math.comb(len(pos_idx) + len(neg_idx), len(pos_idx))
            obs_aucs = {
                "baseline": obs[0],
                "fold": auc([zscore(fold_all)[k] for k in pos_idx],
                            [zscore(fold_all)[k] for k in neg_idx]),
                "induced": obs[-1],
            }
            best_t = max(range(len(THETAS)), key=lambda t: obs[t])

            perm_max = []
            perm_at_best = []
            universe = pos_idx + neg_idx
            for combo in itertools.combinations(range(len(universe)), len(pos_idx)):
                sel = set(combo)
                p = [universe[i] for i in range(len(universe)) if i in sel]
                n_ = [universe[i] for i in range(len(universe)) if i not in sel]
                pp = profile(p, n_)
                perm_max.append(max(pp))
                perm_at_best.append(pp[best_t])
            p_fw = sum(1 for v in perm_max if v >= obs[best_t]) / N
            p_pt = sum(1 for v in perm_at_best if v >= obs[best_t]) / N

            ind_beats_both = obs_aucs["induced"] >= max(
                obs_aucs["baseline"], obs_aucs["fold"]
            ) - 1e-12
            summary["grid"][mod][rule_name] = {
                "auc_baseline": round(obs_aucs["baseline"], 4),
                "auc_fold": round(obs_aucs["fold"], 4),
                "auc_induced": round(obs_aucs["induced"], 4),
                "theta_at_max_deg": THETAS[best_t],
                "auc_at_max": round(obs[best_t], 4),
                "auc_span_over_grid": round(max(obs) - min(obs), 4),
                "R2_flat_profile": bool(max(obs) - min(obs) <= 0.10),
                "p_pointwise_at_max": round(p_pt, 4),
                "p_familywise_maxstat": round(p_fw, 4),
                "R4_significant_family_wise": bool(p_fw < 0.05),
                "induced_ge_both_rivals": bool(ind_beats_both),
            }
            for t, th in enumerate(THETAS):
                rows_out.append(
                    {
                        "module": mod,
                        "group_rule": rule_name,
                        "theta_deg": th,
                        "auc": round(obs[t], 4),
                    }
                )

    # R3 verdict across the whole grid
    violations = [
        f"{mod}/{rule}"
        for mod, d in summary["grid"].items()
        for rule, v in d.items()
        if not v["induced_ge_both_rivals"]
    ]
    n_cells = sum(len(d) for d in summary["grid"].values())
    summary["R3_verdict"] = {
        "n_cells_tested": n_cells,
        "n_violations": len(violations),
        "violations": violations,
        "supported": len(violations) == 0,
    }

    with open(os.path.join(OUT, "contrast_grid_v2.tsv"), "w", encoding="utf-8", newline="") as fh:
        cols = ["module", "group_rule", "theta_deg", "auc"]
        fh.write("\t".join(cols) + "\n")
        for r in rows_out:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(os.path.join(OUT, "contrast_profile_v2.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False)

    print("units:", summary["units"])
    print("groups:", {k: len(v) for k, v in groups.items()})
    print()
    for rule in rules:
        ps = summary["permutation_space"][rule]
        print(f"{rule:24s} n={ps['n_pos']}v{ps['n_neg']} assignments={ps['n_assignments']} "
              f"min_p={ps['min_achievable_p']} head<= {ps['max_headcount_for_p_below_0.05']}")
    print()
    hdr = f"{'module':24s} {'rule':24s} {'base':>6s} {'fold':>6s} {'induc':>6s} {'max':>6s} " \
          f"{'span':>6s} {'p_fw':>7s} {'ind>both':>8s}"
    print(hdr)
    for mod, d in summary["grid"].items():
        for rule, v in d.items():
            print(f"{mod:24s} {rule:24s} {v['auc_baseline']:6.3f} {v['auc_fold']:6.3f} "
                  f"{v['auc_induced']:6.3f} {v['auc_at_max']:6.3f} {v['auc_span_over_grid']:6.3f} "
                  f"{v['p_familywise_maxstat']:7.4f} {str(v['induced_ge_both_rivals']):>8s}")
    print()
    print("R3:", summary["R3_verdict"])
    print()
    for mod, d in summary["degeneracy"].items():
        print(f"degeneracy {mod:24s} rho={d['spearman_baseline_vs_induced_n11']:.3f} "
              f"slope={d['ols_slope_induced_on_baseline_n11']:.3f} R2={d['ols_r2_n11']:.3f} "
              f"sd(f)/sd(b)={d['sd_fold_over_sd_baseline_n11']:.3f} R1={d['R1_degenerate']}")
    print()
    for mod, d in summary["ordinal_all_units"].items():
        print(f"ordinal   {mod:24s} base={d['spearman_baseline_vs_grade']:+.3f} "
              f"fold={d['spearman_fold_vs_grade']:+.3f} induced={d['spearman_induced_vs_grade']:+.3f}")
    print()
    print("wrote:", os.path.join(OUT, "contrast_grid_v2.tsv"))
    print("wrote:", os.path.join(OUT, "contrast_profile_v2.json"))


if __name__ == "__main__":
    main()
