#!/usr/bin/env python
"""W14b: organ tropism with tissue fixed effects.

The pooled W14 model mixes two very different sources of variation:
  * between-tissue differences (what "tropism" means), and
  * within-tissue time courses (what the source paper analysed).
This script separates them:

  1. between-tissue: median viral load per tissue vs median score (Spearman)
  2. within-tissue:  viral load ~ score, with tissue fixed effects
  3. paired:         per animal x day, rank tissues by viral load and correlate
  4. monocyte check: does a marker-based myeloid score explain viral load
     within tissue (the source paper's own conclusion)?

Outputs: analysis/w14_organ_20260924/fixedeffects_*.tsv /.txt
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis" / "w14_organ_20260924"

FEATURES = ["S_AT_lectin", "S_AT_broad", "S_EM_all", "S_EM_protease",
            "S_RE_interferon", "S_ISG_control", "S_MONO_marker"]


def main() -> None:
    samples = pd.read_csv(OUT / "sample_table.tsv", sep="\t")
    work = samples.dropna(subset=["log10_viral"]).copy()
    work["tissue"] = work["tissue_short"].astype("category")
    work["day_num"] = pd.to_numeric(
        work["day"].str.replace("D", "", regex=False), errors="coerce")
    print(f"samples with viral load: {len(work)}; "
          f"tissues: {work['tissue'].nunique()}")

    # ---- 1. between-tissue ----
    rows = []
    for feat in FEATURES:
        block = work.groupby("tissue", observed=True).agg(
            viral=("log10_viral", "median"), score=(feat, "median"),
            n=("log10_viral", "size"))
        block = block[block["n"] >= 5]
        if len(block) < 5:
            continue
        rho, p = spearmanr(block["score"], block["viral"])
        rows.append({"level": "between_tissue", "feature": feat,
                     "n": len(block), "estimate": round(float(rho), 4),
                     "p_value": round(float(p), 4)})
    between = pd.DataFrame(rows)
    between.to_csv(OUT / "fixedeffects_between_tissue.tsv", sep="\t",
                   index=False)

    # ---- 2. within-tissue fixed effects ----
    rows = []
    for feat in FEATURES:
        data = work.dropna(subset=[feat, "day_num"])
        data = data.copy()
        data["z"] = (data[feat] - data[feat].mean()) / data[feat].std()
        try:
            fit = smf.ols("log10_viral ~ z + C(tissue) + day_num",
                          data=data).fit()
            rows.append({"level": "within_tissue_fe", "feature": feat,
                         "n": int(fit.nobs),
                         "estimate": round(float(fit.params["z"]), 4),
                         "p_value": round(float(fit.pvalues["z"]), 5),
                         "r2": round(float(fit.rsquared), 4)})
        except Exception as exc:
            rows.append({"level": "within_tissue_fe", "feature": feat,
                         "n": len(data), "estimate": np.nan,
                         "p_value": np.nan, "r2": np.nan,
                         "error": str(exc)})
    within = pd.DataFrame(rows)
    within.to_csv(OUT / "fixedeffects_within_tissue.tsv", sep="\t", index=False)

    # a joint model: entry + machinery + restriction + monocyte
    data = work.dropna(subset=FEATURES + ["day_num"]).copy()
    for feat in FEATURES:
        data["z_" + feat] = (data[feat] - data[feat].mean()) / data[feat].std()
    joint = smf.ols(
        "log10_viral ~ z_S_AT_lectin + z_S_EM_all + z_S_RE_interferon + "
        "z_S_MONO_marker + C(tissue) + day_num", data=data).fit()
    lines = ["joint within-tissue model (tissue fixed effects + day):",
             f"  n={int(joint.nobs)}  R2={joint.rsquared:.3f}  "
             f"adjR2={joint.rsquared_adj:.3f}"]
    for term in ("z_S_AT_lectin", "z_S_EM_all", "z_S_RE_interferon",
                 "z_S_MONO_marker"):
        lines.append(f"  {term:18s} beta={joint.params[term]: .3f} "
                     f"p={joint.pvalues[term]:.4f}")

    # ---- 3. paired across tissues within animal x day ----
    paired = []
    for (animal, day), block in work.groupby(["animal", "day"], observed=True):
        block = block.dropna(subset=["day_num"])
        block = block[block["day_num"] > 0]
        if len(block) < 6:
            continue
        for feat in FEATURES:
            rho, p = spearmanr(block[feat], block["log10_viral"])
            paired.append({"animal": animal, "day": day, "n_tissues": len(block),
                           "feature": feat, "rho": float(rho)})
    paired_frame = pd.DataFrame(paired)
    if not paired_frame.empty:
        summary = (paired_frame.groupby("feature")["rho"]
                   .agg(n="size", median="median", mean="mean")
                   .round(4).reset_index())
    else:
        summary = pd.DataFrame()
    paired_frame.to_csv(OUT / "fixedeffects_paired_per_animal.tsv", sep="\t",
                        index=False)
    summary.to_csv(OUT / "fixedeffects_paired_summary.tsv", sep="\t",
                   index=False)

    pd.set_option("display.width", 200)
    print("\n--- between tissue ---")
    print(between.to_string(index=False))
    print("\n--- within tissue (fixed effects) ---")
    print(within.to_string(index=False))
    print("\n--- joint ---")
    print("\n".join(lines))
    print("\n--- paired per animal x day ---")
    print(summary.to_string(index=False))
    (OUT / "fixedeffects_joint_model.txt").write_text("\n".join(lines),
                                                     encoding="utf-8")


if __name__ == "__main__":
    main()
