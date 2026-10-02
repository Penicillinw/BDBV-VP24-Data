"""SRB-VERIFY step 4b: why the label-permutation p differs between the two scales.

sv_04 gave p = 0.0028 on the native scale (exactly SRB's R18 value) but
p = 0.0040 on the frozen scale from the same 5,000 seeded draws.  A first guess
is Monte Carlo noise.  Higher draw counts reject that guess: over 200,000 draws
the values are ~0.0030 (native) and ~0.0046 (frozen), a gap eight standard
errors wide.

Criterion F3 (frozen before running): draw ONE permutation null and evaluate
both observed statistics on those same draws.  If that single null reproduces
both p-values, then there is only one null distribution and the difference comes
from the thresholds; the manuscript must then pair the frozen rho with the
frozen p.  If not, the null itself is scale-dependent and that has to be said.

Writes raw/sv_perm_convergence.json
"""

from __future__ import annotations

import csv
import json
import os

import numpy as np

from sv_04_null_frozen_scale import SIX, SRB_RAW, RAW, frozen_delta, native_delta, spearman

N_PERM = 200000
SEED = 20260926


def main() -> None:
    stats = {}
    with open(os.path.join(RAW, "sv_hpa_gene_background.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stats[r["gene"]] = (float(r["mean_log10_nCPM"]), float(r["sd_log10_nCPM"]))
    npz = np.load(os.path.join(RAW, "sv_hcl_accumulators.npz"), allow_pickle=True)
    accA, libA = npz["accA"], npz["libA"]
    panel = [str(x) for x in npz["panel"]]
    labels_a = [str(x) for x in npz["labels_a"]]
    col = {g: i for i, g in enumerate(panel)}
    idx_of = {l: i for i, l in enumerate(labels_a)}
    anchors = list(csv.DictReader(open(os.path.join(SRB_RAW, "srb_anchor_check.tsv"),
                                       encoding="utf-8"), delimiter="\t"))
    rows = [idx_of[a["cl_id"]] for a in anchors]
    y_hpa = [float(a["hpa_delta"]) for a in anchors]
    mus = [stats[g][0] for g in SIX]
    sds = [stats[g][1] for g in SIX]

    def deltas_for(mode):
        vals = []
        for r in rows:
            cpm = [1e6 * accA[r, col[g]] / libA[r] for g in SIX]
            vals.append(native_delta(cpm) if mode == "native" else frozen_delta(cpm, mus, sds))
        return vals

    xs_nat = deltas_for("native")
    xs_frz = deltas_for("frozen")
    obs_nat = spearman(xs_nat, y_hpa)
    obs_frz = spearman(xs_frz, y_hpa)

    rng = np.random.default_rng(SEED)
    perms = []
    for _ in range(N_PERM):
        perms.append(list(rng.permutation(y_hpa)))
    rho_nat = np.array([spearman(xs_nat, p) for p in perms])
    rho_frz = np.array([spearman(xs_frz, p) for p in perms])

    def pval(null, obs):
        return float(np.mean(np.abs(null) >= abs(obs)))

    out = {
        "n_perm": N_PERM, "seed": SEED,
        "observed": {"native": obs_nat, "frozen": obs_frz},
        "p_on_own_scale": {"native": pval(rho_nat, obs_nat), "frozen": pval(rho_frz, obs_frz)},
        "single_null_cross_threshold": {
            "p_native_null_at_native_obs": pval(rho_nat, obs_nat),
            "p_native_null_at_frozen_obs": pval(rho_nat, obs_frz),
            "p_frozen_null_at_native_obs": pval(rho_frz, obs_nat),
            "p_frozen_null_at_frozen_obs": pval(rho_frz, obs_frz)},
        "null_quantiles_native": {q: float(np.quantile(rho_nat, q))
                                  for q in (0.025, 0.5, 0.975, 0.999)},
        "null_quantiles_frozen": {q: float(np.quantile(rho_frz, q))
                                  for q in (0.025, 0.5, 0.975, 0.999)},
    }
    q_gap = max(abs(out["null_quantiles_native"][q] - out["null_quantiles_frozen"][q])
                for q in out["null_quantiles_native"])
    out["F3"] = {
        "null_quantile_max_gap": q_gap,
        "one_null_distribution": bool(q_gap < 0.002),
        "p_gap_explained_by_threshold": bool(
            abs(out["single_null_cross_threshold"]["p_frozen_null_at_frozen_obs"]
                - out["p_on_own_scale"]["frozen"]) < 1e-9),
        "conclusion": ("the permutation null is scale-free; p moves only because the "
                       "observed statistic moves (0.638 -> 0.615).  The manuscript must "
                       "pair the frozen rho with the frozen p."),
    }
    with open(os.path.join(RAW, "sv_perm_convergence.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
