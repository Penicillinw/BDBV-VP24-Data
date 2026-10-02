"""SRB-VERIFY step 5b: exact permutation p for the six-label TS/HCL comparison.

sv_05 recovered rho = 0.485714 exactly but a permutation p of 0.3632 against
SRB's 0.3524.  With n = 6 there are only 6! = 720 distinct permutations, so the
exact p can be enumerated instead of sampled; that removes the Monte Carlo
question entirely and gives the number the manuscript should quote.

Writes raw/sv_ts_exact_perm.json
"""

from __future__ import annotations

import csv
import itertools
import json
import os

from sv_04_null_frozen_scale import frozen_delta, spearman
from sv_03_verify import SIX

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
RAW = os.path.join(ROOT, "analysis", "srb_verify_20260926", "raw")


def main() -> None:
    stats = {}
    with open(os.path.join(RAW, "sv_hpa_gene_background.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stats[r["gene"]] = (float(r["mean_log10_nCPM"]), float(r["sd_log10_nCPM"]))
    mus = [stats[g][0] for g in SIX]
    sds = [stats[g][1] for g in SIX]
    ts = {r["label"]: r for r in csv.DictReader(
        open(os.path.join(RAW, "sv_ts_pseudobulk.tsv"), encoding="utf-8"), delimiter="\t")}
    hcl = {r["label"]: r for r in csv.DictReader(
        open(os.path.join(RAW, "sv_hcl_spleen_pseudobulk.tsv"), encoding="utf-8"), delimiter="\t")}
    # the HCL side must meet the same 50-cell floor as the members (SRB's rule);
    # without it a 22-cell HCL monocyte label would enter the intersection
    hcl = {u: r for u, r in hcl.items() if int(r["n_cells"]) >= 50}
    shared = sorted(set(ts) & set(hcl))
    xs = [frozen_delta([float(ts[u]["cpm_" + g]) for g in SIX], mus, sds) for u in shared]
    ys = [frozen_delta([float(hcl[u][g]) for g in SIX], mus, sds) for u in shared]
    obs = spearman(xs, ys)
    perms = list(itertools.permutations(range(len(shared))))
    rhos = [spearman(xs, [ys[i] for i in p]) for p in perms]
    n = len(rhos)
    ge = sum(1 for v in rhos if abs(v) >= abs(obs))
    n_ge_one_sided = sum(1 for v in rhos if v >= obs)
    out = {"n_matched": len(shared), "matched": shared, "observed_spearman": obs,
           "n_permutations_exact": n, "n_abs_ge_observed": ge,
           "exact_two_sided_p": ge / n, "exact_one_sided_p": n_ge_one_sided / n,
           "distinct_rho_values": len(set(round(v, 12) for v in rhos))}
    with open(os.path.join(RAW, "sv_ts_exact_perm.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
