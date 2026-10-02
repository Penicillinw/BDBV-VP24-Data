"""SRB-VERIFY step 4: the null distributions on the FROZEN scale (task item 2.7).

Why this exists
    SRB's R18 reported the anchor agreement as the *native-unit* Spearman
    (0.6376) but the manuscript sentence quotes the *frozen-transform* Spearman
    (0.615).  A p-value computed on one scale cannot be printed next to a
    correlation computed on the other.  This script recomputes both nulls on the
    frozen scale, and reproduces the native-scale null as a positive control on
    the same code path.

Pre-frozen verdict rule (item 2.7)
    F1  control: this pipeline must reproduce SRB's native-scale null
        (observed 0.6376; 0/200 random modules reach it; permutation p 0.0028).
    F2  the frozen-scale null must be reported as its own set of numbers
        (target statistic, fraction of random modules reached, permutation p).

Writes raw/sv_null_scale.json
"""

from __future__ import annotations

import csv
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
RAW = os.path.join(ROOT, "analysis", "srb_verify_20260926", "raw")
SRB_RAW = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926", "raw")

SIX = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]
MODS = [(0, 1), (2, 3), (4, 5)]
N_PERM = 5000
SEED = 20260926


def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(a, b):
    ra, rb = rank(list(a)), rank(list(b))
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    da = math.sqrt(sum((v - ma) ** 2 for v in ra))
    db = math.sqrt(sum((v - mb) ** 2 for v in rb))
    return num / (da * db) if da and db else float("nan")


def native_delta(cpm):
    return sum((math.log10(cpm[a] + 1.0) + math.log10(cpm[b] + 1.0)) / 2.0
               for a, b in MODS) / 3.0


def frozen_delta(cpm, mus, sds):
    tot = 0.0
    for a, b in MODS:
        za = (math.log10(cpm[a] + 1.0) - mus[a]) / sds[a]
        zb = (math.log10(cpm[b] + 1.0) - mus[b]) / sds[b]
        tot += (za + zb) / 2.0
    return tot / 3.0


def pct(values, q):
    s = sorted(values)
    return s[min(len(s) - 1, int(q * len(s)))]


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

    def deltas_for(genes, mode):
        mus = [stats[g][0] for g in genes]
        sds = [stats[g][1] for g in genes]
        vals = []
        for r in rows:
            cpm = [1e6 * accA[r, col[g]] / libA[r] for g in genes]
            vals.append(native_delta(cpm) if mode == "native" else frozen_delta(cpm, mus, sds))
        return vals

    modules = json.load(open(os.path.join(RAW, "sv_modules.json"), encoding="utf-8"))["modules"]
    usable = [m for m in modules if all(g in col and g in stats for g in m)]
    n_skipped = len(modules) - len(usable)

    out = {"n_modules": len(modules), "n_modules_usable": len(usable),
           "n_modules_skipped_missing_gene": n_skipped, "n_anchors": len(rows)}

    for mode in ("native", "frozen"):
        obs = spearman(deltas_for(SIX, mode), y_hpa)
        null = [spearman(deltas_for(m, mode), y_hpa) for m in usable]
        frac = sum(1 for v in null if v >= obs) / len(null)
        rng = np.random.default_rng(SEED)
        xs = deltas_for(SIX, mode)
        perm = []
        for _ in range(N_PERM):
            perm.append(spearman(xs, list(rng.permutation(y_hpa))))
        p_perm = float(np.mean([abs(v) >= abs(obs) for v in perm]))
        out["R18_" + mode] = {
            "observed_anchor_spearman": obs,
            "null_median": pct(null, 0.5), "null_p95": pct(null, 0.95), "null_max": max(null),
            "fraction_random_modules_ge_observed": frac,
            "label_permutation_p": p_perm, "n_perm": N_PERM,
            "specific_at_p95": bool(obs > pct(null, 0.95)),
        }
        print(f"R18 [{mode:6s}] observed {obs:.6f}; random-module null median {pct(null, 0.5):.3f}, "
              f"p95 {pct(null, 0.95):.3f}, max {max(null):.3f}; fraction >= observed {frac:.3f}; "
              f"label-permutation p {p_perm:.4f}")

    # formal check that the permutation null is scale-invariant (both vectors have no ties)
    n = out["R18_native"]["label_permutation_p"]
    f = out["R18_frozen"]["label_permutation_p"]
    out["F1_control_reproduces_SRB_native"] = bool(
        abs(out["R18_native"]["observed_anchor_spearman"] - 0.637593984962406) < 1e-9
        and out["R18_native"]["fraction_random_modules_ge_observed"] == 0.0
        and abs(n - 0.0028) < 0.001)
    out["permutation_null_scale_invariant"] = bool(abs(n - f) < 1e-12)
    out["tie_check"] = {"hpa_delta_unique": len(set(y_hpa)) == len(y_hpa),
                        "hcl_anchor_delta_unique_native": len(set(deltas_for(SIX, "native"))) == len(rows),
                        "hcl_anchor_delta_unique_frozen": len(set(deltas_for(SIX, "frozen"))) == len(rows)}

    with open(os.path.join(RAW, "sv_null_scale.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if not k.startswith("R18_")}, indent=1))


if __name__ == "__main__":
    main()
