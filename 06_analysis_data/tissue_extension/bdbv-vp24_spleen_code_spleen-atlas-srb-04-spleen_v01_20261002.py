"""SRB step 4 - register the spleen panel on the frozen scale and place it (R10-R13).

Step 3 produced two things that must be kept apart: the *ordering* transfers from
HCL to the HPA panel (anchor Spearman 0.615), but the *absolute scale* does not -
the frozen transform over-estimates delta on the anchors by a constant offset.
This step quantifies that offset, propagates it, and reports where the spleen
compartment would sit, under three conventions:

  C1  frozen transform, uncorrected   (no fitted parameters; also the biased one)
  C2  affine anchor calibration       y = a + b * x fitted on the 20 anchors
  C3  leave-one-anchor-out calibration, to show how much the placement depends
      on any single anchor

Uncertainty on the spleen members themselves comes from a cell-level bootstrap
(2,000 resamples of the spleen cells, pseudobulk recomputed each time), which is
independent of the cross-resource issue.

Placement is reported against the 26 compartment medians of the frozen coverage
audit (analysis/b6_coverage_20260926/compartment_coverage.tsv), in which the
spleen currently has no member at all.

Criteria frozen before the numbers were read:
  R10  the anchor offset is reported as mean/median/SD of (transfer - HPA);
  R11  the calibrated spleen value is reported with the LOO range of the median
       over spleen members, and the calibration is declared usable only if that
       LOO range is narrower than the gap between adjacent compartment medians at
       the spleen's position;
  R12  the spleen is placed only as an interval, never as a single rank;
  R13  every member keeps its cell count, and members below 100 cells are flagged.

Writes raw/srb_placement.json, raw/srb_spleen_members.tsv, raw/srb_bootstrap.tsv.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys

import h5py
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926")
RAW = os.path.join(OUT, "raw")
B6 = os.path.join(ROOT, "analysis", "b6_coverage_20260926", "compartment_coverage.tsv")
H2RAW = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926", "raw")
HCL = os.path.join(H2RAW, "atlas_h5ad", "2adb1f8a_Construction_of_a_human_cell_landscape_at_single.h5ad")

SPLEEN_LOW, SPLEEN_HIGH = 191895, 207700
BLOCK = 20000
N_BOOT = 2000
SEED = 20260926
SIX = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]
MODS = [("IFNAR1", "IFNAR2"), ("IFNLR1", "IL10RB"), ("ISG15", "MX1")]


def components(cpm, stats=None):
    out = []
    for genes in MODS:
        vals = []
        for g in genes:
            lv = math.log10(cpm.get(g, 0.0) + 1.0)
            if stats is not None:
                lv = (lv - stats[g]["mean_log10_nCPM"]) / stats[g]["sd_log10_nCPM"]
            vals.append(lv)
        out.append(sum(vals) / len(vals))
    return out


def ols(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((v - mx) ** 2 for v in x)
    sxy = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    b = sxy / sxx if sxx else float("nan")
    return b, my - b * mx


def scan_spleen():
    sys.path.insert(0, os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926"))
    from b5_common import decode_h5_node

    with h5py.File(HCL, "r") as f:
        node = f["raw/X"]
        genes = [str(x) for x in decode_h5_node(f["raw/var"]["feature_name"])]
        pos = {}
        for k, g in enumerate(SIX):
            hit = [i for i, nm in enumerate(genes) if nm == g]
            if hit:
                pos[hit[0]] = k
        lab = decode_h5_node(f["obs"]["cell_type_ontology_term_id"]).astype(str)[SPLEEN_LOW:SPLEEN_HIGH + 1]
        lo, hi = SPLEEN_LOW, SPLEEN_HIGH + 1
        vals = np.zeros((hi - lo, len(SIX)), dtype=np.float64)
        lib = np.zeros(hi - lo, dtype=np.float64)
        data, indices, indptr = node["data"], node["indices"], node["indptr"]
        lut = np.full(int(node.attrs["shape"][1]), -1, dtype=np.int16)
        for gidx, col in pos.items():
            lut[gidx] = col
        for i0 in range(lo, hi, BLOCK):
            i1 = min(hi, i0 + BLOCK)
            ip = indptr[i0:i1 + 1]
            d0, d1 = int(ip[0]), int(ip[-1])
            if d1 <= d0:
                continue
            d = data[d0:d1].astype(np.float64)
            ix = indices[d0:d1]
            rel = np.asarray(ip - d0, dtype=np.int64)
            rows = np.repeat(np.arange(i1 - i0), np.diff(rel))
            lib[i0 - lo:i1 - lo] = np.bincount(rows, weights=d, minlength=i1 - i0)
            col_of = lut[ix]
            m = col_of >= 0
            if m.any():
                vals[i0 - lo + rows[m], col_of[m]] = d[m]
    return vals, lib, lab


def delta_for(mask, vals, lib, stats):
    sums = vals[mask].sum(axis=0)
    lsum = float(lib[mask].sum())
    if lsum <= 0:
        return None, None
    cpm = {g: 1e6 * sums[k] / lsum for k, g in enumerate(SIX)}
    return sum(components(cpm)) / 3.0, sum(components(cpm, stats)) / 3.0


def median(xs):
    s = sorted(xs)
    return s[len(s) // 2]


def main() -> None:
    stats = json.load(open(os.path.join(RAW, "hpa_gene_stats_verified.json"), encoding="utf-8"))["gene_stats"]
    members = [r for r in csv.DictReader(open(os.path.join(RAW, "hcl_spleen_pseudobulk.tsv"), encoding="utf-8"),
                                         delimiter="\t") if r["route"] == "cl"]
    anchors = list(csv.DictReader(open(os.path.join(RAW, "srb_anchor_check.tsv"), encoding="utf-8"), delimiter="\t"))

    ax = [float(a["delta_transfer_whole"]) for a in anchors]
    ay = [float(a["hpa_delta"]) for a in anchors]
    res = [ax[i] - ay[i] for i in range(len(ax))]
    mean_r = sum(res) / len(res)
    sd_r = math.sqrt(sum((v - mean_r) ** 2 for v in res) / (len(res) - 1))
    sres = sorted(res)
    print(f"R10 anchors n = {len(ax)}: offset (transfer - HPA) mean {mean_r:+.3f}, median {median(res):+.3f}, "
          f"SD {sd_r:.3f}, range {sres[0]:+.3f} .. {sres[-1]:+.3f}")

    v = [float(m["delta_transfer"]) for m in members]
    med_frozen = median(v)
    b, a = ols(ax, ay)
    cal = [a + b * x for x in v]
    med_cal = median(cal)
    print(f"spleen members n = {len(v)}: frozen-transform median {med_frozen:+.3f}; "
          f"affine-calibrated median (y = {a:+.3f} + {b:.3f}x) {med_cal:+.3f}")

    loo = []
    for i in range(len(ax)):
        xs = [ax[j] for j in range(len(ax)) if j != i]
        ys = [ay[j] for j in range(len(ax)) if j != i]
        bb, aa = ols(xs, ys)
        loo.append({"left_out": anchors[i]["hpa_cell_type"], "slope": bb, "intercept": aa,
                    "spleen_median": median([aa + bb * x for x in v])})
    vals_loo = sorted(x["spleen_median"] for x in loo)
    print(f"R11 LOO over anchors: calibrated spleen median ranges {vals_loo[0]:+.3f} .. {vals_loo[-1]:+.3f}")

    # stratified bootstrap: resample cells *within each label*, so the member set
    # (and hence the composition of the median) is held fixed and only sampling
    # noise in each member's pseudobulk is propagated.
    vals, lib, lab = scan_spleen()
    rng = np.random.default_rng(SEED)
    member_idx = [(u, np.where(lab == u)[0]) for u in sorted(set(lab.tolist())) if (lab == u).sum() >= 50]
    boot = []
    for _ in range(N_BOOT):
        ds = []
        for _u, rows_ in member_idx:
            pick = rows_[rng.integers(0, len(rows_), len(rows_))]
            _raw, d_tr = delta_for(np.ones(len(pick), dtype=bool), vals[pick], lib[pick], stats)
            if d_tr is not None:
                ds.append(d_tr)
        boot.append(median(ds))
    bs = sorted(boot)
    lo_b, hi_b = bs[int(.025 * len(bs))], bs[int(.975 * len(bs))]
    print(f"bootstrap ({len(bs)}/{N_BOOT} usable): spleen frozen-transform median {lo_b:+.3f} .. {hi_b:+.3f} (2.5-97.5%)")

    comps = [r for r in csv.DictReader(open(B6, encoding="utf-8"), delimiter="\t") if r["median_delta"]]
    comps.sort(key=lambda r: -float(r["median_delta"]))
    above = [c for c in comps if float(c["median_delta"]) > med_cal]
    below = [c for c in comps if float(c["median_delta"]) <= med_cal]
    gap = (min(float(c["median_delta"]) for c in above) - max(float(c["median_delta"]) for c in below)) if above and below else float("nan")
    print(f"R12 placement: calibrated median {med_cal:+.3f} lies above {below[0]['compartment'] if below else '-'} "
          f"({float(below[0]['median_delta']):+.3f}) and below {above[-1]['compartment'] if above else '-'} "
          f"({float(above[-1]['median_delta']):+.3f}); adjacent-median gap = {gap:.3f}; "
          f"LOO span {vals_loo[-1]-vals_loo[0]:.3f} -> calibration {'usable' if (vals_loo[-1]-vals_loo[0]) < gap else 'NOT usable'}")

    with open(os.path.join(RAW, "srb_spleen_members.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["label", "label_name", "n_cells", "delta_transfer", "delta_raw",
                                           "calibrated", "under_100_cells"], delimiter="\t")
        w.writeheader()
        for m, c in zip(members, cal):
            w.writerow({"label": m["label"], "label_name": m["label_name"], "n_cells": m["n_cells"],
                        "delta_transfer": m["delta_transfer"], "delta_raw": m["delta_raw"],
                        "calibrated": round(c, 5), "under_100_cells": int(m["n_cells"]) < 100})
    with open(os.path.join(RAW, "srb_bootstrap.tsv"), "w", encoding="utf-8", newline="") as fh:
        fh.write("boot_median_frozen\n")
        for x in bs:
            fh.write(f"{x:.5f}\n")

    summary = {
        "R10": {"n_anchors": len(ax), "offset_mean": mean_r, "offset_median": median(res), "offset_sd": sd_r,
                "offset_min": sres[0], "offset_max": sres[-1],
                "reading": "the frozen transform over-estimates delta on HCL by about this offset"},
        "R11": {"slope": b, "intercept": a, "spleen_median_frozen": med_frozen,
                "spleen_median_calibrated": med_cal, "loo_range": [vals_loo[0], vals_loo[-1]], "loo": loo},
        "R12": {"compartment_above": above[-1]["compartment"] if above else None,
                "compartment_below": below[0]["compartment"] if below else None,
                "adjacent_gap": gap, "n_compartment_medians": len(comps),
                "calibration_usable": bool((vals_loo[-1] - vals_loo[0]) < gap)},
        "bootstrap": {"n": len(bs), "median_2p5": lo_b, "median_97p5": hi_b},
        "members": [{"label": m["label"], "label_name": m["label_name"], "n_cells": int(m["n_cells"])} for m in members],
    }
    json.dump(summary, open(os.path.join(RAW, "srb_placement.json"), "w", encoding="utf-8"), indent=1)
    print("wrote raw/srb_placement.json, raw/srb_spleen_members.tsv, raw/srb_bootstrap.tsv")


if __name__ == "__main__":
    main()
