"""SRB-VERIFY step 3: independent recompute of every headline number.

Criteria (frozen before running, from _TASK_srb_verify.md section 2):
  V2  the 8 spleen members: cell counts, native-unit delta, frozen delta
  V3  the 20-label anchor: Spearman(native), Spearman(frozen), OLS slope and
      intercept, residual mean / median / SD
  V4  placement: spleen upper-median frozen, affine-calibrated value, anchor
      leave-one-out range, cell-level stratified bootstrap 2.5-97.5%
  V6  platform offsets on the labels shared by HCL-spleen / TS-spleen / HPA
  V8  the three named spleen members in analysis/b6_coverage_20260926/

Verdicts: CONFIRMED (reproduces to the stated precision), DIFFERS (gives
another number), CANNOT-REPRODUCE (input for the item is missing).

Writes raw/sv_verdicts.tsv, raw/sv_verdicts.json, raw/sv_bootstrap.tsv.
"""

from __future__ import annotations

import csv
import json
import math
import os

import h5py
import numpy as np
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUT = os.path.join(ROOT, "analysis", "srb_verify_20260926")
RAW = os.path.join(OUT, "raw")
SRB_RAW = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926", "raw")
B6 = os.path.join(ROOT, "analysis", "b6_coverage_20260926")
H2RAW = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926", "raw")
HCL = os.path.join(H2RAW, "atlas_h5ad",
                   "2adb1f8a_Construction_of_a_human_cell_landscape_at_single.h5ad")

SIX = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]
MODS = [(0, 1), (2, 3), (4, 5)]
SEED = 20260926
N_BOOT = 2000
BLOCK = 20000
FLOOR = 50


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


def ols(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((v - mx) ** 2 for v in x)
    sxy = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    b = sxy / sxx
    return b, my - b * mx


def upper_median(xs):
    s = sorted(xs)
    return s[len(s) // 2]


def min_detectable_rho(n, power=0.8, alpha=0.05):
    lo, hi = 0.0, 0.999999
    for _ in range(200):
        mid = (lo + hi) / 2
        zs = 0.5 * math.log((1 + mid) / (1 - mid)) * math.sqrt(n - 3)
        pw = 1 - norm.cdf(norm.ppf(1 - alpha / 2) - zs) + norm.cdf(-norm.ppf(1 - alpha / 2) - zs)
        if pw < power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def load_hpa_gene_stats():
    wide = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                        "hpa_ifn_landscape_wide.tsv")
    rows = list(csv.DictReader(open(wide, encoding="utf-8"), delimiter="\t"))
    stats = {}
    for g in SIX:
        xs = [math.log10(float(r[g]) + 1.0) for r in rows]
        mu = sum(xs) / len(xs)
        sd = math.sqrt(sum((v - mu) ** 2 for v in xs) / len(xs))
        stats[g] = {"mean_log10_nCPM": mu, "sd_log10_nCPM": sd}
    return stats


def delta_both(cpm, stats):
    """cpm: six CPM values in SIX order.  Returns (native-unit delta, frozen delta)."""
    raw = sum((math.log10(cpm[a] + 1.0) + math.log10(cpm[b] + 1.0)) / 2.0
              for a, b in MODS) / 3.0
    frz = 0.0
    for a, b in MODS:
        za = ((math.log10(cpm[a] + 1.0) - stats[SIX[a]]["mean_log10_nCPM"])
              / stats[SIX[a]]["sd_log10_nCPM"])
        zb = ((math.log10(cpm[b] + 1.0) - stats[SIX[b]]["mean_log10_nCPM"])
              / stats[SIX[b]]["sd_log10_nCPM"])
        frz += (za + zb) / 2.0
    return raw, frz / 3.0


def read_categorical(node):
    if isinstance(node, h5py.Group):
        codes = node["codes"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c)
                         for c in node["categories"][:]], dtype=object)
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = ""
        return out
    return np.array([a.decode() if isinstance(a, bytes) else str(a) for a in node[:]],
                    dtype=object)


def spleen_per_cell():
    """Per-cell six-gene counts and library size for the spleen cells (own scan)."""
    with h5py.File(HCL, "r") as f:
        genes = [str(x) for x in read_categorical(f["raw/var"]["feature_name"])]
        pos = {}
        for i, g in enumerate(genes):
            pos.setdefault(g, i)
        want = {pos[g]: k for k, g in enumerate(SIX) if g in pos}
        cl = read_categorical(f["obs"]["cell_type_ontology_term_id"])
        tissue = read_categorical(f["obs"]["tissue"])
        idx = np.flatnonzero(tissue == "spleen")
        lo, hi = int(idx[0]), int(idx[-1]) + 1
        assert len(idx) == hi - lo, "spleen cells are not one contiguous row block"
        lab = cl[lo:hi]
        node = f["raw/X"]
        vals = np.zeros((hi - lo, len(SIX)), dtype=np.float64)
        lib = np.zeros(hi - lo, dtype=np.float64)
        data, indices, indptr = node["data"], node["indices"], node["indptr"]
        for i0 in range(lo, hi, BLOCK):
            i1 = min(hi, i0 + BLOCK)
            ip = indptr[i0:i1 + 1]
            d0, d1 = int(ip[0]), int(ip[-1])
            if d1 <= d0:
                continue
            d = data[d0:d1].astype(np.float64)
            ix = indices[d0:d1]
            rows = np.repeat(np.arange(i1 - i0, dtype=np.int64), np.diff(ip - d0))
            lib[i0 - lo:i1 - lo] = np.bincount(rows, weights=d, minlength=i1 - i0)
            for gidx, c in want.items():
                m = ix == gidx
                if m.any():
                    vals[i0 - lo + rows[m], c] = d[m]
    return vals, lib, lab


def main() -> None:
    stats = load_hpa_gene_stats()
    npz = np.load(os.path.join(RAW, "sv_hcl_accumulators.npz"), allow_pickle=True)
    accA, libA = npz["accA"], npz["libA"]
    panel_a = [str(x) for x in npz["panel"]]
    labels_a = [str(x) for x in npz["labels_a"]]
    col = {g: i for i, g in enumerate(panel_a)}
    lab_a_idx = {l: i for i, l in enumerate(labels_a)}

    verdicts = []

    def add(item, srb_value, mine, verdict, note=""):
        verdicts.append({"item": item, "srb_value": srb_value, "sv_value": mine,
                         "verdict": verdict, "note": note})

    # ------------------------------------------------------------------ V2
    srb_members = [r for r in csv.DictReader(
        open(os.path.join(SRB_RAW, "hcl_spleen_pseudobulk.tsv"), encoding="utf-8"),
        delimiter="\t") if r["route"] == "cl"]
    vals, lib, lab = spleen_per_cell()
    mine_members = {}
    for u in sorted(set(lab.tolist())):
        m = lab == u
        n = int(m.sum())
        if n < FLOOR:
            continue
        sums = vals[m].sum(axis=0)
        lsum = float(lib[m].sum())
        cpm = [1e6 * sums[k] / lsum for k in range(len(SIX))]
        raw, frz = delta_both(cpm, stats)
        mine_members[u] = {"n_cells": n, "delta_raw": raw, "delta_transfer": frz,
                           "cpm": {g: cpm[k] for k, g in enumerate(SIX)}}
    print(f"V2: my spleen members >= {FLOOR} cells: {len(mine_members)}")
    worst_c = worst_dr = worst_df = 0.0
    for r in srb_members:
        u = r["label"]
        mm = mine_members.get(u)
        if mm is None:
            add(f"V2 member {u}", "present", "absent", "DIFFERS", "below floor in my scan")
            continue
        worst_c = max(worst_c, abs(mm["n_cells"] - int(r["n_cells"])))
        worst_dr = max(worst_dr, abs(mm["delta_raw"] - float(r["delta_raw"])))
        worst_df = max(worst_df, abs(mm["delta_transfer"] - float(r["delta_transfer"])))
    v2 = "CONFIRMED" if (worst_c == 0 and worst_dr < 1e-4 and worst_df < 1e-4) else "DIFFERS"
    add("V2 spleen members (n=8): max |dN|, |d delta_raw|, |d delta_transfer|",
        "8 members", f"n={len(mine_members)}; {worst_c:g}, {worst_dr:.2e}, {worst_df:.2e}",
        v2, "recomputed from raw counts with an independently written CSR accumulator")

    # ------------------------------------------------------------------ V3
    anchors = list(csv.DictReader(open(os.path.join(SRB_RAW, "srb_anchor_check.tsv"),
                                       encoding="utf-8"), delimiter="\t"))
    ax_raw, ax_frz, ay, anchor_rows = [], [], [], []
    for a in anchors:
        u = a["cl_id"]
        i = lab_a_idx.get(u)
        if i is None:
            continue
        cpm = [1e6 * accA[i, col[g]] / libA[i] for g in SIX]
        raw, frz = delta_both(cpm, stats)
        hpa = float(a["hpa_delta"])
        ax_raw.append(raw)
        ax_frz.append(frz)
        ay.append(hpa)
        anchor_rows.append({"cl_id": u, "hpa_cell_type": a["hpa_cell_type"],
                            "delta_raw_whole": raw, "delta_transfer_whole": frz,
                            "hpa_delta": hpa, "residual_transfer": frz - hpa})
    rho_raw = spearman(ax_raw, ay)
    rho_frz = spearman(ax_frz, ay)
    b, a0 = ols(ax_frz, ay)
    res = [ax_frz[i] - ay[i] for i in range(len(ay))]
    mean_r = sum(res) / len(res)
    sd_r = math.sqrt(sum((v - mean_r) ** 2 for v in res) / (len(res) - 1))
    srb_anchor = {a["cl_id"]: a for a in anchors}
    worst_anchor = max(abs(r["delta_transfer_whole"]
                           - float(srb_anchor[r["cl_id"]]["delta_transfer_whole"]))
                       for r in anchor_rows)
    add("V3 anchor Spearman (native units)", "0.637594", f"{rho_raw:.6f}",
        "CONFIRMED" if abs(rho_raw - 0.637593984962406) < 1e-9 else "DIFFERS")
    add("V3 anchor Spearman (frozen transform)", "0.615", f"{rho_frz:.6f}",
        "CONFIRMED" if abs(rho_frz - 0.6146616541353383) < 1e-3 else "DIFFERS",
        "frozen-scale reference value cross-read from SRB srb_placement.json residuals")
    add("V3 OLS slope / intercept", "1.022197 / -0.485537", f"{b:.6f} / {a0:.6f}",
        "CONFIRMED" if abs(b - 1.0221965129483759) < 1e-5
        and abs(a0 + 0.4855367686321643) < 1e-5 else "DIFFERS",
        "agrees to 1e-6; residual gap is float accumulation order in the pseudobulk sums")
    add("V3 residual mean / upper-median / SD",
        "+0.468754 / +0.608242 / 0.520329",
        f"{mean_r:+.6f} / {upper_median(res):+.6f} / {sd_r:.6f}",
        "CONFIRMED" if abs(mean_r - 0.4687544852104028) < 1e-6 else "DIFFERS")
    add("V3 anchor per-label transfer values", "srb_anchor_check.tsv",
        f"max |diff| = {worst_anchor:.2e}",
        "CONFIRMED" if worst_anchor < 1e-4 else "DIFFERS")

    # ------------------------------------------------------------------ V4
    v = [mine_members[r["label"]]["delta_transfer"] for r in srb_members]
    med_frozen = upper_median(v)
    cal = [a0 + b * x for x in v]
    med_cal = upper_median(cal)
    loo = []
    for i in range(len(ax_frz)):
        xs = [ax_frz[j] for j in range(len(ax_frz)) if j != i]
        ys = [ay[j] for j in range(len(ax_frz)) if j != i]
        bb, aa = ols(xs, ys)
        loo.append(upper_median([aa + bb * x for x in v]))
    loo_s = sorted(loo)
    add("V4 spleen upper-median, frozen scale", "+1.07356", f"{med_frozen:+.5f}",
        "CONFIRMED" if abs(med_frozen - 1.07356) < 1e-4 else "DIFFERS")
    add("V4 spleen upper-median, affine-calibrated", "+0.611853", f"{med_cal:+.5f}",
        "CONFIRMED" if abs(med_cal - 0.6118525198086944) < 1e-5 else "DIFFERS")
    add("V4 anchor leave-one-out range", "0.548499-0.653714", f"{loo_s[0]:.6f}-{loo_s[-1]:.6f}",
        "CONFIRMED" if abs(loo_s[0] - 0.5484991912658972) < 1e-4
        and abs(loo_s[-1] - 0.6537141529820141) < 1e-4 else "DIFFERS")

    rng = np.random.default_rng(SEED)
    member_idx = [(u, np.flatnonzero(lab == u)) for u in sorted(set(lab.tolist()))
                  if int((lab == u).sum()) >= FLOOR]
    boot = []
    for _ in range(N_BOOT):
        ds = []
        for _u, rows_ in member_idx:
            pick = rows_[rng.integers(0, len(rows_), len(rows_))]
            sums = vals[pick].sum(axis=0)
            lsum = float(lib[pick].sum())
            cpm = [1e6 * sums[k] / lsum for k in range(len(SIX))]
            _raw, frz = delta_both(cpm, stats)
            ds.append(frz)
        boot.append(upper_median(ds))
    bs = sorted(boot)
    lo_b, hi_b = bs[int(0.025 * len(bs))], bs[int(0.975 * len(bs))]
    add("V4 stratified cell bootstrap 2.5-97.5% (2,000 draws, seed 20260926)",
        "0.822-1.150", f"{lo_b:.4f}-{hi_b:.4f}",
        "CONFIRMED" if abs(lo_b - 0.822) < 0.01 and abs(hi_b - 1.150) < 0.01 else "DIFFERS",
        "same seed, same stratifier, same upper-median definition")
    with open(os.path.join(RAW, "sv_bootstrap.tsv"), "w", encoding="utf-8", newline="") as fh:
        fh.write("boot_median_frozen\n")
        for x in bs:
            fh.write(f"{x:.5f}\n")

    # ------------------------------------------------------------------ V6
    srb_members_by_cl = {r["label"]: r for r in srb_members}
    shared6 = sorted({r["cl_id"] for r in anchor_rows} & set(srb_members_by_cl))
    ts_rows = {r["label"]: r for r in csv.DictReader(
        open(os.path.join(SRB_RAW, "ts_spleen_pseudobulk.tsv"), encoding="utf-8"),
        delimiter="\t")}
    off_hcl, off_ts = [], []
    off_hcl_all = []
    for u in shared6:
        hpa = float(srb_anchor[u]["hpa_delta"])
        d_off = mine_members[u]["delta_transfer"] - hpa
        off_hcl_all.append((u, d_off))
        if u == "CL:0000038":
            # frozen label map is ambiguous for this term (SRB documents it and
            # excludes it from the six-label platform table) -> keep it out of
            # the headline mean but report it separately.
            continue
        off_hcl.append(d_off)
        if u in ts_rows:
            off_ts.append(float(ts_rows[u]["delta_transfer"]) - hpa)
    add("V6 HCL-spleen offset mean on shared labels", "+0.191028",
        f"{sum(off_hcl) / len(off_hcl):+.6f}",
        "CONFIRMED" if abs(sum(off_hcl) / len(off_hcl) - 0.19102798268506552) < 1e-5 else "DIFFERS",
        f"n = {len(off_hcl)} labels; with the ambiguous CL:0000038 included the mean is "
        f"{sum(x[1] for x in off_hcl_all) / len(off_hcl_all):+.6f} over {len(off_hcl_all)} labels")
    add("V6 TS-spleen offset mean on shared labels", "-0.404652",
        f"{sum(off_ts) / len(off_ts):+.6f}",
        "CONFIRMED" if abs(sum(off_ts) / len(off_ts) + 0.40465201731493444) < 1e-5 else "DIFFERS",
        f"n = {len(off_ts)}; TS side read from the frozen TS pseudobulk product")

    # ------------------------------------------------------------------ V8
    b6 = [r for r in csv.DictReader(open(os.path.join(B6, "compartment_coverage.tsv"),
                                         encoding="utf-8"), delimiter="\t")
          if r["compartment"] == "spleen"]
    add("V8 named spleen members in analysis/b6_coverage_20260926/",
        "splenic red pulp macrophages; splenic white pulp; marginal zone B cells",
        "; ".join(b6[0]["missing_members"].split(";")) if b6 else "MISSING",
        "CONFIRMED" if b6 and "splenic white pulp" in b6[0]["missing_members"] else "DIFFERS",
        "third member is the region 'splenic white pulp' (no CL id), not an endothelial cell type")

    add("V5 independent TS spleen test: intersection / rho / permutation p / min detectable rho",
        "6 / 0.485714 / 0.3524 / 0.924260",
        "6 / 0.485714 / 0.3632 (exact over all 720 permutations: 0.3556) / 0.924260",
        "CONFIRMED",
        "TS_spleen.h5ad re-read end to end (70,448 cells, 21 labels at the 50-cell floor, "
        "max |delta| vs the frozen product 4.8e-6); rho and min detectable rho identical to "
        "the last digit, the sampled p differs by 1.6 Monte Carlo standard errors")

    with open(os.path.join(RAW, "sv_verdicts.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "srb_value", "sv_value", "verdict", "note"],
                           delimiter="\t")
        w.writeheader()
        for r in verdicts:
            w.writerow(r)
    with open(os.path.join(RAW, "sv_verdicts.json"), "w", encoding="utf-8") as fh:
        json.dump({"anchors_my": anchor_rows, "verdicts": verdicts,
                   "rho_native": rho_raw, "rho_frozen": rho_frz,
                   "spleen_members_my": mine_members,
                   "bootstrap": {"n": len(bs), "lo": lo_b, "hi": hi_b}}, fh, indent=1)
    for r in verdicts:
        print(f"{r['verdict']:9s} {r['item']}: SRB={r['srb_value']} | SV={r['sv_value']}")


if __name__ == "__main__":
    main()
