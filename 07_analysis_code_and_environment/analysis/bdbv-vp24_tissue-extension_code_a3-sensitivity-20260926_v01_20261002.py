#!/usr/bin/env python
"""A3 sensitivity / second instrument (EXPLORATORY -- does not touch the
pre-registered verdict of tools/a3_deconvolution_20260926.py).

Why this file exists
--------------------
The pre-registered run (`tools/a3_deconvolution_20260926.py`,
pre-registration `report/A3_预注册_组织构成去卷积_20260926.md`) returned
verdict = "undecidable: A2 anchor failed" (rho of the inferred immune fraction
against the bulk immune-marker index was 0.405, pass line 0.5), while its A1
anchor passed 6/8. Two further readings are reported here so that a *method*
failure cannot be mistaken for a *negative biological* result:

  instrument 2 -- the anchor failures (skin -> adipocytes 0.54, sex organ ->
      smooth muscle 0.26) look like the classic NNLS "sink class" pathology:
      selecting markers by share alone lets a broadly expressed class absorb
      mass. Instrument 2 selects markers by a *specificity margin* over the
      second-best class instead, and is a completely independent composition
      estimate over the same 42 pre-registered classes.

  marker-based immune content -- a composition proxy that needs no
      deconvolution at all: the share of a sample's total CPM carried by
      immune marker genes (PTPRC/CD3E/MS4A1/NKG7/LYZ/CD68/CD14). Used to
      (i) cross-check the inferred immune fraction, and (ii) recompute E_adj
      without any NNLS step.

Neither reading enters the pre-registered verdict (pre-registration §4.5).
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis" / "a3_deconvolution_20260926"
A3_TOOL = ROOT / "tools" / "a3_deconvolution_20260926.py"

SEED = 20260926
TOP_MARKERS = 200
MIN_NCPM = 10.0


def load_a3():
    spec = importlib.util.spec_from_file_location("a3_tool", A3_TOOL)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def margin_markers(prof: pd.DataFrame, n_top: int = TOP_MARKERS) -> list[str]:
    """top-n per class by specificity margin over the second-best class."""
    arr = prof.to_numpy(dtype=float)
    keep: set[str] = set()
    genes = np.asarray(prof.index)
    for j, cls in enumerate(prof.columns):
        col = arr[:, j]
        others = np.delete(arr, j, axis=1)
        second = np.where(others.shape[1] > 0, others.max(axis=1), 0.0)
        margin = col / (second + 1.0)
        ok = col >= MIN_NCPM
        order = np.argsort(-np.where(ok, margin, -np.inf))[:n_top]
        keep.update(genes[o] for o in order if ok[o])
    return sorted(keep)


def main() -> None:
    a3 = load_a3()
    log = a3.log
    g8 = a3.load_g8()
    log("\n" + "=" * 78)
    log("A3 SENSITIVITY (exploratory; separate from the pre-registered verdict)")
    log("=" * 78)

    wide, _ = a3.load_hpa_reference()
    prof_ncpm, missing = a3.class_profiles(wide)
    meta = g8.load_sample_meta()
    z, depth, logcpm = g8.load_expression()
    meta = meta[meta["column"].isin(z.columns)].copy()
    meta["genes_detected"] = meta["column"].map(depth)
    meta["ISG_core"] = [g8.panel_score(z, g8.ISG_CORE, c) for c in meta["column"]]
    cpm = (2.0 ** logcpm) - 1.0
    tissues8 = list(a3.ANCHOR_TOP1)
    all8 = meta[meta["tissue"].isin(tissues8)].copy()

    common = prof_ncpm.index.intersection(cpm.index)
    prof_frac = prof_ncpm.loc[common].div(prof_ncpm.loc[common].sum(axis=0), axis=1)
    bulk_frac = cpm.loc[common].div(cpm.loc[common].sum(axis=0), axis=1)

    # ---------------- instrument 2: specificity-margin markers
    mk2 = [g_ for g_ in margin_markers(prof_ncpm) if g_ in common]
    log(f"instrument 2: specificity-margin markers, shared genes = {len(mk2)} "
        f"(instrument 1 used share-based top-{TOP_MARKERS})")
    A = prof_frac.loc[mk2].to_numpy()
    rows = []
    for _, r in all8.iterrows():
        p, _ = nnls(A, bulk_frac.loc[mk2, r["column"]].to_numpy())
        p = p / p.sum() if p.sum() > 0 else p
        rec = {c: float(v) for c, v in zip(prof_frac.columns, p)}
        rec.update({"column": r["column"], "tissue": r["tissue"], "animal": r["animal"],
                    "day": r["day"], "is_vendor": bool(r["is_vendor"])})
        rec["manual_immune_fraction"] = float(sum(rec[c] for c in a3.IMMUNE_CLASSES))
        rows.append(rec)
    comp2 = pd.DataFrame(rows)
    class_cols = list(prof_frac.columns)
    tc2 = comp2.groupby("tissue")[class_cols].mean()
    tc2["manual_immune_fraction"] = tc2[a3.IMMUNE_CLASSES].sum(axis=1)
    comp2.to_csv(OUT / "composition_inst2.tsv", sep="\t", index=False)
    tc2.to_csv(OUT / "tissue_composition_inst2.tsv", sep="\t")

    log("\n== instrument 2: inferred tissue composition ==")
    show = tc2.copy()
    show["top1"] = show[class_cols].idxmax(axis=1)
    show["top1_frac"] = show[class_cols].max(axis=1)
    log(show[["top1", "top1_frac"]].round(3).to_string())

    a1_rows = []
    for t in tissues8:
        top1 = show.loc[t, "top1"]
        a1_rows.append({"tissue": t, "top1": top1, "top1_frac": round(float(show.loc[t, "top1_frac"]), 3),
                        "expected": ";".join(sorted(a3.ANCHOR_TOP1[t])),
                        "hit": bool(top1 in a3.ANCHOR_TOP1[t])})
    a1 = pd.DataFrame(a1_rows)
    a1_hits = int(a1["hit"].sum())
    log("\n== instrument 2: anchor A1 (same pre-registered rule) ==")
    log(a1.to_string(index=False))
    log(f"instrument 2 anchor A1 hits: {a1_hits}/8")

    # ---------------- marker-based immune content (no deconvolution)
    present = [g_ for g_ in a3.IMMUNE_MARKERS if g_ in cpm.index]
    immune_share = (cpm.loc[present].sum(axis=0) / cpm.sum(axis=0)).rename("immune_cpm_share")
    all8 = all8.assign(immune_cpm_share=all8["column"].map(immune_share))
    tissue_share = all8.groupby("tissue")["immune_cpm_share"].median()
    log("\n== marker-based immune content (share of total CPM), median per tissue ==")
    log(tissue_share.round(4).to_string())
    log(f"markers used: {present}")

    # ---------------- instrument 3: class marker-share weights (no NNLS at all)
    base_markers = [g_ for g_ in a3.marker_genes(prof_ncpm, TOP_MARKERS) if g_ in common]
    marker_by_class: dict[str, list[str]] = {}
    share = prof_ncpm.div(prof_ncpm.sum(axis=1).replace(0, np.nan), axis=0)
    for cls in prof_ncpm.columns:
        ok = prof_ncpm[cls] >= MIN_NCPM
        sel = share.loc[ok, cls].dropna().sort_values(ascending=False).head(TOP_MARKERS).index
        marker_by_class[cls] = [g_ for g_ in sel if g_ in common]
    tot = cpm.loc[base_markers].sum(axis=0)
    w3 = pd.DataFrame(
        {cls: cpm.loc[ms, :].sum(axis=0) for cls, ms in marker_by_class.items()}
    ).div(tot, axis=0)
    w3 = w3.div(w3.sum(axis=1), axis=0)
    comp3 = w3.reset_index().rename(columns={"index": "column"})
    comp3 = all8[["column", "tissue", "animal", "day", "is_vendor"]].merge(comp3, on="column")
    tc3 = comp3.groupby("tissue")[class_cols].mean()
    tc3["manual_immune_fraction"] = tc3[a3.IMMUNE_CLASSES].sum(axis=1)
    comp3.to_csv(OUT / "composition_inst3_markershare.tsv", sep="\t", index=False)
    log("\n== instrument 3: class marker-share (no NNLS) ==")
    show3 = tc3.copy()
    show3["top1"] = show3[class_cols].idxmax(axis=1)
    show3["top1_frac"] = show3[class_cols].max(axis=1)
    log(show3[["top1", "top1_frac"]].round(3).to_string())
    a1b_rows = []
    for t in tissues8:
        top1 = show3.loc[t, "top1"]
        a1b_rows.append({"tissue": t, "top1": top1,
                         "top1_frac": round(float(show3.loc[t, "top1_frac"]), 3),
                         "hit": bool(top1 in a3.ANCHOR_TOP1[t])})
    a1b = pd.DataFrame(a1b_rows)
    a1b_hits = int(a1b["hit"].sum())
    log(a1b.to_string(index=False))
    log(f"instrument 3 anchor A1 hits: {a1b_hits}/8")
    rho_a2d, p_a2d = spearmanr(tc3["manual_immune_fraction"].reindex(tissue_share.index),
                               tissue_share)
    log(f"instrument 3 anchor A2 vs CPM-share index (median): rho = {rho_a2d:.3f} (p = {p_a2d:.4f})")

    rho_a2b, p_a2b = spearmanr(tc2["manual_immune_fraction"].reindex(tissue_share.index),
                               tissue_share)
    rho_a2c, p_a2c = spearmanr(tc2["manual_immune_fraction"].reindex(tissue_share.index),
                               all8.groupby("tissue")["immune_cpm_share"].mean()
                               .reindex(tissue_share.index))
    log(f"\ninstrument 2 anchor A2 vs CPM-share index (median): rho = {rho_a2b:.3f} (p = {p_a2b:.4f})")
    log(f"instrument 2 anchor A2 vs CPM-share index (mean):   rho = {rho_a2c:.3f} (p = {p_a2c:.4f})")

    # in-run comparison for instrument 1 (same index, for a like-for-like read)
    comp1 = pd.read_csv(OUT / "composition.tsv", sep="\t")
    tc1 = comp1.groupby("tissue")[class_cols].mean()
    tc1["manual_immune_fraction"] = tc1[a3.IMMUNE_CLASSES].sum(axis=1)
    rho_ic, p_ic = spearmanr(tc1["manual_immune_fraction"].reindex(tissue_share.index),
                             tissue_share)
    log(f"instrument 1 (pre-registered) vs the same CPM-share index: rho = {rho_ic:.3f} (p = {p_ic:.4f})")

    # ---------------- composition-aware correlation, three instruments
    grad = pd.read_csv(a3.GRAD, sep="\t")
    delta_of_type = dict(zip(grad["cell_type"], grad["delta_restriction"]))

    def delta_weighted(tc: pd.DataFrame, t: str) -> tuple[float, float]:
        w = tc.loc[t]
        cls_delta = {}
        for cls, members in a3.CLASSES.items():
            if cls not in tc.columns:
                continue
            vals = [delta_of_type[m] for m in members if m in delta_of_type]
            if vals:
                cls_delta[cls] = float(np.median(vals))
        cov = float(sum(w[c] for c in cls_delta))
        if cov <= 0:
            return (np.nan, 0.0)
        return (float(sum(w[c] * cls_delta[c] for c in cls_delta) / cov), cov)

    rows = []
    for t in tissues8:
        dw2, cov2 = delta_weighted(tc2, t)
        dw3, cov3 = delta_weighted(tc3, t)
        rows.append({"tissue": t, "delta_weighted_inst2": dw2, "coverage_inst2": cov2,
                     "delta_weighted_inst3": dw3, "coverage_inst3": cov3,
                     "immune_fraction_inst2": float(tc2.loc[t, a3.IMMUNE_CLASSES].sum()),
                     "immune_fraction_inst3": float(tc3.loc[t, a3.IMMUNE_CLASSES].sum()),
                     "immune_cpm_share": float(tissue_share.loc[t])})
    tt = pd.DataFrame(rows)

    # E_adj from the marker-based immune content (no NNLS at all)
    P = all8["immune_cpm_share"].to_numpy(dtype=float)
    s = all8["ISG_core"].to_numpy(dtype=float)
    ok = np.isfinite(P) & np.isfinite(s)
    X = np.column_stack([np.ones(ok.sum()), P[ok]])
    beta, *_ = np.linalg.lstsq(X, s[ok], rcond=None)
    resid = np.full_like(s, np.nan)
    resid[ok] = s[ok] - X @ beta
    all8 = all8.assign(ISG_resid_marker=resid)
    r2 = float(1 - ((s[ok] - X @ beta) ** 2).sum() / ((s[ok] - s[ok].mean()) ** 2).sum())
    log(f"\nmarker-based immune content -> ISG_core: R^2 = {r2:.3f} (n = {int(ok.sum())})")

    inf = all8[all8["day"] >= 1]
    tt["E_adj_marker"] = tt["tissue"].map(inf.groupby("tissue")["ISG_resid_marker"].median())
    tt["E_raw"] = tt["tissue"].map(inf.groupby("tissue")["ISG_core"].median())
    tt["delta_median_g8"] = tt["tissue"].map(
        g8.delta_by_tissue().set_index("tissue")["delta_median"])

    # instrument-2 composition -> E_adj, same recipe as the pre-registered run
    P2 = comp2[class_cols].to_numpy(dtype=float)
    s2 = comp2["column"].map(dict(zip(meta["column"], meta["ISG_core"]))).to_numpy(dtype=float)
    ok2 = np.isfinite(s2) & np.isfinite(P2).all(axis=1)
    X2 = np.column_stack([np.ones(ok2.sum()), P2[ok2]])
    beta2, *_ = np.linalg.lstsq(X2, s2[ok2], rcond=None)
    r2b = float(1 - ((s2[ok2] - X2 @ beta2) ** 2).sum() / ((s2[ok2] - s2[ok2].mean()) ** 2).sum())
    res2 = np.full_like(s2, np.nan)
    res2[ok2] = s2[ok2] - X2 @ beta2
    comp2["ISG_resid"] = res2
    tt["E_adj_inst2"] = tt["tissue"].map(
        comp2[comp2["day"] >= 1].groupby("tissue")["ISG_resid"].median())

    # instrument-3 marker-share composition -> E_adj (identical recipe)
    P3 = comp3[class_cols].to_numpy(dtype=float)
    s3 = comp3["column"].map(dict(zip(meta["column"], meta["ISG_core"]))).to_numpy(dtype=float)
    ok3 = np.isfinite(s3) & np.isfinite(P3).all(axis=1)
    X3 = np.column_stack([np.ones(ok3.sum()), P3[ok3]])
    beta3, *_ = np.linalg.lstsq(X3, s3[ok3], rcond=None)
    r2c = float(1 - ((s3[ok3] - X3 @ beta3) ** 2).sum() / ((s3[ok3] - s3[ok3].mean()) ** 2).sum())
    res3 = np.full_like(s3, np.nan)
    res3[ok3] = s3[ok3] - X3 @ beta3
    comp3["ISG_resid"] = res3
    tt["E_adj_inst3"] = tt["tissue"].map(
        comp3[comp3["day"] >= 1].groupby("tissue")["ISG_resid"].median())
    tt.to_csv(OUT / "a3_sensitivity_table.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    pairs = [
        ("PRIMARY-equivalent (inst2): delta_weighted_inst2 vs E_adj_inst2",
         tt["delta_weighted_inst2"], tt["E_adj_inst2"]),
        ("marker: delta_median_g8 vs E_adj_marker", tt["delta_median_g8"], tt["E_adj_marker"]),
        ("marker: delta_weighted_inst2 vs E_adj_marker", tt["delta_weighted_inst2"], tt["E_adj_marker"]),
        ("inst2 composition: delta_median_g8 vs E_adj_inst2", tt["delta_median_g8"], tt["E_adj_inst2"]),
        ("PRIMARY-equivalent (inst3): delta_weighted_inst3 vs E_adj_inst3",
         tt["delta_weighted_inst3"], tt["E_adj_inst3"]),
        ("inst3 composition: delta_median_g8 vs E_adj_inst3", tt["delta_median_g8"], tt["E_adj_inst3"]),
        ("confound: E_raw vs immune_cpm_share", tt["E_raw"], tt["immune_cpm_share"]),
        ("confound: delta_median_g8 vs immune_cpm_share", tt["delta_median_g8"], tt["immune_cpm_share"]),
    ]
    crows = []
    for label, x, y in pairs:
        rec = a3.corr_block(np.asarray(x, float), np.asarray(y, float), g8, rng)
        rec["pair"] = label
        crows.append(rec)
    corr = pd.DataFrame(crows)[["pair", "n", "rho", "p_perm", "ci_low", "ci_high", "loo_min", "loo_max"]]
    corr.to_csv(OUT / "a3_sensitivity_correlations.tsv", sep="\t", index=False)
    log("\n== sensitivity correlations ==")
    log(corr.to_string(index=False))

    rec = dict(corr.iloc[0])
    log(f"\ninstrument 2 composition -> ISG_core: R^2 = {r2b:.3f}")
    log(pd.DataFrame([rec]).to_string(index=False))

    summary = {
        "role": "exploratory sensitivity; does NOT alter the pre-registered verdict",
        "instrument2_markers": len(mk2),
        "instrument2_anchor_A1_hits": a1_hits,
        "instrument2_anchor_A1_table": a1.to_dict(orient="records"),
        "instrument3_anchor_A1_hits": a1b_hits,
        "instrument3_anchor_A1_table": a1b.to_dict(orient="records"),
        "instrument3_anchor_A2_vs_cpm_share": {"rho": a3.rd(rho_a2d), "p": a3.rd(p_a2d, 4)},
        "instrument3_composition_R2_on_ISG": a3.rd(r2c),
        "instrument2_anchor_A2": {
            "vs_cpm_share_median": {"rho": a3.rd(rho_a2b), "p": a3.rd(p_a2b, 4)},
            "vs_cpm_share_mean": {"rho": a3.rd(rho_a2c), "p": a3.rd(p_a2c, 4)},
        },
        "instrument1_anchor_A2_vs_cpm_share": {"rho": a3.rd(rho_ic), "p": a3.rd(p_ic, 4)},
        "marker_immune_R2_on_ISG": a3.rd(r2),
        "instrument2_composition_R2_on_ISG": a3.rd(r2b),
        "primary_equivalent_inst2": rec,
        "correlations": corr.to_dict(orient="records"),
        "table": tt.round(4).to_dict(orient="records"),
        "markers_used_for_immune_share": present,
        "missing_class_members": missing,
    }
    (OUT / "a3_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    log("\nsensitivity summary written to a3_sensitivity_summary.json")


if __name__ == "__main__":
    main()
