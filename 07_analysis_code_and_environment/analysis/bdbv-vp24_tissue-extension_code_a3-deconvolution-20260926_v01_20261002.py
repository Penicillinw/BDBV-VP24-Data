#!/usr/bin/env python
"""A3: is G8's tissue-level conflict a *tissue composition* artefact?

Pre-registration: report/A3_预注册_组织构成去卷积_20260926.md
(written before any new number in this script was computed).

Question
--------
G8 compared the tissue-median Delta against the bulk-tissue ISG level E_raw and got
rho = -0.333 (n = 8 tissues, permutation p = 0.43), with 4/8 tissues conflicting.
Both sides can be composition-confounded:
  * E_raw = whole-tissue ISG level -> immune-rich tissues score high regardless of
    the IFN response of their parenchyma;
  * Delta_median = median over *named* cell types -> ignores how abundant those
    cells actually are in the tissue.
This script estimates the composition of every macaque tissue sample by NNLS against
a human single-cell reference (HPA `rna_single_cell_type.tsv`), then rebuilds both sides
composition-aware:
  * Delta_weighted(t): abundance-weighted class Delta, coverage reported;
  * E_adj(t): residual of the ISG_core score after regressing out the inferred
    composition (fit on all samples of the 8 tissues, residual taken per sample,
    tissue median over infected samples).

All statistics reuse the frozen G8 implementations (imported from
tools/t4_g8_tissue_ifn_env_20260926.py) so that the reproduction check of
rho(Delta_median, E_raw) = -0.333 is a genuine pipeline regression test.

Hard boundaries (see pre-registration §4): human reference vs macaque tissue vs
EBOV Makona (not BDBV); relative proportions only; NNLS with `len(CLASSES)` classes is
under-determined; this cannot test the multiplier and must not alter the frozen
Delta values or ranking.

Class-count deviation (registered 2026-09-26)
---------------------------------------------
Pre-registration §2.1 fixed the class table at 40 classes, but the executed `CLASSES`
constant is 42 (it adds `megakaryocytes_platelets` and `erythrocytes`). The deviation is
recorded in `report/A3_修订记录1_细胞类数_20260926.md`, together with a re-run of the
whole pipeline at exactly the pre-registered 40 classes
(`tools/a3_classcount40_sensitivity_20260926.py`): no gate, no verdict and no sign is
changed, but the A2 anchor moves from 0.405 to 0.262 (both below the 0.5 pass line).

Outputs (analysis/a3_deconvolution_20260926/)
    composition.tsv, tissue_composition.tsv, a3_tissue_table.tsv,
    a3_correlations.tsv, anchors.json, a3_summary.json, a3_log.txt
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
OUT.mkdir(parents=True, exist_ok=True)
G8_TOOL = ROOT / "tools" / "t4_g8_tissue_ifn_env_20260926.py"
HPA_SC = ROOT / "data" / "hpa" / "rna_single_cell_type.tsv"
GRAD = ROOT / "analysis" / "t4_ifn_landscape_20260925" / "restriction_gradient.tsv"

SEED = 20260926
N_PERM = 5000
N_BOOT = 5000
TOP_MARKERS = 200
MARKER_SENSITIVITY = (100, 200, 400)

# ---------------------------------------------------------------- cell classes
# Frozen in the pre-registration BEFORE any composition was estimated.
CLASSES: dict[str, list[str]] = {
    "hepatocytes": ["hepatocytes"],
    "cholangiocytes": ["cholangiocytes"],
    "hepatic_stellate": ["hepatic stellate cells"],
    "macrophages": ["kupffer cells", "macrophages", "hofbauer cells"],
    "kidney_tubule": [
        "proximal tubule cells", "distal convoluted tubule cells",
        "loop of henle epithelial cells", "renal collecting duct intercalated cells",
        "renal collecting duct principal cells", "renal connecting tubule cells",
        "papillary tip epithelial cells",
    ],
    "podocytes": ["podocytes"],
    "neurons": ["brain excitatory neurons", "brain inhibitory neurons", "other brain neurons"],
    "astrocytes": ["astrocytes", "bergmann glia"],
    "oligodendrocytes": ["oligodendrocytes", "oligodendrocyte progenitor cells"],
    "microglia": ["microglia"],
    "ependymal": ["ependymal cells", "choroid plexus epithelial cells"],
    "adrenal_cortex": ["adrenal cortex cells"],
    "adrenal_medulla": ["adrenal medulla cells"],
    "t_cells": ["t-cells"],
    "b_cells": ["b-cells"],
    "nk_cells": ["nk-cells"],
    "plasma_cells": ["plasma cells"],
    "innate_lymphoid": ["innate lymphoid cells"],
    "monocytes": ["monocytes", "monocyte progenitors"],
    "neutrophils": ["neutrophils", "neutrophil progenitors"],
    "mast_cells": ["mast cells"],
    "dendritic_cells": ["cdc", "pdcs"],
    "lymphatic_endothelial": ["lymphatic endothelial cells"],
    "vascular_endothelial": ["vascular endothelial cells"],
    "fibroblasts": ["fibroblasts", "fibro-adipogenic progenitors"],
    "smooth_muscle": ["smooth muscle cells", "vascular smooth muscle cells", "pericytes"],
    "keratinocytes": ["basal keratinocytes", "suprabasal keratinocytes"],
    "melanocytes": ["melanocytes"],
    "alveolar_epithelium": [
        "alveolar cells type 1", "alveolar cells type 2", "transitional alveolar cells",
    ],
    "airway_epithelium": [
        "respiratory basal cells", "respiratory ciliated cells",
        "respiratory deuterosomal cells", "respiratory ionocytes",
        "respiratory secretory cells",
    ],
    "germ_cells": [
        "undifferentiated spermatogonia", "differentiating spermatogonia",
        "early primary spermatocytes", "late primary spermatocytes",
        "early spermatids", "late spermatids",
    ],
    "sertoli": ["sertoli cells"],
    "leydig": ["leydig cells"],
    "peritubular_myoid": ["peritubular myoid cells"],
    "granulosa": ["granulosa cells"],
    "oocytes": ["oocytes"],
    "ovarian_stroma": ["ovarian stromal cells"],
    "epididymal_epithelium": [
        "epididymal basal cells", "epididymal clear cells", "epididymal principal cells",
        "epididymal efferent duct absorptive cells", "epididymal efferent duct ciliated cells",
    ],
    "adipocytes": ["adipocytes"],
    "hematopoietic_progenitors": [
        "hematopoietic stem cells", "megakaryocyte progenitors",
        "megakaryocyte-erythroid progenitors",
    ],
    "megakaryocytes_platelets": ["megakaryocytes", "platelets"],
    "erythrocytes": ["erythrocytes", "erythrocyte progenitors"],
}

# Immune classes, used only for the A2 anchor and the confound-existence check.
IMMUNE_CLASSES = [
    "macrophages", "t_cells", "b_cells", "nk_cells", "plasma_cells",
    "innate_lymphoid", "monocytes", "neutrophils", "mast_cells",
    "dendritic_cells", "hematopoietic_progenitors",
]
IMMUNE_MARKERS = ["PTPRC", "CD3E", "MS4A1", "CD20", "NKG7", "LYZ", "CD68", "CD14"]

ANCHOR_TOP1 = {
    "liver": {"hepatocytes"},
    "kidney": {"kidney_tubule"},
    "brain": {"neurons", "astrocytes", "oligodendrocytes", "microglia", "ependymal"},
    "adrenal gland": {"adrenal_cortex", "adrenal_medulla"},
    "lymph node": {"t_cells", "b_cells", "nk_cells", "plasma_cells", "innate_lymphoid"},
    "skin": {"keratinocytes", "fibroblasts"},
    "sex organ": {"germ_cells", "sertoli", "leydig", "epididymal_epithelium"},
    "lung": {"alveolar_epithelium", "airway_epithelium"},
}
ANCHOR_MIN_HITS = 6
ANCHOR_MIN_RHO_IMMUNE = 0.5


def log(msg: str) -> None:
    print(msg)
    with open(OUT / "a3_log.txt", "a", encoding="utf-8") as fh:
        fh.write(msg + "\n")


# ------------------------------------------------------------- g8 pipeline reuse

def load_g8():
    spec = importlib.util.spec_from_file_location("g8_tool", G8_TOOL)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    # route the imported module's logging into our own log file
    mod.log = log
    return mod


# ------------------------------------------------------------------ reference

def load_hpa_reference() -> tuple[pd.DataFrame, pd.DataFrame]:
    """gene symbol x cell type, linear nCPM (mean over duplicate symbols)."""
    hpa = pd.read_csv(HPA_SC, sep="\t", usecols=["Gene name", "Cell type", "nCPM"])
    hpa = hpa.dropna(subset=["Gene name", "Cell type"])
    # HPA leaves a cell type blank for a gene that was not detected there.
    # Pre-registration §2.2 fixes the convention: not detected = 0 nCPM.
    hpa["nCPM"] = pd.to_numeric(hpa["nCPM"], errors="coerce")
    n_blank = int(hpa["nCPM"].isna().sum())
    wide = hpa.pivot_table(
        index="Gene name", columns="Cell type", values="nCPM",
        aggfunc="mean", observed=True,
    )
    wide = wide.fillna(0.0)
    log(f"HPA reference: {n_blank} (gene, cell type) pairs are blank -> set to 0 nCPM "
        f"(pre-registration convention)")
    return wide, hpa


def resolve_member(wide: pd.DataFrame, name: str) -> str | None:
    if name in wide.columns:
        return name
    # tolerant match for the mueller-glia label (encoding mangling)
    if name == "müller glia":
        hits = [c for c in wide.columns if isinstance(c, str) and c.endswith("ller glia")]
        return hits[0] if hits else None
    return None


def class_profiles(wide: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    cols, missing = {}, []
    for cls, members in CLASSES.items():
        resolved = []
        for m in members:
            r = resolve_member(wide, m)
            if r is None:
                missing.append(f"{cls}:{m}")
            else:
                resolved.append(r)
        if resolved:
            cols[cls] = wide[resolved].mean(axis=1)
    return pd.DataFrame(cols), missing


def marker_genes(prof: pd.DataFrame, n_top: int, min_ncpm: float = 10.0) -> list[str]:
    share = prof.div(prof.sum(axis=1).replace(0, np.nan), axis=0)
    keep: set[str] = set()
    for cls in prof.columns:
        ok = prof[cls] >= min_ncpm
        s = share.loc[ok, cls].dropna().sort_values(ascending=False)
        keep.update(s.head(n_top).index)
    return sorted(keep)


# ------------------------------------------------------------------ statistics

def rd(x, nd=3):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), nd)


def corr_block(x: np.ndarray, y: np.ndarray, g8, rng) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 4:
        return {"n": int(mask.sum()), "rho": None, "p_perm": None,
                "ci_low": None, "ci_high": None, "loo_min": None, "loo_max": None}
    x, y = x[mask], y[mask]
    ci = g8.boot_ci(x, y, rng)
    return {
        "n": int(mask.sum()),
        "rho": rd(spearmanr(x, y).statistic),
        "p_perm": rd(g8.perm_p(x, y, rng), 4),
        "ci_low": rd(ci[0]), "ci_high": rd(ci[1]),
        **{k: rd(v) for k, v in g8.loo(x, y).items()},
    }


def main() -> None:
    (OUT / "a3_log.txt").write_text("", encoding="utf-8")
    g8 = load_g8()

    # ---- reference
    wide, hpa_long = load_hpa_reference()
    prof_ncpm, missing_members = class_profiles(wide)
    log(f"HPA single-cell reference: {wide.shape[0]} genes x {wide.shape[1]} cell types")
    log(f"pre-registered classes with >=1 resolved member: {prof_ncpm.shape[1]} / {len(CLASSES)}")
    log(f"members not found in HPA (reported, not substituted): {missing_members}")

    # ---- bulk (identical gene mapping / filtering to G8)
    meta = g8.load_sample_meta()
    z, depth, logcpm = g8.load_expression()
    meta = meta[meta["column"].isin(z.columns)].copy()
    meta["ISG_core"] = [g8.panel_score(z, g8.ISG_CORE, c) for c in meta["column"]]
    meta["genes_detected"] = meta["column"].map(depth)
    cpm = (2.0 ** logcpm) - 1.0            # logcpm = log2(cpm+1)

    tissues8 = list(ANCHOR_TOP1)
    all8 = meta[meta["tissue"].isin(tissues8)].copy()
    log(f"samples in the 8 primary tissues (all conditions): {len(all8)}; "
        f"infected: {(all8['day'] >= 1).sum()}; controls/vendor: {(all8['day'] <= 0).sum()}")

    # ---- composition by NNLS, all samples of the 8 tissues
    common = prof_ncpm.index.intersection(cpm.index)
    log(f"genes shared between HPA reference and macaque matrix: {len(common)}")
    prof_frac = prof_ncpm.loc[common].div(prof_ncpm.loc[common].sum(axis=0), axis=1)
    bulk_frac = cpm.loc[common].div(cpm.loc[common].sum(axis=0), axis=1)

    marker_sets = {n: marker_genes(prof_ncpm, n) for n in MARKER_SENSITIVITY}
    markers = [g for g in marker_sets[TOP_MARKERS] if g in common]
    log(f"marker genes (top-{TOP_MARKERS} per class, shared): {len(markers)}")

    comp_rows = []
    for _, r in all8.iterrows():
        col = r["column"]
        A = prof_frac.loc[markers].to_numpy()
        b = bulk_frac.loc[markers, col].to_numpy()
        p, _ = nnls(A, b)
        tot = p.sum()
        p = p / tot if tot > 0 else p
        rec = {cls: float(v) for cls, v in zip(prof_frac.columns, p)}
        rec.update({"column": col, "tissue": r["tissue"], "animal": r["animal"],
                    "day": r["day"], "is_vendor": bool(r["is_vendor"]),
                    "manual_immune_fraction": float(sum(rec[c] for c in IMMUNE_CLASSES))})
        comp_rows.append(rec)
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "composition.tsv", sep="\t", index=False)

    class_cols = list(prof_frac.columns)
    tissue_comp = comp.groupby("tissue")[class_cols].mean()
    tissue_comp.to_csv(OUT / "tissue_composition.tsv", sep="\t")

    log("\n== inferred tissue composition (mean over all samples of the tissue) ==")
    show = tissue_comp.copy()
    show["top1"] = show[class_cols].idxmax(axis=1)
    show["top1_frac"] = show[class_cols].max(axis=1)
    log(show[["top1", "top1_frac"] + IMMUNE_CLASSES].round(3).to_string())

    # ---- anchor A1: top-1 class
    anchor_rows = []
    for t in tissues8:
        top1 = show.loc[t, "top1"]
        hit = top1 in ANCHOR_TOP1[t]
        exp = show.loc[t, class_cols].drop(top1, errors="ignore").idxmax()
        anchor_rows.append({"tissue": t, "top1": top1, "top1_frac": rd(show.loc[t, "top1_frac"]),
                            "expected": ";".join(sorted(ANCHOR_TOP1[t])), "hit": bool(hit),
                            "runner_up": exp, "runner_up_frac": rd(show.loc[t, exp])})
    a1 = pd.DataFrame(anchor_rows)
    a1_hits = int(a1["hit"].sum())
    log("\n== anchor A1 (pre-registered top-1 class) ==")
    log(a1.to_string(index=False))
    log(f"anchor A1 hits: {a1_hits}/{len(tissues8)} (pass >= {ANCHOR_MIN_HITS})")

    # ---- anchor A2: inferred immune fraction vs bulk immune-marker index
    present_markers = [g for g in IMMUNE_MARKERS if g in z.index]
    log(f"immune markers present in the macaque matrix: {present_markers} "
        f"(absent: {[g for g in IMMUNE_MARKERS if g not in z.index]})")
    zmeta = meta.set_index("column")
    immune_index = {}
    for _, r in all8.iterrows():
        col = r["column"]
        immune_index[col] = float(z.loc[present_markers, col].mean())
    comp["immune_marker_index"] = comp["column"].map(immune_index)

    tissue_immune = comp.groupby("tissue")[["manual_immune_fraction", "immune_marker_index"]].mean()
    rho_a2, p_a2 = spearmanr(tissue_immune["manual_immune_fraction"], tissue_immune["immune_marker_index"])
    log("\n== anchor A2: inferred immune fraction vs bulk immune-marker index (8 tissues) ==")
    log(tissue_immune.round(3).to_string())
    log(f"anchor A2: rho = {rho_a2:.3f} (perm-free Spearman p = {p_a2:.4f}); "
        f"pass >= {ANCHOR_MIN_RHO_IMMUNE}")

    a3_depth = {}
    for t in tissues8:
        d = comp[comp["tissue"] == t]
        if d.shape[0] >= 3 and d["manual_immune_fraction"].nunique() >= 2:
            rho, p = spearmanr(d["manual_immune_fraction"], d["column"].map(depth))
            a3_depth[t] = {"rho_vs_depth": rd(rho), "p": rd(p, 4), "n": int(d.shape[0])}
    log("\n== anchor A3 (not a pass/fail gate): inferred immune fraction vs genes_detected ==")
    log(json.dumps(a3_depth, ensure_ascii=False))

    gates = {
        "A1_top1_hits": a1_hits,
        "A1_required": ANCHOR_MIN_HITS,
        "A1_pass": bool(a1_hits >= ANCHOR_MIN_HITS),
        "A2_rho_immune": rd(rho_a2),
        "A2_required": ANCHOR_MIN_RHO_IMMUNE,
        "A2_pass": bool(np.isfinite(rho_a2) and rho_a2 >= ANCHOR_MIN_RHO_IMMUNE),
    }
    gates["all_pass"] = bool(gates["A1_pass"] and gates["A2_pass"])

    # ---- Delta_weighted
    grad = pd.read_csv(GRAD, sep="\t")
    delta_of_type = dict(zip(grad["cell_type"], grad["delta_restriction"]))
    class_delta, class_n = {}, {}
    for cls, members in CLASSES.items():
        if cls not in tissue_comp.columns:
            continue
        vals = [delta_of_type[m] for m in members if m in delta_of_type]
        if vals:
            class_delta[cls] = float(np.median(vals))
            class_n[cls] = len(vals)
    log(f"\nclasses with a frozen Delta value: {len(class_delta)} / {tissue_comp.shape[1]}")

    rows = []
    for t in tissues8:
        w = tissue_comp.loc[t]
        cov = float(sum(w[c] for c in class_delta))
        dw = float(sum(w[c] * class_delta[c] for c in class_delta) / cov) if cov > 0 else np.nan
        rows.append({"tissue": t, "delta_weighted": dw, "coverage": rd(cov),
                     "delta_median_g8": float(g8.delta_by_tissue()
                                              .set_index("tissue").loc[t, "delta_median"]),
                     "immune_fraction": float(w[IMMUNE_CLASSES].sum()),
                     "depth_median_infected": float(
                         all8[(all8["tissue"] == t) & (all8["day"] >= 1)]["genes_detected"].median())})
    tt = pd.DataFrame(rows)

    # ---- E_adj: residualise the sample-level ISG score on the inferred composition
    P = comp[class_cols].to_numpy()
    s = comp["column"].map(dict(zip(meta["column"], meta["ISG_core"]))).to_numpy()
    ok = np.isfinite(s) & np.isfinite(P).all(axis=1)
    X = np.column_stack([np.ones(ok.sum()), P[ok]])
    beta, *_ = np.linalg.lstsq(X, s[ok], rcond=None)
    fit = X @ beta
    resid = np.full_like(s, np.nan)
    resid[ok] = s[ok] - fit
    ss_tot = float(((s[ok] - s[ok].mean()) ** 2).sum())
    r2 = float(1.0 - ((s[ok] - fit) ** 2).sum() / ss_tot) if ss_tot > 0 else np.nan
    comp["ISG_core"] = s
    comp["ISG_resid"] = resid
    log(f"\ncomposition -> ISG_core regression: R^2 = {r2:.3f} "
        f"(fit on {int(ok.sum())} samples of the 8 tissues, all conditions)")

    inf = comp[comp["day"] >= 1]
    e_raw = inf.groupby("tissue")["ISG_core"].median()
    e_adj = inf.groupby("tissue")["ISG_resid"].median()
    tt["E_raw"] = tt["tissue"].map(e_raw)
    tt["E_adj"] = tt["tissue"].map(e_adj)
    tt["n_infected"] = tt["tissue"].map(inf.groupby("tissue").size())
    tt.to_csv(OUT / "a3_tissue_table.tsv", sep="\t", index=False)

    log("\n== tissue table ==")
    log(tt.round(3).to_string(index=False))

    # ---- correlations
    rng = np.random.default_rng(SEED)
    pairs = [
        ("PRIMARY: delta_weighted vs E_adj", tt["delta_weighted"], tt["E_adj"], "primary"),
        ("reproduction: delta_median_g8 vs E_raw", tt["delta_median_g8"], tt["E_raw"], "check"),
        ("descriptive: delta_weighted vs E_raw", tt["delta_weighted"], tt["E_raw"], "descriptive"),
        ("descriptive: delta_median_g8 vs E_adj", tt["delta_median_g8"], tt["E_adj"], "descriptive"),
        ("confound check: E_raw vs immune_fraction", tt["E_raw"], tt["immune_fraction"], "confound"),
        ("confound check: delta_median vs immune_fraction", tt["delta_median_g8"],
         tt["immune_fraction"], "confound"),
        ("confound check: delta_weighted vs immune_fraction", tt["delta_weighted"],
         tt["immune_fraction"], "confound"),
    ]
    crows = []
    for label, x, y, kind in pairs:
        rec = corr_block(np.asarray(x, dtype=float), np.asarray(y, dtype=float), g8, rng)
        rec.update({"pair": label, "kind": kind})
        # depth-residualised version (both axes, as in G8)
        dv = tt["depth_median_infected"].to_numpy(dtype=float)
        xr = g8.residualise(np.asarray(x, dtype=float), dv)
        yr = g8.residualise(np.asarray(y, dtype=float), dv)
        m = np.isfinite(xr) & np.isfinite(yr)
        rec["rho_depth_residual"] = rd(spearmanr(xr[m], yr[m]).statistic) if m.sum() >= 5 else None
        crows.append(rec)
    corr = pd.DataFrame(crows)[
        ["kind", "pair", "n", "rho", "p_perm", "ci_low", "ci_high",
         "loo_min", "loo_max", "rho_depth_residual"]
    ]
    corr.to_csv(OUT / "a3_correlations.tsv", sep="\t", index=False)
    log("\n== correlations ==")
    log(corr.to_string(index=False))

    # ---- marker-set sensitivity for the composition itself
    sens = {}
    for n in MARKER_SENSITIVITY:
        if n == TOP_MARKERS:
            sens[str(n)] = {c: 1.0 for c in class_cols}
            continue
        mk = [g for g in marker_sets[n] if g in common]
        Pn = []
        for col in comp["column"]:
            p, _ = nnls(prof_frac.loc[mk].to_numpy(), bulk_frac.loc[mk, col].to_numpy())
            Pn.append(p / p.sum() if p.sum() > 0 else p)
        Pn = pd.DataFrame(Pn, columns=prof_frac.columns)
        vals = {}
        for c in class_cols:
            rho, _ = spearmanr(tissue_comp[c].to_numpy(),
                               Pn.assign(tissue=comp["tissue"].to_numpy())
                                 .groupby("tissue")[c].mean()
                                 .reindex(tissue_comp.index).to_numpy())
            vals[c] = rd(rho)
        sens[str(n)] = vals
    log("\n== marker-set sensitivity (per-class Spearman of tissue composition, top-N vs top-200) ==")
    sens_df = pd.DataFrame(sens)
    log(sens_df.round(3).to_string())

    # ---- verdict (pre-registered)
    primary = corr[corr["kind"] == "primary"].iloc[0]
    repro = corr[corr["kind"] == "check"].iloc[0]
    if not gates["all_pass"]:
        verdict = "不可判定：去卷积未通过锚点（A1/A2 至少一条未过）"
    else:
        rho, p = primary["rho"], primary["p_perm"]
        resid_rho = primary["rho_depth_residual"]
        if rho is not None and rho >= 0.5 and p is not None and p < 0.05 and (
            resid_rho is None or np.sign(resid_rho) == np.sign(rho)
        ):
            verdict = "构成混杂可解释 G8 的失败：构成校正后 Δ 与组织 IFN 环境转为一致"
        elif rho is not None and rho <= 0 and (p is None or p > 0.05):
            verdict = "不是构成问题：构成校正后仍不一致，组织级外推应继续撤回"
        else:
            verdict = "混合/不可判定：只能报告两个量各自方向，不得给组织级结论"

    repro_ok = (repro["rho"] is not None and abs(repro["rho"] - (-0.333)) < 0.02
                and int(repro["n"]) == 8)
    if not repro_ok:
        verdict = "管线一致性检查未通过（G8 复现量偏离 −0.333）→ 本件整体作废"

    summary = {
        "preregistration": "report/A3_预注册_组织构成去卷积_20260926.md",
        "seed": SEED,
        "reference": "data/hpa/rna_single_cell_type.tsv (HPA single-cell nCPM, human, uninfected)",
        "bulk": "GSE226106 (Normandin 2023, PMID 38169842) macaque EBOV Makona, 8 tissues",
        "n_classes": int(tissue_comp.shape[1]),
        "n_marker_genes": len(markers),
        "missing_class_members": missing_members,
        "anchors": gates,
        "anchor_A1_table": a1.to_dict(orient="records"),
        "anchor_A2_table": tissue_immune.round(4).reset_index().to_dict(orient="records"),
        "anchor_A3_not_a_gate": a3_depth,
        "composition_R2_on_ISG": rd(r2),
        "marker_sensitivity": sens,
        "pipeline_reproduction": {"pair": "delta_median_g8 vs E_raw", **repro.to_dict()},
        "reproduction_ok": bool(repro_ok),
        "tissue_table": tt.round(4).to_dict(orient="records"),
        "correlations": corr.to_dict(orient="records"),
        "verdict": verdict,
        "limitations": [
            "human HPA reference vs macaque tissue vs EBOV Makona (not BDBV): cannot test the multiplier",
            f"relative proportions only; NNLS with {len(CLASSES)} classes is under-determined "
            f"({len(CLASSES)} executed vs 40 pre-registered; see report/A3_修订记录1_细胞类数_20260926.md)",
            "frozen Delta values and ranking are not recomputed here (decision D-011)",
        ],
    }
    (OUT / "anchors.json").write_text(
        json.dumps({"anchors": gates, "A1": a1.to_dict(orient="records"),
                    "A2": tissue_immune.round(4).reset_index().to_dict(orient="records"),
                    "A3": a3_depth}, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8")
    (OUT / "a3_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8")

    log(f"\nVERDICT: {verdict}")
    log(f"primary rho = {primary['rho']} (p_perm = {primary['p_perm']}, "
        f"depth-residual rho = {primary['rho_depth_residual']})")


if __name__ == "__main__":
    main()
