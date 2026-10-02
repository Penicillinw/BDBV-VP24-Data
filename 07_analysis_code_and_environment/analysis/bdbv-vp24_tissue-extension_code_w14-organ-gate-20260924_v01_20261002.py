#!/usr/bin/env python
"""W14: does an entry x post-entry host-factor combination explain organ
tropism, in vivo?

Why this is the test the project was missing
-------------------------------------------
W13 showed that in whole blood nothing at the transcript level separates
infected from bystander cells once sequencing depth is matched, but blood has
no organ structure.  GSE226106 / Normandin 2023 (Cell Genomics, PMID 38169842)
provides exactly the missing resource and it is fully open:

    host RNA-seq   17 tissues x 21 EBOV-infected rhesus macaques (308 samples)
    gene map       GitHub gene_conversion.csv (row index -> Ensembl -> symbol)
    viral load     GitHub qPCR_table.csv (Kulesh EBOV assay per biosample)
    cell mixtures  GitHub MuSiC/TeRNAdecov compositions per sample

Hypothesis under test: the *combination* of an entry arm and a post-entry
programme in a tissue predicts how much virus accumulates in that tissue.
The competing explanation, stated by the authors themselves, is that viral
load simply tracks monocyte content.  Every comparison below is therefore run
with and without the deconvolved monocyte/macrophage fraction.

Outputs (analysis/w14_organ_20260924/)
  sample_table.tsv        per sample: tissue, animal, day, viral load, scores
  tissue_table.tsv        per tissue: median viral load + scores
  tissue_correlations.tsv Spearman of each score vs viral load
  within_animal.tsv       per animal: rank correlation across tissues
  combination_2x2.tsv     entry x machinery median split vs viral load
  regression_summary.txt  OLS models (viral load ~ monocyte +/- layers)
  w14_summary.json
"""

from __future__ import annotations

import io
import json
import re
import urllib.request
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "analysis" / "ext_resource_scan_20260924"
GIT = SCAN / "github"
OUT = ROOT / "analysis" / "w14_organ_20260924"
OUT.mkdir(parents=True, exist_ok=True)
GIT.mkdir(parents=True, exist_ok=True)
SCAN.mkdir(parents=True, exist_ok=True)

RAW = "https://raw.githubusercontent.com/broadinstitute/EbolaNaturalHistory/main/"

DEMO_DECONV = ["liver", "spleen", "lymph_node", "lung", "kidney", "brain",
               "skin", "adrenal"]

MODULES = {
    "AT_lectin": ["CD209", "CLEC10A", "CLEC4G", "CD209L2"],
    "AT_broad": ["HAVCR1", "AXL", "TYRO3", "TIMD4", "MERTK"],
    "EM_core": ["NPC1", "NPC2"],
    "EM_protease": ["CTSB", "CTSL"],
    "EM_traffic": ["TPCN1", "TPCN2", "CD63", "VPS11", "VPS16", "VPS18",
                   "VPS33A", "VPS39", "VPS41", "VPS29", "RAB7A", "PLEKHM1",
                   "GNPTAB", "SLC39A9", "PIK3C3", "TMEM106B"],
    "RE_interferon": ["BST2", "SERINC3", "SERINC5", "MX1", "MX2", "ISG15",
                      "TRIM6"],
    "ISG_control": ["IFIT1", "IFIT2", "IFIT3", "OAS1", "OAS2", "OASL",
                    "RSAD2", "IFI6", "IFI44", "IFI44L"],
    "MONO_marker": ["CD14", "LYZ", "FCN1", "S100A8", "S100A9", "VCAN",
                    "CSF1R", "CD68", "ITGAM", "FCGR3A"],
}
MODULES["EM_all"] = (MODULES["EM_core"] + MODULES["EM_protease"]
                     + MODULES["EM_traffic"])


def fetch_bytes(url: str, timeout: int = 600) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "bdbv/1.0"})
    return urllib.request.urlopen(req, timeout=timeout).read()


def ensure(url: str, local: Path) -> Path:
    if not local.exists() or local.stat().st_size < 100:
        local.write_bytes(fetch_bytes(url))
        print(f"  downloaded {local.name}: {local.stat().st_size:,} bytes")
    return local


def build_sample_table() -> pd.DataFrame:
    counts_path = ensure(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226106/suppl/"
        "GSE226106_20230121_counts_submission.txt.gz",
        SCAN / "GSE226106_counts.txt.gz")
    gene_path = ensure(RAW + "10-TernaDecov/comparisons/music/ebov_with_music/"
                       "01-EBOV_Blood_Deconvolution_Music_with_Seqwel/"
                       "gene_conversion.csv", GIT / "gene_conversion.csv")
    qpcr_path = ensure(RAW + "02-qPCR_Analysis/qPCR_table.csv",
                       GIT / "qPCR_table.csv")

    genes = pd.read_csv(gene_path)
    genes["idx"] = pd.to_numeric(genes.iloc[:, 0], errors="coerce")
    idx_to_symbol = dict(zip(genes["idx"], genes["external_gene_name"]))

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.index = pd.to_numeric(counts.index, errors="coerce")
    counts.columns = [c.strip().strip('"') for c in counts.columns]
    print(f"  counts: {counts.shape[0]} rows x {counts.shape[1]} samples")
    symbols = pd.Series({i: idx_to_symbol.get(i, "") for i in counts.index})
    mapped = int((symbols != "").sum())
    print(f"  gene symbols mapped: {mapped}/{len(symbols)}")

    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    logged = np.log2(cpm + 1)
    logged.index = symbols.values
    logged = logged.groupby(level=0).mean()          # collapse duplicate symbols

    scores = pd.DataFrame(index=logged.columns)
    present: dict[str, list[str]] = {}
    for name, members in MODULES.items():
        found = [g for g in members if g in logged.index]
        present[name] = found
        scores["S_" + name] = logged.loc[found].mean(axis=0) if found else np.nan
        scores["D_" + name] = (logged.loc[found] > 0).mean(axis=0) if found else np.nan
    missing = {k: [g for g in v if g not in logged.index]
               for k, v in MODULES.items()}

    # ---- sample annotation from the series matrix ----
    annot = pd.read_csv(SCAN / "GSE226106_samples.tsv", sep="\t", index_col=0)
    titles = annot.loc["Sample_title"].tolist()
    tissues = annot.loc["Sample_source_name_ch1"].tolist()
    times = annot.loc["Sample_characteristics_ch1"].tolist()
    meta = pd.DataFrame({"column": counts.columns})
    meta["title"] = [t.split(" [")[0] for t in titles]
    meta["biosample"] = [re.sub(r"\[.*\]", "", t).strip() for t in titles]
    meta["tissue"] = [t.split(":")[-1] for t in tissues]
    meta["tissue_short"] = [t.split(":")[0] for t in tissues]
    meta["day"] = [t.replace("time:", "").strip() for t in times]
    meta["animal"] = [m.split("_")[0] for m in meta["biosample"]]
    meta = meta.set_index("column")

    scores = scores.join(meta)

    # ---- viral load ----
    qpcr = pd.read_csv(qpcr_path)
    qpcr.columns = [c.strip() for c in qpcr.columns]
    qpcr["Quant.Kulesh"] = pd.to_numeric(qpcr["Quant.Kulesh"], errors="coerce")
    qpcr["Quant.18s"] = pd.to_numeric(qpcr["Quant.18s"], errors="coerce")
    qpcr["viral_per_18s"] = qpcr["Quant.Kulesh"] / qpcr["Quant.18s"]
    qpcr["log10_viral"] = np.log10(qpcr["viral_per_18s"].clip(lower=1e-12))
    key = qpcr.set_index("Biosample")[["viral_per_18s", "log10_viral",
                                       "Quant.Kulesh", "Tissue"]]
    scores = scores.join(key, on="biosample", how="left")
    print(f"  qPCR merged: {scores['log10_viral'].notna().sum()} / "
          f"{len(scores)} samples")

    # ---- deconvolved monocyte fraction ----
    frac = pd.DataFrame(index=scores.index)
    for tissue in DEMO_DECONV:
        url = (RAW + "10-TernaDecov/deconvolve_with_fascicularis_incl_blood/"
               f"output/{tissue}/{tissue}__sample_quantile_compositions.csv")
        local = GIT / f"{tissue}__sample_quantile_compositions.csv"
        try:
            ensure(url, local)
        except Exception as exc:
            print(f"  deconv {tissue} failed: {exc}")
            continue
        comp = pd.read_csv(local, index_col=0)
        comp = comp[comp["quantile"] == 1]
        pivot = comp.pivot_table(index="sample", columns="celltype",
                                 values="values")
        cell_cols = [c for c in ("Mono", "Kupffer", "DC") if c in pivot.columns]
        myeloid = pivot[cell_cols].sum(axis=1) if cell_cols else pd.Series(dtype=float)
        frac.loc[frac.index.intersection(myeloid.index),
                 "myeloid_fraction"] = myeloid
    scores = scores.join(frac)
    return scores, present, missing


def correlations(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    features = [c for c in table.columns if c.startswith("S_")]
    for feat in features:
        sub = table[[feat, "log10_viral"]].dropna()
        if len(sub) < 5:
            continue
        rho, p = spearmanr(sub[feat], sub["log10_viral"])
        rows.append({"feature": feat, "n": len(sub),
                     "spearman_rho": round(float(rho), 4),
                     "p_value": round(float(p), 5)})
    return pd.DataFrame(rows).sort_values("spearman_rho", ascending=False)


def main() -> None:
    print("building sample table ...")
    samples, present, missing = build_sample_table()
    samples.to_csv(OUT / "sample_table.tsv", sep="\t")

    tissue = (samples.groupby("tissue_short")
              .agg(n=("log10_viral", "size"),
                   median_log10_viral=("log10_viral", "median"),
                   S_AT_lectin=("S_AT_lectin", "median"),
                   S_AT_broad=("S_AT_broad", "median"),
                   S_EM_all=("S_EM_all", "median"),
                   S_RE_interferon=("S_RE_interferon", "median"),
                   S_ISG_control=("S_ISG_control", "median"),
                   S_MONO_marker=("S_MONO_marker", "median"),
                   myeloid_fraction=("myeloid_fraction", "median"))
              .sort_values("median_log10_viral", ascending=False))
    tissue.to_csv(OUT / "tissue_table.tsv", sep="\t")
    print(tissue.to_string())

    corr = correlations(samples)
    corr.to_csv(OUT / "tissue_correlations.tsv", sep="\t", index=False)
    print()
    print(corr.to_string(index=False))

    # within-animal rank correlation across tissues
    rows = []
    for (animal, day), block in samples.groupby(["animal", "day"], observed=True):
        block = block.dropna(subset=["log10_viral"])
        if len(block) < 6:
            continue
        for feat in ["S_AT_lectin", "S_AT_broad", "S_EM_all",
                     "S_RE_interferon", "S_ISG_control", "S_MONO_marker"]:
            rho, p = spearmanr(block[feat], block["log10_viral"])
            rows.append({"animal": animal, "day": day, "n_tissues": len(block),
                         "feature": feat, "rho": round(float(rho), 4),
                         "p": round(float(p), 4)})
    within = pd.DataFrame(rows)
    within.to_csv(OUT / "within_animal.tsv", sep="\t", index=False)

    # 2x2 combination classes (median split over samples)
    work = samples.dropna(subset=["log10_viral", "S_AT_lectin", "S_EM_all"]).copy()
    work["entry_high"] = work["S_AT_lectin"] > work["S_AT_lectin"].median()
    work["machine_high"] = work["S_EM_all"] > work["S_EM_all"].median()
    combo = (work.groupby(["entry_high", "machine_high"], observed=True)
             .agg(n=("log10_viral", "size"),
                  median_log10_viral=("log10_viral", "median"),
                  mean_log10_viral=("log10_viral", "mean"),
                  median_myeloid=("myeloid_fraction", "median"))
             .reset_index())
    combo.to_csv(OUT / "combination_2x2.tsv", sep="\t", index=False)

    # regressions with and without the myeloid covariate
    lines = []
    reg = samples.dropna(subset=["log10_viral", "S_AT_lectin", "S_EM_all",
                                 "S_RE_interferon", "S_MONO_marker"]).copy()
    reg["entry"] = (reg["S_AT_lectin"] - reg["S_AT_lectin"].mean()) / reg["S_AT_lectin"].std()
    reg["machine"] = (reg["S_EM_all"] - reg["S_EM_all"].mean()) / reg["S_EM_all"].std()
    reg["restrict"] = (reg["S_RE_interferon"] - reg["S_RE_interferon"].mean()) / reg["S_RE_interferon"].std()
    reg["mono"] = (reg["S_MONO_marker"] - reg["S_MONO_marker"].mean()) / reg["S_MONO_marker"].std()
    reg["interaction"] = reg["entry"] * reg["machine"]
    models = {
        "m0_mono": ["mono"],
        "m1_entry": ["entry"],
        "m2_machine": ["machine"],
        "m3_restrict": ["restrict"],
        "m4_entry_machine": ["entry", "machine"],
        "m5_combo_int": ["entry", "machine", "interaction"],
        "m6_mono_entry_machine": ["mono", "entry", "machine"],
        "m7_mono_all_int": ["mono", "entry", "machine", "restrict", "interaction"],
    }
    y = reg["log10_viral"]
    for name, feats in models.items():
        x = sm.add_constant(reg[feats])
        fit = sm.OLS(y, x).fit()
        lines.append(f"{name}: n={int(fit.nobs)} R2={fit.rsquared:.3f} "
                     f"adjR2={fit.rsquared_adj:.3f} AIC={fit.aic:.1f}")
        for term in feats:
            lines.append(f"    {term:12s} beta={fit.params[term]: .3f} "
                         f"p={fit.pvalues[term]:.4f}")
    (OUT / "regression_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print()
    print("\n".join(lines))

    summary = {
        "date": "2026-09-24",
        "resource": {
            "geo": "GSE226106",
            "paper": "Normandin et al. 2023 Cell Genomics, PMID 38169842, "
                     "doi 10.1016/j.xgen.2023.100440",
            "n_samples": int(len(samples)),
            "n_tissues": int(samples["tissue_short"].nunique()),
            "n_animals": int(samples["animal"].nunique()),
            "gene_rows": 35405,
            "gene_symbols_mapped": "see log",
            "viral_load_source": "qPCR_table.csv Quant.Kulesh / Quant.18s",
            "composition_source": "TeRNAdecov/MuSiC sample_quantile_compositions",
        },
        "panels_present": present,
        "panels_missing": missing,
        "tissue_correlations": corr.to_dict(orient="records"),
        "combination_2x2": combo.to_dict(orient="records"),
        "limitations": [
            "bulk tissue: scores reflect tissue composition as much as per-cell "
            "state, hence the monocyte/myeloid covariate",
            "qPCR viral load is from a different assay than RNA-seq; both come "
            "from the same biosamples but are normalised differently",
            "tissue-specific panels are not adjusted for RNA composition "
            "(e.g. haemoglobin in spleen) beyond CPM",
            "deconvolution is only available for 8 of the tissues",
        ],
    }
    (OUT / "w14_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8")
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
