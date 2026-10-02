#!/usr/bin/env python
"""T4-G8: tissue IFN environment during an in-vivo ebolavirus infection vs the
frozen cell-type "capacity" score Delta.

Pre-registration: report/T4_G8_预注册_组织IFN环境与Δ对照_20260926.md
(written before any number in this script was computed).

Question
--------
Delta is built from *uninfected* human cell-type expression, i.e. it measures
capacity. This script asks whether the tissue-level ordering of Delta agrees
with the IFN environment of the same tissues during a real in-vivo ebolavirus
infection, and which tissues conflict.

Resource (local): GSE226106 / Normandin et al. 2023 Cell Genomics
(PMID 38169842) -- 17 tissues x 21 EBOV-infected rhesus macaques, 308 samples
Plus the study's own uninfected controls (RA1082/RA1819/RA1856, day 0) and the
Zyagen vendor uninfected tissues shipped with the same deposition.

Limits that must be carried into any downstream sentence
--------------------------------------------------------
* macaque, not human; EBOV (Makona), not BDBV -> this cannot test the
  "multiplier" itself, only whether Delta is a usable tissue-level proxy for
  the IFN environment.
* bulk tissue averages, not cell-type resolution (W14 limitation).

Outputs (analysis/t4_g8_tissue_ifn_env_20260926/)
    sample_scores.tsv, tissue_env.tsv, tissue_delta.tsv,
    g8_correlations.tsv, g8_concordance.tsv, g8_bootstrap.tsv,
    g8_depth_control.tsv, g8_summary.json, g8_log.txt
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "analysis" / "ext_resource_scan_20260924"
GIT = SCAN / "github"
WIDE = ROOT / "analysis" / "t4_ifn_landscape_20260925" / "hpa_ifn_landscape_wide.tsv"
GRAD = ROOT / "analysis" / "t4_ifn_landscape_20260925" / "restriction_gradient.tsv"
OUT = ROOT / "analysis" / "t4_g8_tissue_ifn_env_20260926"
OUT.mkdir(parents=True, exist_ok=True)
# True row names of the GSE226106 count matrix, recovered from the study's own
# 01-data/se.qc.rds (rdata parser). The submission file itself has integer
# row indices only, and github/gene_conversion.csv is a biomaRt output whose
# first column is a write.csv row index -- NOT a row id. See report for the
# audit that established this.
ROWNAMES = OUT / "ref" / "se_qc_rownames.txt"

SEED = 20260926
N_PERM = 5000
N_BOOT = 5000
LATE_DAY = 4

# ISG_core must be identical to Delta's ISG priming component (frozen in
# tools/t4_restriction_gradient_20260925.py).
ISG_CORE = ["ISG15", "MX1"]
# ISG_broad reuses the RESTRICT/ISG panel frozen in tools/w17_innate_baseline_20260924.py
ISG_BROAD = [
    "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT2", "IFIT3",
    "IFI6", "IFI27", "IFI44", "BST2", "IFITM1", "IFITM2", "IFITM3", "GBP1",
    "GBP2", "GBP5", "SAMHD1", "TRIM5", "CH25H", "RSAD2", "XAF1", "USP18",
    "STAT2", "IRF7", "DDX58", "IFIH1", "EIF2AK2", "RNASEL", "PLSCR1", "SP100",
    "PML", "LY6E",
]
IFN_REC = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB"]

# ---------------------------------------------------------------- sample meta

TISSUE_PATTERNS = [
    (r"^Liver", "liver"),
    (r"^Spleen", "spleen"),
    (r"^Kidney", "kidney"),
    (r"^Brain", "brain"),
    (r"^Adrenal", "adrenal gland"),
    (r"^LN(\b|[-_])", "lymph node"),
    (r"^Skin", "skin"),
    (r"^Sex[-_]?Organ", "sex organ"),
    (r"^Ovary", "ovary"),
    (r"^Testis", "testis"),
    (r"^Lung", "lung"),
    (r"^SpinalCord", "spinal cord"),
    (r"^WholeBlood", "whole blood"),
    (r"^PBMC", "PBMC"),
]

# Explicit cell-type -> tissue map (frozen before any correlation was computed).
# Only cell types that have an anatomically matching tissue in the macaque
# panel are mapped; nothing is mapped by a posteriori convenience.
CELLTYPE_TISSUE = {
    "liver": ["hepatocytes", "kupffer cells", "cholangiocytes", "hepatic stellate cells"],
    "kidney": [
        "proximal tubule cells", "distal convoluted tubule cells",
        "loop of henle epithelial cells", "podocytes",
        "renal collecting duct intercalated cells",
        "renal collecting duct principal cells", "renal connecting tubule cells",
    ],
    "brain": [
        "astrocytes", "brain excitatory neurons", "brain inhibitory neurons",
        "oligodendrocytes", "oligodendrocyte progenitor cells", "microglia",
        "bergmann glia", "other brain neurons", "ependymal cells",
        "choroid plexus epithelial cells",
    ],
    "adrenal gland": ["adrenal cortex cells", "adrenal medulla cells"],
    "lymph node": [
        "b-cells", "t-cells", "nk-cells", "plasma cells", "innate lymphoid cells",
        "lymphatic endothelial cells", "mast cells",
    ],
    "skin": [
        "basal keratinocytes", "suprabasal keratinocytes", "melanocytes", "fibroblasts",
    ],
    "sex organ": [
        "sertoli cells", "leydig cells", "peritubular myoid cells",
        "undifferentiated spermatogonia", "differentiating spermatogonia",
        "early primary spermatocytes", "late primary spermatocytes",
        "early spermatids", "late spermatids", "granulosa cells", "oocytes",
        "ovarian stromal cells", "epididymal basal cells", "epididymal clear cells",
        "epididymal principal cells", "epididymal efferent duct absorptive cells",
        "epididymal efferent duct ciliated cells",
    ],
    "lung": [
        "alveolar cells type 1", "alveolar cells type 2", "transitional alveolar cells",
        "respiratory basal cells", "respiratory ciliated cells",
        "respiratory deuterosomal cells", "respiratory ionocytes",
        "respiratory secretory cells",
    ],
    "blood": [
        "erythrocytes", "erythrocyte progenitors", "platelets", "megakaryocytes",
        "megakaryocyte progenitors", "megakaryocyte-erythroid progenitors",
        "neutrophils", "neutrophil progenitors", "monocytes", "monocyte progenitors",
        "hematopoietic stem cells", "t-cells", "b-cells", "nk-cells",
    ],
}


def log(msg: str) -> None:
    print(msg)
    with open(OUT / "g8_log.txt", "a", encoding="utf-8") as fh:
        fh.write(msg + "\n")


def symbols_from_se_qc() -> tuple[dict[str, str], int]:
    """row-index (as str) -> gene symbol, using the recovered matrix rownames."""
    ids = [l.strip() for l in ROWNAMES.read_text(encoding="utf-8").splitlines() if l.strip()]
    conv = pd.read_csv(GIT / "gene_conversion.csv")
    conv.columns = ["idx", "ens_version", "ens", "symbol"]
    ens2sym = {e: s for e, s in zip(conv["ens_version"], conv["symbol"]) if isinstance(s, str)}
    mapping = {str(i + 1): ens2sym.get(e, "") for i, e in enumerate(ids)}
    return mapping, sum(1 for v in mapping.values() if v)


def tissue_key(raw: str) -> str | None:
    for pat, key in TISSUE_PATTERNS:
        if re.search(pat, raw):
            return key
    return None


def load_sample_meta() -> pd.DataFrame:
    meta = pd.read_csv(SCAN / "GSE226106_samples.tsv", sep="\t", low_memory=False)
    meta = meta.set_index("field").T.reset_index().rename(columns={"index": "column"})
    rows = []
    for _, r in meta.iterrows():
        col = str(r["column"]).strip().strip('"')
        title = str(r.get("Sample_title", ""))
        source = str(r.get("Sample_source_name_ch1", ""))
        chars = str(r.get("Sample_characteristics_ch1", ""))
        # the count matrix is keyed by the library name inside the last [ ... ]
        brackets = re.findall(r"\[([^\]]+)\]", title)
        library = brackets[-1] if brackets else title.split(" ")[0]
        animal = title.split("_")[0].split(".")[0].strip()
        raw_tissue = title[len(animal):].lstrip("_. ")
        raw_tissue = raw_tissue.split(" [")[0].strip().lstrip(".") or source
        raw_tissue = re.sub(r"^D-?\d+_", "", raw_tissue)          # Zyagen_D000_Adrenal...
        raw_tissue = re.sub(r"[-_]D(m)?-?\d+.*$", "", raw_tissue)  # WholeBlood-Dm30b
        raw_tissue = re.sub(r"[-_](BL|NEC|D0)$", "", raw_tissue)
        m = re.search(r"D(-?\d+)", chars)
        day = int(m.group(1)) if m else np.nan
        rows.append(
            {
                "column": library,
                "sample_column": col,
                "gsm": r.get("Sample_geo_accession", ""),
                "title": title,
                "animal": animal,
                "raw_tissue": raw_tissue,
                "source": source,
                "day": day,
                "tissue": tissue_key(raw_tissue),
                "is_vendor": animal.lower().startswith("zyagen"),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- expression


def load_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    sym, n_mapped = symbols_from_se_qc()
    with gzip.open(SCAN / "GSE226106_counts.txt.gz", "rt") as fh:
        counts = pd.read_csv(fh, sep="\t", index_col=0)
    counts.columns = [c.strip().strip('"') for c in counts.columns]
    counts.index = counts.index.astype(str)
    counts = counts[counts.sum(axis=1) > 0]

    cpm = counts / counts.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1.0)
    logcpm.index = [sym.get(i, "") for i in logcpm.index]
    logcpm = logcpm[logcpm.index != ""]
    logcpm = logcpm.groupby(level=0).mean()
    log(f"gene symbols recovered from se.qc.rds rownames: {n_mapped} rows mapped; "
        f"{logcpm.shape[0]} unique symbols enter scoring")

    # Offset is estimated on high-coverage genes (the sample's global level),
    # then applied to every gene so that low-baseline ISGs (ISG15, MX1) survive.
    # W17 lesson: remove each sample's global level before scoring, otherwise
    # library complexity masquerades as an ISG signal.
    detectable = logcpm >= np.log2(1 + 0.5)
    common = logcpm.index[detectable.sum(axis=1) >= 0.5 * logcpm.shape[1]]
    offset = logcpm.loc[common].median(axis=0)
    centred = logcpm.sub(offset, axis=1)
    expressed = logcpm.index[detectable.sum(axis=1) >= 0.1 * logcpm.shape[1]]
    z = centred.loc[expressed].sub(
        centred.loc[expressed].mean(axis=1), axis=0
    ).div(centred.loc[expressed].std(axis=1).replace(0, np.nan), axis=0)
    missing = [g for g in set(ISG_CORE + ISG_BROAD + IFN_REC) if g not in z.index]
    log(f"panel genes absent from macaque annotation (expression filter): {missing}")
    log(f"genes retained for scoring: {z.shape[0]} of {logcpm.shape[0]}")
    depth = pd.Series((counts > 0).sum(axis=0), index=counts.columns, name="genes_detected")
    return z, depth, logcpm


def panel_score(z: pd.DataFrame, genes: list[str], col: str) -> float:
    present = [g for g in genes if g in z.index and np.isfinite(z.loc[g, col])]
    return float(z.loc[present, col].mean()) if present else np.nan


# ------------------------------------------------------------------- Delta(t)


def delta_by_tissue() -> pd.DataFrame:
    grad = pd.read_csv(GRAD, sep="\t")
    wide = pd.read_csv(WIDE, sep="\t")
    cts = set(wide["cell_type"])
    rows = []
    for tissue, names in CELLTYPE_TISSUE.items():
        matched = [c for c in cts if c in names]
        for c in names:  # tolerant match for names mangled by encoding (e.g. Muller glia)
            if c.endswith("glia"):
                matched += [x for x in cts if x.endswith("ller glia")]
        matched = sorted(set(matched))
        sub = grad[grad["cell_type"].isin(matched)]
        rows.append(
            {
                "tissue": tissue,
                "n_cell_types": len(matched),
                "delta_median": float(sub["delta_restriction"].median()) if len(sub) else np.nan,
                "cell_types": ";".join(matched),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------- statistics


def perm_p(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> float:
    if len(x) < 4:
        return float("nan")
    obs = abs(spearmanr(x, y).statistic)
    cnt = sum(
        abs(spearmanr(x, rng.permutation(y)).statistic) >= obs for _ in range(N_PERM)
    )
    return (cnt + 1) / (N_PERM + 1)


def boot_ci(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    if len(x) < 4:
        return (np.nan, np.nan)
    vals = []
    n = len(x)
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        if len(np.unique(x[idx])) < 3 or len(np.unique(y[idx])) < 3:
            continue
        vals.append(spearmanr(x[idx], y[idx]).statistic)
    if not vals:
        return (np.nan, np.nan)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def loo(x: np.ndarray, y: np.ndarray) -> dict:
    vals = []
    for i in range(len(x)):
        xs, ys = np.delete(x, i), np.delete(y, i)
        if len(np.unique(xs)) < 3 or len(np.unique(ys)) < 3:
            continue
        vals.append(spearmanr(xs, ys).statistic)
    if not vals:
        return {"loo_min": np.nan, "loo_max": np.nan, "loo_n_negative": 0}
    return {
        "loo_min": float(np.min(vals)),
        "loo_max": float(np.max(vals)),
        "loo_n_negative": int(sum(1 for v in vals if v < 0)),
    }


def residualise(values: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    mask = np.isfinite(values) & np.isfinite(covariate)
    out = np.full_like(np.asarray(values, dtype=float), np.nan)
    if mask.sum() < 4:
        return out
    beta = np.polyfit(covariate[mask], values[mask], 1)
    out[mask] = values[mask] - np.polyval(beta, covariate[mask])
    return out


def main() -> None:
    (OUT / "g8_log.txt").write_text("", encoding="utf-8")
    meta = load_sample_meta()
    z, depth, logcpm = load_expression()
    meta = meta[meta["column"].isin(z.columns)].copy()
    meta["ISG_core"] = [panel_score(z, ISG_CORE, c) for c in meta["column"]]
    meta["ISG_broad"] = [panel_score(z, ISG_BROAD, c) for c in meta["column"]]
    meta["IFN_rec"] = [panel_score(z, IFN_REC, c) for c in meta["column"]]
    meta["genes_detected"] = meta["column"].map(depth)
    for gene in ("ISG15", "MX1"):
        meta[gene] = z.loc[gene, meta["column"]].to_numpy() if gene in z.index else np.nan

    # annotation audit: canonical markers must land where biology says they do
    audit = []
    zmeta = meta.set_index("column")
    for gene, tissue, direction in (("ALB", "liver", "+"), ("GFAP", "brain", "+"),
                                    ("ISG15", None, "+"), ("MX1", None, "+"),
                                    ("GAPDH", None, "0")):
        if gene not in z.index:
            audit.append({"gene": gene, "check": "missing from matrix"})
            continue
        series = z.loc[gene]
        rec = {"gene": gene, "tissue": tissue or "infection-induction"}
        if tissue:
            sel = zmeta.index[(zmeta["day"] >= 1) & (zmeta["tissue"] == tissue)]
            inf_t = series[sel].median()
            oth = series[zmeta.index[(zmeta["day"] >= 1) & (zmeta["tissue"] != tissue)]].median()
            rec["value"] = round(float(inf_t - oth), 3)
            rec["expect"] = direction
        else:
            a = series[zmeta.index[zmeta["day"] >= 4]].mean()
            b = series[zmeta.index[zmeta["day"] <= 0]].mean()
            rec["value"] = round(float(a - b), 3)
            rec["expect"] = direction
        audit.append(rec)
    log("\n== annotation audit (capped z-scores) ==")
    log(pd.DataFrame(audit).to_string(index=False))

    meta.to_csv(OUT / "sample_scores.tsv", sep="\t", index=False)

    design = (meta[meta["tissue"].notna()]
              .groupby(["tissue", "is_vendor"])
              .agg(n=("column", "size"),
                   n_infected=("day", lambda s: int((s >= 1).sum())),
                   n_control=("day", lambda s: int((s <= 0).sum())),
                   n_animals=("animal", "nunique"))
              .reset_index())
    log("== sample design (tissue x vendor) ==")
    log(design.to_string(index=False))

    # ------------------------------------------------------------- E(t) x4
    inf = meta[(meta["day"] >= 1) & meta["tissue"].notna()].copy()
    ctrl = meta[(meta["day"] <= 0) & meta["tissue"].notna()].copy()
    blood = {"whole blood", "PBMC"}
    inf_tissue = inf[~inf["tissue"].isin(blood)]
    ctrl_tissue = ctrl[~ctrl["tissue"].isin(blood)]

    def agg(df: pd.DataFrame, name: str) -> pd.Series:
        return df.groupby("tissue")[name].median()

    # (b) day-adjusted tissue coefficients from OLS with tissue dummies
    def day_adjusted(df: pd.DataFrame, name: str) -> pd.Series:
        d = df[["tissue", "day", name]].dropna()
        if d["tissue"].nunique() < 2 or d.shape[0] < 6:
            return pd.Series(dtype=float)
        dummies = pd.get_dummies(d["tissue"], prefix="t", drop_first=False).astype(float)
        X = np.column_stack([dummies.to_numpy(), d["day"].to_numpy()])
        beta, *_ = np.linalg.lstsq(X, d[name].to_numpy(), rcond=None)
        mean_day = float(d["day"].mean())
        return pd.Series(
            beta[: dummies.shape[1]] + beta[-1] * mean_day, index=dummies.columns
        ).rename(lambda c: c[2:])

    # (d) within-animal relative score
    def within_animal(df: pd.DataFrame, name: str) -> pd.Series:
        d = df[["tissue", "animal", name]].dropna().copy()
        d["rel"] = d[name] - d.groupby("animal")[name].transform("median")
        return d.groupby("tissue")["rel"].median()

    env = pd.DataFrame(
        {
            "E_raw": agg(inf_tissue, "ISG_core"),
            "E_broad": agg(inf_tissue, "ISG_broad"),
            "E_rec": agg(inf_tissue, "IFN_rec"),
            "E_dayadj": day_adjusted(inf_tissue, "ISG_core"),
            "E_late": agg(inf_tissue[inf_tissue["day"] >= LATE_DAY], "ISG_core"),
            "E_within_animal": within_animal(inf_tissue, "ISG_core"),
            "E_ctrl": agg(ctrl_tissue, "ISG_core"),
            "E_broad_ctrl": agg(ctrl_tissue, "ISG_broad"),
            "n_infected": inf_tissue.groupby("tissue")["column"].size(),
            "n_late": inf_tissue[inf_tissue["day"] >= LATE_DAY].groupby("tissue")["column"].size(),
            "n_control": ctrl_tissue.groupby("tissue")["column"].size(),
            "depth_median": inf_tissue.groupby("tissue")["genes_detected"].median(),
        }
    )
    env["E_induced"] = env["E_late"] - env["E_ctrl"]
    env["E_induced_broad"] = (
        agg(inf_tissue[inf_tissue["day"] >= LATE_DAY], "ISG_broad") - env["E_broad_ctrl"]
    )

    delta = delta_by_tissue()
    merged = (delta.set_index("tissue")
                    .join(env, how="inner")
                    .reset_index())
    merged.to_csv(OUT / "tissue_env.tsv", sep="\t", index=False)
    delta.to_csv(OUT / "tissue_delta.tsv", sep="\t", index=False)

    log("\n== tissue-level merged table ==")
    log(merged[["tissue", "n_cell_types", "delta_median", "E_raw", "E_late",
                "E_ctrl", "E_induced", "E_within_animal", "n_infected",
                "n_control", "depth_median"]].to_string(index=False))

    # --------------------------------------------------------- correlations
    rng = np.random.default_rng(SEED)
    outcomes = {
        "E_raw": merged["E_raw"].to_numpy(),
        "E_dayadj": merged["E_dayadj"].to_numpy(),
        "E_late": merged["E_late"].to_numpy(),
        "E_within_animal": merged["E_within_animal"].to_numpy(),
        "E_induced": merged["E_induced"].to_numpy(),
        "E_induced_broad": merged["E_induced_broad"].to_numpy(),
    }
    delta_v = merged["delta_median"].to_numpy()
    rows = []
    for name, y in outcomes.items():
        mask = np.isfinite(delta_v) & np.isfinite(y)
        if mask.sum() < 4:
            rows.append({"outcome": name, "n": int(mask.sum()), "rho": np.nan,
                         "p_perm": np.nan, "ci_low": np.nan, "ci_high": np.nan})
            continue
        x, yy = delta_v[mask], y[mask]
        rho = float(spearmanr(x, yy).statistic)
        ci = boot_ci(x, yy, rng)
        rows.append(
            {
                "outcome": name,
                "n": int(mask.sum()),
                "rho": round(rho, 3),
                "p_perm": round(perm_p(x, yy, rng), 4),
                "ci_low": round(ci[0], 3),
                "ci_high": round(ci[1], 3),
                **loo(x, yy),
            }
        )
    corr = pd.DataFrame(rows)
    corr.to_csv(OUT / "g8_correlations.tsv", sep="\t", index=False)
    log("\n== Delta vs tissue IFN environment ==")
    log(corr.to_string(index=False))

    # ------------------------------------------------------- depth controls
    depth_rows = []
    for col in ("E_raw", "E_late", "E_induced"):
        d = merged[[col, "depth_median"]].dropna()
        if len(d) >= 4:
            rho, p = spearmanr(d[col], d["depth_median"])
            depth_rows.append({"score": col, "rho_vs_depth": round(float(rho), 3),
                               "p": round(float(p), 4), "n": int(len(d))})
    dep = pd.DataFrame(depth_rows)
    dep.to_csv(OUT / "g8_depth_control.tsv", sep="\t", index=False)

    resid_rows = []
    for name, y in outcomes.items():
        mask = np.isfinite(delta_v) & np.isfinite(y) & np.isfinite(merged["depth_median"].to_numpy())
        if mask.sum() < 5:
            continue
        y_res = residualise(y, merged["depth_median"].to_numpy())
        x_res = residualise(delta_v, merged["depth_median"].to_numpy())
        m2 = np.isfinite(x_res) & np.isfinite(y_res)
        rho = float(spearmanr(x_res[m2], y_res[m2]).statistic)
        resid_rows.append({"outcome": name, "n": int(m2.sum()),
                           "rho_depth_residual": round(rho, 3),
                           "p_spearman": round(float(spearmanr(x_res[m2], y_res[m2]).pvalue), 4)})
    resid = pd.DataFrame(resid_rows)
    resid.to_csv(OUT / "g8_depth_residual.tsv", sep="\t", index=False)
    log("\n== depth control (E vs genes_detected) ==")
    log(dep.to_string(index=False) if len(dep) else "(no depth data)")
    log("\n== after residualising both axes on tissue depth ==")
    log(resid.to_string(index=False) if len(resid) else "(insufficient n)")

    # ------------------------------------------------------- concordance table
    def tertile(series: pd.Series) -> pd.Series:
        s = series.dropna()
        if s.nunique() < 3:
            return pd.Series("na", index=series.index)
        q = pd.qcut(s, 3, labels=["low", "mid", "high"])
        return q.reindex(series.index).astype(object).fillna("na")

    conc = merged[["tissue", "n_cell_types", "delta_median", "E_raw", "E_late",
                   "E_induced", "E_within_animal", "n_infected", "n_control"]].copy()
    conc["delta_tertile"] = tertile(conc["delta_median"])
    conc["env_tertile"] = tertile(conc["E_raw"])
    conc["env_induced_tertile"] = tertile(conc["E_induced"])
    conc["delta_rank"] = conc["delta_median"].rank(ascending=False)
    conc["env_rank"] = conc["E_raw"].rank(ascending=False)

    flip = {"high": "low", "low": "high", "mid": "mid", "na": "na"}
    conc["conflict_vs_E_raw"] = [
        (a in ("high", "low")) and (b == flip[a]) for a, b in zip(conc["delta_tertile"], conc["env_tertile"])
    ]
    conc["concordant_vs_E_raw"] = [
        (a in ("high", "low")) and (b == a) for a, b in zip(conc["delta_tertile"], conc["env_tertile"])
    ]
    conc.to_csv(OUT / "g8_concordance.tsv", sep="\t", index=False)
    log("\n== concordance / conflict (tertiles) ==")
    log(conc.to_string(index=False))

    primary = corr[corr["outcome"] == "E_raw"].iloc[0]
    induced = corr[corr["outcome"] == "E_induced"].iloc[0]
    residual_primary = (resid[resid["outcome"] == "E_raw"]["rho_depth_residual"].iloc[0]
                        if (resid["outcome"] == "E_raw").any() else np.nan)

    def verdict() -> str:
        n = int(primary["n"])
        if n < 6:
            return "功效不足：组织数 < 6，仅描述性排序"
        rho = primary["rho"]
        if np.isnan(rho):
            return "不可判定"
        if rho >= 0.5 and primary["p_perm"] < 0.05 and (
            np.isnan(residual_primary) or np.sign(residual_primary) == np.sign(rho)
        ):
            return "支持：Δ 可作为感染组织 IFN 环境的组织级代理"
        if rho <= 0:
            return "冲突：Δ 排序与感染组织 IFN 环境排序方向不一致"
        return "混合/不可判定：Δ 只作能力轴，预测限定在两者一致的组织"

    summary = {
        "preregistration": "report/T4_G8_预注册_组织IFN环境与Δ对照_20260926.md",
        "seed": SEED,
        "resource": "GSE226106 (Normandin et al. 2023, PMID 38169842)",
        "gene_annotation": {
            "source": "row names recovered from the study's 01-data/se.qc.rds "
                      "(rdata parser); see analysis/t4_g8_tissue_ifn_env_20260926/ANNOTATION_FIX.md",
            "warning": "github/gene_conversion.csv is a biomaRt write.csv output and must NOT be "
                       "used as a row-id map; all pre-2026-09-26 gene-level results on this "
                       "matrix (W14 module scores, W17 panels) used a mismatched mapping.",
            "n_rows": len(symbols_from_se_qc()[0]),
        },
        "annotation_audit": audit,
        "limitations": [
            "macaque not human",
            "EBOV (Makona) not BDBV -> cannot test the multiplier, only the tissue-level proxy",
            "bulk tissue average, not cell-type resolution",
            "sampling day is not orthogonal to animal course",
        ],
        "not_examined": [
            "blood samples excluded from the primary tissue ranking (reported separately)",
        ],
        "tissues_in_primary_test": merged["tissue"].tolist(),
        "n_tissues": int(merged.shape[0]),
        "correlations": corr.to_dict(orient="records"),
        "depth_control": dep.to_dict(orient="records"),
        "depth_residual": resid.to_dict(orient="records"),
        "conflict_tissues_vs_E_raw": conc.loc[conc["conflict_vs_E_raw"], "tissue"].tolist(),
        "concordant_tissues_vs_E_raw": conc.loc[conc["concordant_vs_E_raw"], "tissue"].tolist(),
        "primary": {"outcome": "E_raw", **primary.to_dict()},
        "induced": {"outcome": "E_induced", **induced.to_dict()},
        "verdict": verdict(),
    }
    (OUT / "g8_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    log(f"\nVERDICT: {summary['verdict']}")
    log(f"conflict tissues vs E_raw: {summary['conflict_tissues_vs_E_raw']}")
    log(f"concordant tissues vs E_raw: {summary['concordant_tissues_vs_E_raw']}")


if __name__ == "__main__":
    main()
