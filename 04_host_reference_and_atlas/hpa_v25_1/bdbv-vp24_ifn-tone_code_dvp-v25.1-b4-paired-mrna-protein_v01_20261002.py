"""B4 / A5 - paired quantitative mRNA-protein audit for the twelve genes of Delta.

Task
----
The manuscript states: "Of the twelve genes it uses, five carry no tissue-level
protein annotation, the antibody for IFNLR1 is rated uncertain, and the KPNA5
and IRF9 protein annotations point in a different direction from their
single-cell mRNA."

That sentence came from a *qualitative* comparison (HPA IHC annotation text vs
HPA single-cell mRNA; see analysis/a2_layer_consistency_20260926/).  The frozen
work plan registers the methodological redo (A5): replace it with a *paired,
quantitative* resource and report each gene's mRNA-protein correlation together
with its percentile in the genome-wide null distribution.

FROZEN CRITERIA (recorded before any value from this resource was inspected)
---------------------------------------------------------------------------
R1  Resource.  Wang et al. 2019, Mol Syst Biol 15:e8503 (PMC6379049).
    Table EV1 sheet "C. Genes"  = protein, gene-level intensity, 29 healthy
                                  human tissues, same donors as the RNA.
    Table EV2 sheet "B. Genes"  = mRNA, gene-level TPM, same 29 tissues.
    Both tables are open-access supplementary material retrieved through the
    Europe PMC supplementaryFiles endpoint.
R2  Focal gene set (frozen, identical to the manuscript's twelve):
    IFNAR1, IFNAR2, IFNLR1, IL10RB, ISG15, MX1                (score, 6)
    KPNA1,  KPNA5,  KPNA6                                     (second axis, 3)
    STAT1,  STAT2,  IRF9                                      (cargo, 3)
    No gene may be substituted or dropped silently.
R3  Universe.  All genes present in BOTH tables, matched on Ensembl gene ID.
R4  Statistic.  Spearman rho between mRNA and protein *across tissues*, per gene.
R5  Primary detection rule.  Listwise: use only tissues where protein intensity
    > 0 AND mRNA TPM > 0.  A gene enters the genome-wide null with rho only if
    at least 10 tissues survive this rule.
R6  Sensitivity rules (reported, not substituted for R5):
    S1 all 29 tissues kept, zeros retained (ties dominate the rank statistic);
    S2 use every tissue with mRNA > 0 regardless of protein detection.
R7  Null.  The genome-wide distribution of rho under R4/R5 over all genes that
    pass R5 provides the percentile for each focal gene.  No randomisation is
    needed because the null is the empirical genome-wide distribution.
R8  Multiplicity.  Benjamini-Hochberg across the twelve focal genes, separately
    inside each rule (R5, S1, S2).
R9  Detectability.  For every focal gene also report the number of the 29
    tissues in which protein was not detected (intensity 0).
R10 Reporting.  Every number in FINDINGS.md must be produced by this script.

Run:
    python analysis/b4_mrna_protein_paired_20260926/b4_paired_mrna_protein.py
"""

from __future__ import annotations

import json
import os
from datetime import date

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw", "epmc", "wang2019")
PROT_XLSX = os.path.join(RAW, "Table_EV1.xlsx")
RNA_XLSX = os.path.join(RAW, "Table_EV2.xlsx")

FOCAL = [
    "IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1",
    "KPNA1", "KPNA5", "KPNA6",
    "STAT1", "STAT2", "IRF9",
]
FOCAL_GROUP = {
    "IFNAR1": "score_typeI", "IFNAR2": "score_typeI",
    "IFNLR1": "score_typeIII", "IL10RB": "score_typeIII",
    "ISG15": "score_priming", "MX1": "score_priming",
    "KPNA1": "second_axis", "KPNA5": "second_axis", "KPNA6": "second_axis",
    "STAT1": "cargo", "STAT2": "cargo", "IRF9": "cargo",
}

ANNOT_COLS = {"Tissue enriched", "Group enriched", "Tissue enhanced",
              "Classification"}
MIN_TISSUES = 10
OUT_TABLES = os.path.join(HERE, "tables")


def load_gene_table(path: str, sheet: str) -> pd.DataFrame:
    """Return a gene x tissue matrix indexed by Ensembl gene ID."""
    cache = os.path.join(RAW, os.path.basename(path).replace(".xlsx", f"_{sheet.split('.')[0].strip()}.tsv"))
    if os.path.exists(cache):
        df = pd.read_csv(cache, sep="\t", index_col=0)
        df.index.name = "ensembl_gene_id"
        return df
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    rows = ws.iter_rows(values_only=True)
    header = [str(h).strip() if h is not None else "" for h in next(rows)]
    id_cols = {"Gene ID", "Gene name"}
    tissue_cols = [h for h in header
                   if h and h not in id_cols and h not in ANNOT_COLS]
    idx_id = header.index("Gene ID")
    idx_name = header.index("Gene name")
    keep = [header.index(c) for c in tissue_cols]
    data = {}
    names = {}
    for row in rows:
        gid = row[idx_id]
        if not gid:
            continue
        gid = str(gid).strip()
        if gid in data:
            continue
        vals = []
        ok = True
        for j in keep:
            v = row[j]
            if v is None or (isinstance(v, str) and not v.strip()):
                vals.append(np.nan)
            else:
                try:
                    vals.append(float(v))
                except (TypeError, ValueError):
                    ok = False
                    break
        if not ok:
            continue
        data[gid] = vals
        names[gid] = str(row[idx_name]).strip() if row[idx_name] else ""
    wb.close()
    df = pd.DataFrame.from_dict(data, orient="index", columns=tissue_cols)
    df.insert(0, "gene_name", pd.Series(names))
    df.index.name = "ensembl_gene_id"
    # harmonise tissue labels
    df.columns = [c.strip() for c in df.columns]
    df.to_csv(cache, sep="\t", float_format="%.6g")
    return df


def spearman_xy(x: np.ndarray, y: np.ndarray):
    """Spearman rho, n, two-sided p.  Returns NaN rho when n < 4 or constant."""
    n = len(x)
    if n < 4 or np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan, n, np.nan
    rho, p = stats.spearmanr(x, y)
    return float(rho), int(n), float(p)


def bh(pvals):
    """Benjamini-Hochberg q values, NaN-safe."""
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    m = len(pv)
    ranked = pv[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def analyse(rule_name: str, mask_fn, prot: pd.DataFrame, rna: pd.DataFrame,
            universe, tissues):
    rows = []
    for gid in universe:
        y = prot.loc[gid, tissues].to_numpy(dtype=float)
        x = rna.loc[gid, tissues].to_numpy(dtype=float)
        mask = mask_fn(x, y)
        rho, n, p = spearman_xy(x[mask], y[mask])
        rows.append({"ensembl_gene_id": gid,
                     "gene_name": prot.loc[gid, "gene_name"],
                     "rule": rule_name, "rho": rho, "n_tissues": n, "p": p})
    df = pd.DataFrame(rows)
    null = df.dropna(subset=["rho"])
    null = null[null["n_tissues"] >= MIN_TISSUES]
    df["percentile_in_genome"] = np.nan
    if len(null):
        vals = null["rho"].to_numpy()
        for i, r in df.iterrows():
            if not np.isnan(r["rho"]):
                df.at[i, "percentile_in_genome"] = float(
                    (vals < r["rho"]).sum() + 0.5 * (vals == r["rho"]).sum()
                ) / len(vals)
    focal = df[df["gene_name"].isin(FOCAL)].copy()
    focal["group"] = focal["gene_name"].map(FOCAL_GROUP)
    focal = focal.sort_values("gene_name")
    focal["q_bh_focal"] = bh(focal["p"].to_numpy())
    return df, focal, {
        "rule": rule_name,
        "genes_in_universe": int(len(df)),
        "genes_in_null": int(len(null)),
        "null_median_rho": float(np.median(null["rho"])) if len(null) else None,
        "null_q25_rho": float(np.percentile(null["rho"], 25)) if len(null) else None,
        "null_q75_rho": float(np.percentile(null["rho"], 75)) if len(null) else None,
        "null_median_n_tissues": float(np.median(null["n_tissues"])) if len(null) else None,
    }


def main():
    os.makedirs(OUT_TABLES, exist_ok=True)
    prot_all = load_gene_table(PROT_XLSX, "C. Genes")
    rna_all = load_gene_table(RNA_XLSX, "B. Genes")

    shared_tissues = [t for t in prot_all.columns
                      if t in set(rna_all.columns) and t != "gene_name"]
    universe = [g for g in prot_all.index if g in set(rna_all.index)]
    prot = prot_all.loc[universe, ["gene_name"] + shared_tissues]
    rna = rna_all.loc[universe, shared_tissues]

    summary = {
        "resource": "Wang et al. 2019 Mol Syst Biol 15:e8503 (PMC6379049)",
        "protein_table": os.path.basename(PROT_XLSX) + " / C. Genes",
        "rna_table": os.path.basename(RNA_XLSX) + " / B. Genes",
        "n_tissues_shared": len(shared_tissues),
        "tissues": shared_tissues,
        "n_genes_in_universe": len(universe),
        "min_tissues_rule": MIN_TISSUES,
        "run_date": str(date.today()),
        "rules": {},
    }

    rules = {
        "R5_primary_listwise_detected":
            lambda x, y: (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y),
        "S1_all_tissues_zeros_kept":
            lambda x, y: np.isfinite(x) & np.isfinite(y),
        "S2_mrna_detected_only":
            lambda x, y: (x > 0) & np.isfinite(x) & np.isfinite(y),
    }

    focal_tables = []
    for name, fn in rules.items():
        all_df, focal, stats_rule = analyse(name, fn, prot, rna, universe,
                                            shared_tissues)
        all_df.to_csv(os.path.join(OUT_TABLES, f"b4_genome_{name}.tsv"),
                      sep="\t", index=False, float_format="%.6g")
        focal_tables.append(focal)
        genes = focal.set_index("gene_name")
        stats_rule["focal"] = {
            g: {
                "rho": None if pd.isna(genes.loc[g, "rho"]) else round(float(genes.loc[g, "rho"]), 4),
                "n_tissues": int(genes.loc[g, "n_tissues"]),
                "p": None if pd.isna(genes.loc[g, "p"]) else float(genes.loc[g, "p"]),
                "q_bh_focal": None if pd.isna(genes.loc[g, "q_bh_focal"]) else float(genes.loc[g, "q_bh_focal"]),
                "percentile_in_genome": None if pd.isna(genes.loc[g, "percentile_in_genome"]) else round(float(genes.loc[g, "percentile_in_genome"]), 4),
            }
            for g in FOCAL if g in genes.index
        }
        missing = [g for g in FOCAL if g not in genes.index]
        stats_rule["focal_genes_missing_from_resource"] = missing
        summary["rules"][name] = stats_rule

    focal_all = pd.concat(focal_tables, ignore_index=True)
    focal_all.to_csv(os.path.join(OUT_TABLES, "b4_focal_genes.tsv"),
                     sep="\t", index=False, float_format="%.6g")

    # detectability of the focal proteins across the 29 tissues
    det_rows = []
    for g in FOCAL:
        match = prot[prot["gene_name"] == g]
        if not len(match):
            det_rows.append({"gene": g, "group": FOCAL_GROUP[g],
                             "ensembl_gene_id": None, "tissues_assayed": 0,
                             "tissues_protein_detected": 0,
                             "tissues_protein_undetected": 0,
                             "fraction_undetected": None,
                             "median_mrna_tpm_detected": None})
            continue
        row = match.iloc[0]
        vals = row[shared_tissues].to_numpy(dtype=float)
        det = np.isfinite(vals) & (vals > 0)
        rvals = rna.loc[match.index[0], shared_tissues].to_numpy(dtype=float)
        det_rows.append({
            "gene": g,
            "group": FOCAL_GROUP[g],
            "ensembl_gene_id": match.index[0],
            "tissues_assayed": len(vals),
            "tissues_protein_detected": int(det.sum()),
            "tissues_protein_undetected": int((~det).sum()),
            "fraction_undetected": float((~det).sum() / len(vals)),
            "median_mrna_tpm_detected": float(np.median(rvals[det])) if det.any() else None,
        })
    det_df = pd.DataFrame(det_rows)
    det_df.to_csv(os.path.join(OUT_TABLES, "b4_detectability.tsv"),
                  sep="\t", index=False, float_format="%.6g")
    summary["detectability"] = det_df.to_dict(orient="records")

    with open(os.path.join(HERE, "b4_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False)[:4000])


if __name__ == "__main__":
    main()
