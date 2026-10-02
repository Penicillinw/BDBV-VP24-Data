"""Step 3: ISG-module levels and folds in GSE342661, type I vs type III IFN.

Frozen before running (see FINDINGS.md for the same text):

  Modules (copied verbatim from the manuscript's frozen definitions):
    ISG17   = ISG15 MX1 MX2 OAS1 OAS2 IFIT1 IFIT2 IFIT3 BST2 RSAD2 USP18
              STAT1 IRF7 IFI27 IFI6 DDX58 IFIH1
    HK      = ACTB GAPDH TBP RPL13A PPIA HPRT1            (negative control)
  Sample score = mean over available module genes of log2(TPM + 1).
  Baseline for a treated sample = mean score of the Ctrl samples with the same
    time point and the same route (n = 5).
  Readings:
    level = score of the treated sample
    fold  = score of the treated sample - its paired Ctrl score (same rep,
            time, route)
  Tests: exact sign-flip (binomial) and Wilcoxon signed-rank (exact) on the
    paired differences, plus paired t as a sensitivity. n and the smallest
    attainable two-sided p are reported for every test.
  Nulls: (i) gene-module null - the type I vs type III difference recomputed on
    2,000 random 17-gene modules drawn from expressed genes, stratified by mean
    expression decile; (ii) label null - treatment labels permuted within
    (time, route) blocks.

Outputs (this directory only): sample_scores.tsv, paired_tests.tsv,
g342_module.json, module_null.tsv
"""

import gzip
import json
import math
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
RAW = HERE / "raw"
MATRIX = RAW / "GSE342661_counts_tpm.matrix.gz"
HGNC = RAW / "hgnc_complete_set.txt"

ISG17 = ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "IFIT1", "IFIT2", "IFIT3",
         "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6",
         "DDX58", "IFIH1"]
HK = ["ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "HPRT1"]
RECEPTOR = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB"]
PRIMING = ["ISG15", "MX1"]

BASELINE_FLOOR = 0.0  # gene must have mean log2(TPM+1) >= 0 in at least one group
N_PERM = 2000
RNG_SEED = 20260926


def load_matrix():
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        df = pd.read_csv(fh, sep="\t", quotechar='"', index_col=0, low_memory=False)
    df.index.name = "ensembl_versioned"
    df = df.astype("float64")
    return df


def ensembl_to_symbol():
    h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
    h = h[h["ensembl_gene_id"].notna() & h["symbol"].notna()]
    # one symbol per Ensembl id: keep Approved rows, drop any remaining ambiguity
    h["rank"] = (h["status"] != "Approved").astype(int)
    h = h.sort_values("rank")
    m = {}
    for ens, sym in zip(h["ensembl_gene_id"], h["symbol"]):
        m.setdefault(ens, sym)
    counts = h.groupby("ensembl_gene_id")["symbol"].nunique()
    ambiguous = set(counts[counts > 1].index)
    # alias / previous symbols -> Ensembl id, so that legacy module names such as
    # DDX58 (now RIGI) resolve to the row the matrix actually carries
    alias = {}
    for col in ("prev_symbol", "alias_symbol"):
        for names, ens in zip(h[col].fillna(""), h["ensembl_gene_id"]):
            for nm in str(names).replace("|", ";").split(";"):
                nm = nm.strip().strip('"')
                if nm:
                    alias.setdefault(nm, ens)
    return m, ambiguous, alias


def sign_flip_p(diffs):
    d = np.asarray(diffs, dtype=float)
    n = len(d)
    pos = int((d > 0).sum())
    neg = int((d < 0).sum())
    k = max(pos, neg)
    p = 2.0 * sum(math.comb(n, i) for i in range(k, n + 1)) / (2 ** n)
    return min(1.0, p), pos, neg


def paired_block(a, b, label, extra=None):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = a - b
    if len(d) < 1:
        return None
    out = {
        "test": label,
        "n_pairs": int(len(d)),
        "mean_a": float(a.mean()),
        "mean_b": float(b.mean()),
        "mean_diff": float(d.mean()),
        "sd_diff": float(d.std(ddof=1)) if len(d) > 1 else float("nan"),
        "cohens_dz": float(d.mean() / d.std(ddof=1)) if len(d) > 1 and d.std(ddof=1) > 0 else float("nan"),
        "min_attainable_two_sided_p_signflip": 2.0 / (2 ** len(d)),
    }
    p_sf, pos, neg = sign_flip_p(d)
    out["sign_flip_p"] = p_sf
    out["n_positive"] = pos
    out["n_negative"] = neg
    try:
        w = stats.wilcoxon(d, method="exact")
        out["wilcoxon_exact_p"] = float(w.pvalue)
    except Exception:
        out["wilcoxon_exact_p"] = float("nan")
    try:
        t = stats.ttest_rel(a, b)
        out["paired_t_p"] = float(t.pvalue)
    except Exception:
        out["paired_t_p"] = float("nan")
    if extra:
        out.update(extra)
    return out


def main() -> None:
    design = pd.read_csv(HERE / "design.tsv", sep="\t")
    assert len(design) == 180
    tpm = load_matrix()
    symbols, ambiguous, alias = ensembl_to_symbol()
    base_ids = pd.Index([i.split(".")[0] for i in tpm.index])
    sym = pd.Index([symbols.get(i, None) for i in base_ids])
    keep = np.asarray(pd.notna(sym)) & ~np.asarray([i in ambiguous for i in base_ids])
    mapped = tpm.loc[keep].copy()
    mapped.index = sym[keep]
    # collapse duplicates by max TPM (report how many)
    dup_syms = int(mapped.index.duplicated().sum())
    mapped = mapped.groupby(level=0).max()
    # symbol-level resolver: approved symbol first, then previous/alias symbols
    current_symbols = set(mapped.index)
    ens_of_symbol = {}
    for base, s in zip(base_ids, sym):
        if s is not None and s in current_symbols:
            ens_of_symbol.setdefault(s, base)

    def resolve(gene):
        """Return (row label present in the matrix, route) or (None, route)."""
        if gene in current_symbols:
            return gene, "approved symbol"
        ens = ens_of_symbol.get(gene) or alias.get(gene)
        if ens is None:
            return None, "unresolved"
        candidate = symbols.get(ens)
        if candidate in current_symbols:
            return candidate, f"alias/previous symbol -> {candidate}"
        return None, "unresolved"

    resolved = {}  # module gene -> row label actually used
    for genes in (ISG17, HK, RECEPTOR, PRIMING):
        for g in genes:
            if g in resolved:
                continue
            resolved[g] = resolve(g)

    coverage = {
        "matrix_rows": int(len(tpm)),
        "rows_with_hgnc_symbol": int(keep.sum()),
        "mapping_rate": float(keep.sum() / len(tpm)),
        "duplicate_symbols_collapsed": dup_syms,
        "ambiguous_ensembl_dropped": int(sum(i in ambiguous for i in base_ids)),
    }
    for name, genes in (("ISG17", ISG17), ("HK", HK), ("RECEPTOR", RECEPTOR),
                        ("PRIMING", PRIMING)):
        present = [g for g in genes if resolved[g][0] is not None]
        coverage[f"{name}_present"] = present
        coverage[f"{name}_missing"] = [g for g in genes if resolved[g][0] is None]
        coverage[f"{name}_alias_routes"] = {
            g: resolved[g][1] for g in genes if resolved[g][0] is not None
            and resolved[g][0] != g
        }

    log2tpm = np.log2(mapped + 1.0)

    def module_score(genes):
        present = [resolved[g][0] for g in genes if resolved[g][0] is not None]
        return log2tpm.loc[present].mean(axis=0), present

    scores = {}
    for name, genes in (("ISG17", ISG17), ("HK", HK), ("RECEPTOR", RECEPTOR),
                        ("PRIMING", PRIMING)):
        s, present = module_score(genes)
        scores[name] = s
    samples = pd.DataFrame(scores)
    samples.index.name = "column"
    samples = samples.reset_index().merge(design, on="column", how="left")
    assert samples["gsm"].notna().all()
    samples.to_csv(HERE / "sample_scores.tsv", sep="\t", index=False)

    # baseline (Ctrl) per time x route x replicate; each replicate has its own
    # untreated well, so the baseline is the paired Ctrl sample, not a group mean
    ctrl_score = (
        samples[samples["treatment"] == "Ctrl"]
        .set_index(["time", "route", "rep"])["ISG17"]
    )
    samples["baseline_ISG17"] = [
        ctrl_score.get((t, r, p), float("nan"))
        for t, r, p in zip(samples["time"], samples["route"], samples["rep"])
    ]
    samples["fold_ISG17"] = samples["ISG17"] - samples["baseline_ISG17"]
    samples.to_csv(HERE / "sample_scores.tsv", sep="\t", index=False)

    treated = samples[samples["treatment"] != "Ctrl"].copy()
    results = []

    # --- induction vs matched control, per treatment x time x route (n = 5) ---
    for (t, r, trt), grp in treated.groupby(["time", "route", "treatment"]):
        grp = grp.sort_values("rep")
        rec = paired_block(
            grp["ISG17"].values,
            grp["baseline_ISG17"].values,
            f"level_induced_vs_ctrl|{t}|{r}|{trt}",
        )
        rec.update({"kind": "induction_level", "time": t, "route": r, "treatment": trt})
        results.append(rec)
        rec2 = paired_block(
            grp["fold_ISG17"].values,
            np.zeros(len(grp)),
            f"fold_gt_zero|{t}|{r}|{trt}",
        )
        rec2.update({"kind": "fold_vs_zero", "time": t, "route": r, "treatment": trt})
        results.append(rec2)

    # --- type I (IFNB1) vs type III isoforms, paired by replicate ---
    wide = treated.pivot_table(
        index=["time", "route", "rep"], columns="treatment", values="ISG17"
    )
    wide_fold = treated.pivot_table(
        index=["time", "route", "rep"], columns="treatment", values="fold_ISG17"
    )
    for iso in ["IFNL1", "IFNL2", "IFNL3", "IFNL4"]:
        for t in sorted(treated["time"].unique()):
            for r in sorted(treated["route"].unique()):
                a, b = wide.loc[(t, r), "IFNB1"], wide.loc[(t, r), iso]
                rec = paired_block(a.values, b.values, f"level_IFNB1_minus_{iso}|{t}|{r}")
                rec.update({"kind": "level_I_vs_III", "time": t, "route": r, "isoform": iso})
                results.append(rec)
    # pooled across time points, per route, n = 15
    for r in sorted(treated["route"].unique()):
        sub = wide.loc[(slice(None), r), :]
        for iso in ["IFNL1", "IFNL2", "IFNL3", "IFNL4"]:
            rec = paired_block(
                sub["IFNB1"].values, sub[iso].values, f"level_IFNB1_minus_{iso}|pooled|{r}"
            )
            rec.update({"kind": "level_I_vs_III_pooled", "time": "pooled", "route": r, "isoform": iso})
            results.append(rec)
        l123 = sub[["IFNL1", "IFNL2", "IFNL3"]].mean(axis=1)
        rec = paired_block(
            sub["IFNB1"].values, l123.values, f"level_IFNB1_minus_IFNL123_mean|pooled|{r}"
        )
        rec.update({"kind": "level_I_vs_III_pooled", "time": "pooled", "route": r,
                    "isoform": "IFNL123_mean"})
        results.append(rec)
        # fold reading, same contrast
        subf = wide_fold.loc[(slice(None), r), :]
        rec = paired_block(
            subf["IFNB1"].values, subf[["IFNL1", "IFNL2", "IFNL3"]].mean(axis=1).values,
            f"fold_IFNB1_minus_IFNL123_mean|pooled|{r}",
        )
        rec.update({"kind": "fold_I_vs_III_pooled", "time": "pooled", "route": r,
                    "isoform": "IFNL123_mean"})
        results.append(rec)

    # --- route gate: apical-only vs apical-basal, paired by rep within time ---
    route_wide = samples.pivot_table(
        index=["time", "rep"], columns=["treatment", "route"], values="ISG17"
    )
    for trt in ["Ctrl", "IFNB1", "IFNL1", "IFNL2", "IFNL3", "IFNL4"]:
        for t in sorted(samples["time"].unique()):
            try:
                a = route_wide.loc[t, (trt, "apical-basal")]
                b = route_wide.loc[t, (trt, "apical-only")]
            except KeyError:
                continue
            rec = paired_block(a.values, b.values, f"level_ab_minus_a|{t}|{trt}")
            rec.update({"kind": "route_gate", "time": t, "route": "ab_vs_a", "treatment": trt})
            results.append(rec)

    tests = pd.DataFrame([r for r in results if r])
    front = ["kind", "test", "time", "route", "treatment", "isoform", "n_pairs",
             "mean_a", "mean_b", "mean_diff", "cohens_dz", "sign_flip_p",
             "n_positive", "n_negative", "wilcoxon_exact_p", "paired_t_p",
             "min_attainable_two_sided_p_signflip"]
    tests = tests[[c for c in front if c in tests.columns] + [c for c in tests.columns if c not in front]]
    tests.to_csv(HERE / "paired_tests.tsv", sep="\t", index=False)

    # --- gene-module null for the primary contrast (pooled apical-basal, n=15) ---
    rng = np.random.default_rng(RNG_SEED)
    expressed = log2tpm.mean(axis=1)
    expressed = expressed[expressed > BASELINE_FLOOR]
    deciles = pd.qcut(expressed, 10, labels=False, duplicates="drop")
    genes_by_decile = {d: list(g.index) for d, g in expressed.groupby(deciles)}

    def fold_frame(score_by_column, design_df):
        s = pd.DataFrame({"column": score_by_column.index, "score": score_by_column.values})
        s = s.merge(design_df, on="column", how="left")
        ctrl = (
            s[s["treatment"] == "Ctrl"]
            .groupby(["time", "route", "rep"])["score"].mean()
        )
        keys = list(zip(s["time"], s["route"], s["rep"]))
        s["baseline"] = [ctrl.get(k, float("nan")) for k in keys]
        s["fold"] = s["score"] - s["baseline"]
        return s

    def primary_stat(design_df, score_by_column):
        s = fold_frame(score_by_column, design_df)
        sub = s[s["treatment"].isin(["IFNB1", "IFNL1", "IFNL2", "IFNL3"])]
        piv = sub.pivot_table(index=["time", "route", "rep"], columns="treatment",
                              values="fold")
        if not {"IFNB1", "IFNL1", "IFNL2", "IFNL3"}.issubset(piv.columns):
            return float("nan")
        return float((piv["IFNB1"] - piv[["IFNL1", "IFNL2", "IFNL3"]].mean(axis=1)).mean())

    def module_stat(genes):
        rows = [g if g in log2tpm.index else resolved.get(g, (None, ""))[0] for g in genes]
        rows = [r for r in rows if r is not None and r in log2tpm.index]
        s = log2tpm.loc[rows].mean(axis=0)
        return primary_stat(design, s)

    observed = module_stat(ISG17)
    nulls = []
    decile_keys = list(genes_by_decile)
    # draw 17 genes in proportion to decile sizes, so the null module matches the
    # expression profile of ISG17 rather than a flat genome-wide draw
    weights = np.array([len(genes_by_decile[d]) for d in decile_keys], float)
    weights = weights / weights.sum()
    for _ in range(N_PERM):
        drawn = []
        for d in rng.choice(decile_keys, size=17, replace=True, p=weights):
            drawn.append(genes_by_decile[d][rng.integers(len(genes_by_decile[d]))])
        drawn = list(dict.fromkeys(drawn))
        if len(drawn) < 17:
            continue
        nulls.append(module_stat(drawn))
    nulls = np.array(nulls)
    pd.DataFrame({"null_type_I_minus_III_fold_diff": nulls}).to_csv(
        HERE / "module_null.tsv", sep="\t", index=False
    )
    pct = float((nulls <= observed).mean() * 100)

    # --- label null: permute treatment labels within (time, route) blocks ---
    isg_score = pd.Series(
        log2tpm.loc[[resolved[g][0] for g in ISG17 if resolved[g][0] is not None]]
        .mean(axis=0)
        .values,
        index=log2tpm.columns,
    )
    label_null = []
    for _ in range(N_PERM):
        shuffled = design.copy()
        for _, grp in shuffled.groupby(["time", "route"]):
            vals = np.array(list(grp["treatment"].astype(object).values), dtype=object)
            rng.shuffle(vals)
            shuffled.loc[grp.index, "treatment"] = vals
        value = primary_stat(shuffled, isg_score)
        if not math.isnan(value):
            label_null.append(value)
    label_null = np.array(label_null)
    lpct = float((label_null <= observed).mean() * 100)

    summary = {
        "coverage": coverage,
        "ISG17_present": coverage["ISG17_present"],
        "observed_primary_statistic": observed,
        "primary_statistic_definition": (
            "mean over all (time, route, replicate) of "
            "fold(IFNB1) - mean fold(IFNL1, IFNL2, IFNL3); "
            "fold = ISG17 level minus matched Ctrl level"
        ),
        "module_null": {
            "n": int(len(nulls)),
            "median": float(np.median(nulls)),
            "q2.5": float(np.percentile(nulls, 2.5)),
            "q97.5": float(np.percentile(nulls, 97.5)),
            "percentile_of_observed": pct,
            "two_sided_empirical_p": float(2 * min((nulls <= observed).mean(), (nulls >= observed).mean())),
        },
        "label_null": {
            "n": int(len(label_null)),
            "median": float(np.median(label_null)),
            "q2.5": float(np.percentile(label_null, 2.5)),
            "q97.5": float(np.percentile(label_null, 97.5)),
            "percentile_of_observed": lpct,
        },
        "ISG17_mean_by_condition": {
            f"{t}|{r}|{trt}": float(v)
            for (t, r, trt), v in samples.groupby(["time", "route", "treatment"])["ISG17"].mean().items()
        },
        "ISG17_fold_mean_by_condition": {
            f"{t}|{r}|{trt}": float(v)
            for (t, r, trt), v in samples.groupby(["time", "route", "treatment"])["fold_ISG17"].mean().items()
        },
        "HK_mean_by_condition": {
            f"{t}|{r}|{trt}": float(v)
            for (t, r, trt), v in samples.groupby(["time", "route", "treatment"])["HK"].mean().items()
        },
    }
    (HERE / "g342_module.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ("ISG17_mean_by_condition", "ISG17_fold_mean_by_condition",
                                   "HK_mean_by_condition")}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
