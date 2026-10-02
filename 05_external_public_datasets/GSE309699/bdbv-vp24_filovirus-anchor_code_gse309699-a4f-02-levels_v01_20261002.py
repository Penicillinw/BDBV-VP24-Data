"""A4-filovirus independent recompute, step 2: panel levels per interferon arm.

Frozen criteria (written before any count value was inspected)
-------------------------------------------------------------
C5  The arm -> column map comes **only** from ``geo_column_map.tsv`` produced by
    step 1 (GEO metadata). No treatment label is ever inferred from expression.
C6  Primary normalisation: log2(CPM + 1) with each library's CPM denominator taken as
    the sum over *all* matrix rows (not over the panel), so the denominator cannot be
    inflated by panel selection. Sensitivity: median-of-ratios size factors
    (DESeq-style), reported side by side, never substituted silently.
C7  The ISG panel is the 33-gene Ensembl list inherited **frozen** from the original
    implementation. It is not modified. A 12-gene core subset is declared here, in
    advance, as a *sensitivity* analysis only; if the two disagree the disagreement is
    reported rather than resolved in favour of either.
C8  Duplicate Ensembl rows inside the panel are reported, never silently collapsed
    (the original implementation built a dict keyed by the Ensembl field, which would
    have kept only the last duplicate).
C9  Uncertainty on every arm-level mean is the library-level bootstrap (2000
    resamples, seed 20260926) resampling *libraries*, not genes.
C10 A gene is "present" if the Ensembl ID occurs in the matrix. Missing panel genes
    are reported; the panel mean is taken over the genes actually present, and the
    count is stated every time a mean is quoted.

Reads: ``geo_column_map.tsv`` (this directory), the deposited count matrix (read-only).
Writes: ``arm_levels.tsv``, ``gene_levels.tsv``, ``comparison.json``.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import os
import random
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"G:\本迪布焦研究"
A2 = os.path.join(ROOT, "analysis", "a2_filovirus_anchor_20260926")
COUNTS = os.path.join(A2, "raw", "GSE309699_raw_counts_All_samples.csv.gz")

# --- frozen panel inherited from the original implementation (33 Ensembl IDs) ---
ISG_PANEL = {
    "ENSG00000187608": "ISG15", "ENSG00000157601": "MX1", "ENSG00000183486": "MX2",
    "ENSG00000089127": "OAS1", "ENSG00000111335": "OAS2", "ENSG00000111331": "OAS3",
    "ENSG00000103197": "IFIT1", "ENSG00000119922": "IFIT2", "ENSG00000119917": "IFIT3",
    "ENSG00000134321": "RSAD2", "ENSG00000130303": "BST2", "ENSG00000184979": "USP18",
    "ENSG00000126709": "IFI6", "ENSG00000165949": "IFI27", "ENSG00000185507": "IRF7",
    "ENSG00000107201": "DDX58", "ENSG00000115267": "IFIH1", "ENSG00000115415": "STAT1",
    "ENSG00000170581": "STAT2", "ENSG00000213928": "IRF9", "ENSG00000132530": "XAF1",
    "ENSG00000138646": "HERC5", "ENSG00000138496": "PARP9", "ENSG00000134326": "CMPK2",
    "ENSG00000132274": "TRIM22", "ENSG00000117228": "GBP1", "ENSG00000162645": "GBP2",
    "ENSG00000154451": "GBP5", "ENSG00000125347": "IRF1", "ENSG00000169245": "CXCL10",
    "ENSG00000138755": "CXCL9", "ENSG00000169248": "CXCL11", "ENSG00000131203": "IDO1",
}
# declared in advance, sensitivity only
CORE12 = {"ENSG00000187608": "ISG15", "ENSG00000157601": "MX1", "ENSG00000089127": "OAS1",
          "ENSG00000103197": "IFIT1", "ENSG00000119917": "IFIT3", "ENSG00000134321": "RSAD2",
          "ENSG00000130303": "BST2", "ENSG00000184979": "USP18", "ENSG00000126709": "IFI6",
          "ENSG00000165949": "IFI27", "ENSG00000185507": "IRF7", "ENSG00000115415": "STAT1"}
TYPE_II = ["CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP2", "GBP5"]
TYPE_I = ["MX1", "RSAD2", "IFIT2"]

SEED = 20260926
ARMS = ["IFN alpha", "IFN beta", "IFN gamma", "IFN lambda", "untreated"]


def read_map():
    rows = []
    with open(os.path.join(HERE, "geo_column_map.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    return rows


def arm_of(treatment):
    t = treatment.lower()
    for a in ("alpha", "beta", "gamma", "lambda"):
        if a in t:
            return f"IFN {a}"
    return "untreated"


def load_matrix():
    """Return (columns, ids, values[list of lists]) from one single pass."""
    with gzip.open(COUNTS, "rt", encoding="utf-8-sig", errors="replace") as fh:
        header = next(fh).rstrip("\n").split(",")
        cols = header[2:]
        ids, vals = [], []
        for line in fh:
            f = line.rstrip("\n").split(",")
            if len(f) < 3:
                continue
            ids.append(f[1])
            vals.append([float(x) for x in f[2:]])
    return cols, ids, vals


def logcpm(vals, ncol):
    totals = [0.0] * ncol
    for v in vals:
        for i, x in enumerate(v):
            totals[i] += x
    out = []
    for v in vals:
        out.append([math.log2(v[i] / totals[i] * 1e6 + 1.0) if totals[i] > 0 else 0.0
                    for i in range(ncol)])
    return out, totals


def size_factors_mor(vals, ncol):
    """DESeq-style median-of-ratios size factors (sensitivity normalisation)."""
    logs = []
    for v in vals:
        if all(x > 0 for x in v):
            logs.append([math.log(x) for x in v])

    def geomean(j):
        return st.fmean(row[j] for row in logs)          # log-scale geometric mean

    gm = [geomean(j) for j in range(ncol)]
    sf = [st.median([row[j] - gm[j] for row in logs]) for j in range(ncol)]
    return [math.exp(x) for x in sf]


def main():
    rng = random.Random(SEED)
    mapping = read_map()
    cols, ids, vals = load_matrix()
    assert len(cols) == len(vals[0]), "header/value width mismatch"
    idx = {c: i for i, c in enumerate(cols)}
    gidx = {g: i for i, g in enumerate(ids)}
    col_arm, col_time, col_batch = {}, {}, {}
    for r in mapping:
        col_arm[r["column"]] = arm_of(r["treatment"])
        col_time[r["column"]] = 6 if r["treatment"].startswith("6 ") else (
            24 if r["treatment"].startswith("24 ") else None)
        col_batch[r["column"]] = r["batch"]

    lc, totals = logcpm(vals, len(cols))
    sf = size_factors_mor(vals, len(cols))
    # MoR-normalised log2: log2(v/sf + 1)
    lc_mor = [[math.log2(v[i] / sf[i] + 1.0) for i in range(len(cols))] for v in vals]

    dup = {}
    for i, g in enumerate(ids):
        dup[g] = dup.get(g, 0) + 1
    dup_panel = {g: n for g, n in dup.items() if n > 1 and g in ISG_PANEL}

    def row_of(gid):
        return lc[gidx[gid]] if gid in gidx else None

    present = [g for g in ISG_PANEL if g in dup]
    missing = [f"{g}:{ISG_PANEL[g]}" for g in ISG_PANEL if g not in dup]
    present_core = [g for g in CORE12 if g in dup]

    arm_cols = {a: sorted([c for c, aa in col_arm.items() if aa == a],
                          key=lambda c: (c[0], int(c[1:]))) for a in ARMS}
    print("arm sizes:", {a: len(v) for a, v in arm_cols.items()})
    print(f"panel present {len(present)}/{len(ISG_PANEL)}; missing {missing}")
    print(f"duplicate Ensembl rows inside panel: {dup_panel or 'none'}")

    def panel_mean(tab, subset_cols, genes):
        per_gene = [st.fmean([tab[gidx[g]][idx[c]] for c in subset_cols]) for g in genes]
        return st.fmean(per_gene)

    def boot(tab, subset_cols, genes, reps=2000):
        means = []
        for _ in range(reps):
            pick = [rng.choice(subset_cols) for _ in subset_cols]
            means.append(panel_mean(tab, pick, genes))
        means.sort()
        return means[int(0.025 * reps)], means[int(0.975 * reps)]

    rows = []
    for a in ARMS:
        cs = arm_cols[a]
        if not cs:
            continue
        for block, sel in (("all", cs),
                           ("6h", [c for c in cs if col_time[c] == 6]),
                           ("24h", [c for c in cs if col_time[c] == 24])):
            if not sel:
                continue
            pm = panel_mean(lc, sel, present)
            lo, hi = boot(lc, sel, present)
            pm_mor = panel_mean(lc_mor, sel, present)
            pm_core = panel_mean(lc, sel, present_core)
            rows.append({"arm": a, "block": block, "n": len(sel), "panel_n": len(present),
                         "panel_mean_log2cpm": round(pm, 3),
                         "boot_lo": round(lo, 3), "boot_hi": round(hi, 3),
                         "panel_mean_MoR_sidecheck": round(pm_mor, 3),
                         "core12_mean_log2cpm": round(pm_core, 3),
                         "batches": ",".join(sorted({col_batch[c] for c in sel}))})

    with open(os.path.join(HERE, "arm_levels.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # per-gene table (all arms, "all" block)
    with open(os.path.join(HERE, "gene_levels.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        cols_out = ["gene", "ensembl"] + ARMS + ["log2FC_top_arm_vs_untreated"]
        fh.write("\t".join(cols_out) + "\n")
        for g in ISG_PANEL:
            v = row_of(g)
            if v is None:
                fh.write("\t".join([ISG_PANEL[g], g] + ["NA"] * (len(ARMS) + 1)) + "\n")
                continue
            means = {a: st.fmean(v[idx[c]] for c in arm_cols[a]) for a in ARMS}
            top = max((a for a in ARMS if a != "untreated"), key=lambda a: means[a])
            fh.write("\t".join([ISG_PANEL[g], g] + [f"{means[a]:.3f}" for a in ARMS]
                               + [f"{means[top] - means['untreated']:.3f}"]) + "\n")

    # signature contrast: type-II minus type-I panel mean, per arm
    sig = {}
    for a in ARMS:
        cs = arm_cols[a]
        ii = panel_mean(lc, cs, [g for g in ISG_PANEL if ISG_PANEL[g] in TYPE_II])
        i1 = panel_mean(lc, cs, [g for g in ISG_PANEL if ISG_PANEL[g] in TYPE_I])
        sig[a] = {"typeII_mean": round(ii, 3), "typeI_mean": round(i1, 3),
                  "typeII_minus_typeI": round(ii - i1, 3)}

    # bootstrap the gamma - lambda difference on the "all" block
    ga, la = arm_cols["IFN gamma"], arm_cols["IFN lambda"]
    diffs = []
    for _ in range(2000):
        gg = [rng.choice(ga) for _ in ga]
        ll = [rng.choice(la) for _ in la]
        diffs.append(panel_mean(lc, gg, present) - panel_mean(lc, ll, present))
    diffs.sort()
    gl = {"observed": round(panel_mean(lc, ga, present) - panel_mean(lc, la, present), 3),
          "boot_lo": round(diffs[50], 3), "boot_hi": round(diffs[1950], 3),
          "share_gamma_above_lambda": round(sum(d > 0 for d in diffs) / len(diffs), 4)}

    # The paper states (sentence 6) that one control batch belongs with beta/gamma and
    # the other with alpha/lambda. Compute the batch-matched induction, because pooling
    # the two control batches would mix two different sequencing dates.
    batch_blocks = {}
    for b in sorted({col_batch[c] for c in cols}):
        sel = [c for c in cols if col_batch[c] == b]
        batch_blocks[b] = {"n": len(sel),
                           "level": round(panel_mean(lc, sel, present), 3),
                           "arms": sorted({col_arm[c] for c in sel})}
    matched = {}
    for arm, control_batch in (("IFN beta", "20200722"), ("IFN gamma", "20200722"),
                               ("IFN alpha", "20200820"), ("IFN lambda", "20200820")):
        ctrl = [c for c in cols if col_batch[c] == control_batch and col_arm[c] == "untreated"]
        matched[arm] = {
            "control_batch": control_batch,
            "control_n": len(ctrl),
            "arm_level": round(panel_mean(lc, arm_cols[arm], present), 3),
            "control_level": round(panel_mean(lc, ctrl, present), 3),
            "delta_vs_batch_matched_control": round(
                panel_mean(lc, arm_cols[arm], present) - panel_mean(lc, ctrl, present), 3)}

    # direct comparison with the original implementation's published numbers
    orig_map = json.load(open(os.path.join(A2, "column_map.json"), encoding="utf-8"))
    orig_levels = orig_map["R2_panel_level"]
    name = {"IFN alpha": "A", "IFN beta": "B", "IFN gamma": "G",
            "IFN lambda": "M", "untreated": None}
    cmp_rows = []
    for a in ARMS:
        if a == "untreated":
            mine = panel_mean(lc, arm_cols[a], present)
            theirs = orig_levels.get("K", None), orig_levels.get("N", None)
            theirs_v = st.fmean([x for x in theirs if x is not None])
        else:
            mine = panel_mean(lc, arm_cols[a], present)
            theirs_v = orig_levels.get(name[a])
        cmp_rows.append({"arm": a, "independent": round(mine, 3),
                         "original": round(theirs_v, 3),
                         "abs_diff": round(abs(mine - theirs_v), 3)})

    out = {"panel_genes_present": len(present), "panel_genes_missing": missing,
           "duplicate_panel_rows": dup_panel,
           "arm_sizes": {a: len(v) for a, v in arm_cols.items()},
           "arm_levels": rows, "signature_contrast": sig,
           "gamma_minus_lambda": gl, "vs_original": cmp_rows,
           "batch_blocks": batch_blocks, "batch_matched_induction": matched}
    with open(os.path.join(HERE, "comparison.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)

    print("\narm levels (log2 CPM+1, ISG panel):")
    for r in rows:
        if r["block"] == "all":
            print(f"  {r['arm']:12s} n={r['n']:2d} level={r['panel_mean_log2cpm']:.3f} "
                  f"[{r['boot_lo']:.3f},{r['boot_hi']:.3f}] "
                  f"MoR={r['panel_mean_MoR_sidecheck']:.3f} core12={r['core12_mean_log2cpm']:.3f}")
    print("\ntime blocks:")
    for r in rows:
        if r["block"] != "all":
            print(f"  {r['arm']:12s} {r['block']:3s} n={r['n']} level={r['panel_mean_log2cpm']:.3f}")
    print("\ntype II - type I signature:", json.dumps(sig, indent=1))
    print("\ngamma - lambda (all block):", gl)
    print("\nbatch blocks:", json.dumps(batch_blocks, indent=1))
    print("batch-matched induction:", json.dumps(matched, indent=1))
    print("\nvs original R2_panel_level:")
    for r in cmp_rows:
        print(f"  {r['arm']:12s} independent={r['independent']:.3f} "
              f"original={r['original']:.3f} diff={r['abs_diff']:.3f}")
    print("\nwrote arm_levels.tsv, gene_levels.tsv, comparison.json")


if __name__ == "__main__":
    main()
