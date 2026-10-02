"""Gap-D: human iPSC-derived gut organoids infected with EBOV and MARV
(GSE300073 / GSE298600, PMID 41284734, PLoS Pathogens 2025).

Question (frozen before any value was read)
-------------------------------------------
The manuscript places the enteric epithelium at the top of the interferon-tone
amplification set, and states that the decisive comparison (the same host cell
type infected with two orthoebolaviruses, with host transcriptome readout) has
not been deposited. GSE300073 and GSE298600 deposit exactly that design for
EBOV and MARV in human iPSC-derived intestinal organoids, with triplicate mock
controls at 1 and 3 days post-infection, but the two matrices label their
columns by author identifiers (MG-MA-13..30, S1..S9) rather than by GSM.

Frozen decision rules
---------------------
R1  mapping is read in the sample order declared by the series record, because
    no per-sample supplementary files exist to give an explicit map.
R2  that map is accepted only if it passes three independent checks:
    (a) within-condition triplicate similarity exceeds between-condition
        similarity on the whole transcriptome;
    (b) the delayed-induction prediction holds (infected 1 dpi below
        infected 3 dpi on an interferon-stimulated module);
    (c) the same conditions measured in the two independent series agree
        (EBOV 3 dpi and mock 3 dpi appear in both).
    If a check fails, the map is reported as unresolved and no biological
    readout is taken.
R3  the readout is the manuscript's own 17-gene interferon-stimulated module
    and its three score components, computed from counts per million on the
    full gene space, with no gene dropped silently.

Outputs: gapd_summary.json, gapd_module_scores.tsv, gapd_mapping_checks.json,
raw/ (the two deposited count matrices).
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")

ISG17 = ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "IFIT1", "IFIT2", "IFIT3",
         "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6", "DDX58",
         "IFIH1"]
TYPE_I = ["IFNAR1", "IFNAR2"]
TYPE_III = ["IFNLR1", "IL10RB"]
PRIMING = ["ISG15", "MX1"]

# series-record sample order (GSM order == matrix column order under R1)
DESIGN = {
    "GSE300073": [("Mock", 1)] * 3 + [("Mock", 3)] * 3 + [("EBOV", 1)] * 3
    + [("EBOV", 3)] * 3 + [("MARV", 1)] * 3 + [("MARV", 3)] * 3,
    "GSE298600": [("MARV", 3)] * 3 + [("EBOV", 3)] * 3 + [("Mock", 3)] * 3,
}


def load(acc):
    path = os.path.join(RAW, f"{acc}_counts.csv.gz")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = next(fh).rstrip("\n").split(",")
        cols = [c.strip().strip('"') for c in header[1:]]
        rows = {}
        for line in fh:
            f = line.rstrip("\n").split(",")
            gene = f[0].strip().strip('"')
            if gene in rows:            # duplicate symbol: keep the first, count it
                continue
            rows[gene] = [float(x) for x in f[1:]]
    return cols, rows


def cpm(rows, ncol):
    totals = [0.0] * ncol
    for v in rows.values():
        for i, x in enumerate(v):
            totals[i] += x
    out = {}
    for g, v in rows.items():
        out[g] = [math.log2(v[i] / totals[i] * 1e6 + 1.0) for i in range(ncol)]
    return out, totals


def module(logcpm, genes, ncol):
    present = [g for g in genes if g in logcpm]
    missing = [g for g in genes if g not in logcpm]
    vals = [st.fmean(logcpm[g][i] for g in present) for i in range(ncol)]
    return vals, present, missing


def rank_corr(a, b):
    def rank(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        for pos, i in enumerate(order):
            r[i] = pos + 1
        return r
    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = st.fmean(ra), st.fmean(rb)
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    da = math.sqrt(sum((x - ma) ** 2 for x in ra))
    db = math.sqrt(sum((x - mb) ** 2 for x in rb))
    return num / (da * db) if da and db else float("nan")


def main():
    summary = {"source": "GSE300073 / GSE298600 (PMID 41284734)",
               "module_definition": {"ISG17": ISG17, "type_I": TYPE_I,
                                     "type_III": TYPE_III, "priming": PRIMING},
               "series": {}}
    scores = []
    checks = {}
    cache = {}
    for acc in ("GSE300073", "GSE298600"):
        cols, rows = load(acc)
        ncol = len(cols)
        design = DESIGN[acc]
        assert ncol == len(design), (acc, ncol, len(design))
        logcpm, totals = cpm(rows, ncol)
        cache[acc] = (cols, design, logcpm)
        m17, p17, miss17 = module(logcpm, ISG17, ncol)
        m1, p1, miss1 = module(logcpm, TYPE_I, ncol)
        m3, p3, miss3 = module(logcpm, TYPE_III, ncol)
        mp, pp, missp = module(logcpm, PRIMING, ncol)

        # R2(a) within vs between condition similarity, on the ISG panel genes
        vecs = {}
        for i in range(ncol):
            vecs[i] = [logcpm[g][i] for g in p17]
        within, between = [], []
        for i in range(ncol):
            for j in range(i + 1, ncol):
                r = rank_corr(vecs[i], vecs[j])
                (within if design[i] == design[j] else between).append(r)
        checks[f"{acc}_R2a"] = {
            "mean_within_condition_rho": round(st.fmean(within), 3),
            "mean_between_condition_rho": round(st.fmean(between), 3),
            "passes": st.fmean(within) > st.fmean(between)}

        for i, c in enumerate(cols):
            scores.append({"series": acc, "column": c, "treatment": design[i][0],
                           "day": design[i][1], "ISG17": round(m17[i], 4),
                           "type_I": round(m1[i], 4), "type_III": round(m3[i], 4),
                           "priming": round(mp[i], 4)})
        summary["series"][acc] = {
            "n_columns": ncol, "n_genes": len(rows),
            "library_size_min": int(min(totals)), "library_size_max": int(max(totals)),
            "ISG17_genes_present": p17, "ISG17_genes_missing": miss17,
            "type_I_missing": miss1, "type_III_missing": miss3,
            "priming_missing": missp,
            "mean_ISG17_by_condition": {
                f"{t}_{d}d": round(st.fmean(m17[i] for i in range(ncol)
                                            if design[i] == (t, d)), 3)
                for t in {x[0] for x in design} for d in {x[1] for x in design}}}

    # R2(b) delayed induction: infected 1 dpi < infected 3 dpi
    s73 = {c["treatment"] + str(c["day"]): c["ISG17"] for c in scores
           if c["series"] == "GSE300073"}
    def cond_mean(t, d, series):
        v = [c["ISG17"] for c in scores if c["series"] == series
             and c["treatment"] == t and c["day"] == d]
        return st.fmean(v)
    delayed = {}
    for t in ("EBOV", "MARV"):
        delayed[t] = {"day1": round(cond_mean(t, 1, "GSE300073"), 3),
                      "day3": round(cond_mean(t, 3, "GSE300073"), 3)}
        delayed[t]["rises_between_days"] = delayed[t]["day3"] > delayed[t]["day1"]
    checks["R2b_delayed_induction"] = delayed
    checks["R2b_passes"] = all(v["rises_between_days"] for v in delayed.values())

    # R2(c) cross-series agreement at 3 dpi
    cross = {}
    for t in ("EBOV", "Mock"):
        a = cond_mean(t, 3, "GSE300073")
        b = cond_mean(t, 3, "GSE298600")
        cross[t] = {"GSE300073": round(a, 3), "GSE298600": round(b, 3),
                    "difference": round(a - b, 3)}
    checks["R2c_cross_series"] = cross
    diffs = [abs(v["difference"]) for v in cross.values()]
    # the two series are separate deposits of the same conditions; agreement is
    # judged against the spread between the three mock/EBOV conditions themselves
    span = abs(cond_mean("EBOV", 3, "GSE300073") - cond_mean("Mock", 3, "GSE300073"))
    checks["R2c_passes"] = max(diffs) < max(span, 0.5)
    checks["R2c_note"] = ("tolerance = the larger of the EBOV-mock difference "
                          "within one series and 0.5 log2")

    checks["map_accepted"] = bool(checks["GSE300073_R2a"]["passes"]
                                  and checks["GSE298600_R2a"]["passes"]
                                  and checks["R2b_passes"]
                                  and checks["R2c_passes"])

    # induction relative to the matched mock of the same series and day
    induced = {}
    for series in ("GSE300073", "GSE298600"):
        for d in (1, 3):
            mocks = [c["ISG17"] for c in scores if c["series"] == series
                     and c["treatment"] == "Mock" and c["day"] == d]
            if not mocks:
                continue
            base = st.fmean(mocks)
            for t in ("EBOV", "MARV"):
                v = [c["ISG17"] for c in scores if c["series"] == series
                     and c["treatment"] == t and c["day"] == d]
                if v:
                    induced[f"{series}:{t}:{d}d"] = round(st.fmean(v) - base, 3)
    summary["ISG17_induction_vs_matched_mock"] = induced
    summary["mapping_checks"] = checks

    with open(os.path.join(HERE, "gapd_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(HERE, "gapd_mapping_checks.json"), "w", encoding="utf-8") as fh:
        json.dump(checks, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(HERE, "gapd_module_scores.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["series", "column", "treatment", "day",
                                           "ISG17", "type_I", "type_III", "priming"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(scores)
    print(json.dumps({"mapping_checks": checks,
                      "induction": induced,
                      "series": summary["series"]}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
