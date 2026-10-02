"""A2 filovirus functional anchor — GSE309699 (Elliff 2025, Viruses 17:1577).

Question (written before any value was inspected)
-------------------------------------------------
The manuscript's "level reading" says the host-side score tracks the *level* a cell
attains after interferon rather than the *fold induction*. GSE309699 provides the first
opportunity to test that reading against a **filovirus** functional readout in a
**normal human epithelial** cell: telomerase-immortalised normal human skin
keratinocytes (NHSK-1) pre-treated with IFN-alpha, -beta, -gamma or -lambda and then
challenged with EBOV-deltaVP30 or rVSV/EBOV-GP (Figure 3 of the paper), with the
matched transcriptome deposited here.

Before any biological number can be read, the column-to-treatment map of the deposited
count matrix must be established from the data, because the supplementary table is
labelled by author well identifiers (A1..N8) rather than by GSM. Two facts are known
independently from the series metadata:
  * group sizes are 8, 8, 5, 8 and 15 (alpha, beta, gamma, lambda, untreated);
  * IFN-gamma carries only five libraries and its column labels are G1, G2, G5, G6, G7,
    i.e. the trailing digit encodes the time block (<=4 = 6 h, >=5 = 24 h).

Pre-registered decision rule for the map
----------------------------------------
R1 (time structure): for a group that is an IFN treatment, the trailing digit must
   split an ISG panel consistently: at least 70% of panel genes must be higher in the
   digit>=5 block than in the digit<=4 block, OR at least 70% lower. A split closer to
   50/50 is what an untreated control split into two arbitrary halves looks like.
R2 (level): the untreated control is the group whose ISG panel level is lowest.
R3 (size): the candidate assignment must reproduce the known group sizes 8/8/5/8/15.
All three rules must agree; if they do not, the map is reported as unresolved and no
biological readout is taken.

Outputs (all inside this directory): column_map.json, module_levels.tsv, FINDINGS.md
companion numbers.
"""
from __future__ import annotations

import gzip
import json
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
COUNTS = os.path.join(RAW, "GSE309699_raw_counts_All_samples.csv.gz")

# ---------------------------------------------------------------- gene panel
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
RECEPTORS = {"ENSG00000142166": "IFNAR1", "ENSG00000159110": "IFNAR2",
             "ENSG00000162892": "IFNLR1", "ENSG00000159128": "IL10RB"}
HOUSEKEEPING = {"ENSG00000075624": "ACTB", "ENSG00000111640": "GAPDH"}


def load_counts():
    with gzip.open(COUNTS, "rt", encoding="utf-8-sig", errors="replace") as fh:
        header = next(fh).rstrip("\n").split(",")
        cols = header[2:]
        want = dict(ISG_PANEL)
        want.update(RECEPTORS)
        want.update(HOUSEKEEPING)
        gene_rows = {}
        n_rows = 0
        for line in fh:
            n_rows += 1
            f = line.rstrip("\n").split(",")
            if f[1] in want:
                gene_rows[want[f[1]]] = [float(x) for x in f[2:]]
    return cols, gene_rows, n_rows, len(want)


def cpm(vecs, cols):
    """Counts -> counts per million per library (the file is already scaled but not
    to a common total, so normalise explicitly)."""
    out = {g: [0.0] * len(cols) for g in vecs}
    totals = [0.0] * len(cols)
    for v in vecs.values():
        for i, x in enumerate(v):
            totals[i] += x
    # totals over the panel only would be dishonest; use the full row space instead
    return out, totals


def load_total_per_library():
    with gzip.open(COUNTS, "rt", encoding="utf-8-sig", errors="replace") as fh:
        header = next(fh).rstrip("\n").split(",")
        ncol = len(header) - 2
        totals = [0.0] * ncol
        for line in fh:
            f = line.rstrip("\n").split(",")
            try:
                vals = [float(x) for x in f[2:]]
            except ValueError:
                continue
            for i, x in enumerate(vals):
                totals[i] += x
    return totals


def group_of(col):
    return col[0]


def digit_of(col):
    return int(col[1:])


def main():
    cols, vecs, n_rows, n_want = load_counts()
    missing = [g for g, _ in [(v, k) for k, v in
                              list(ISG_PANEL.items()) + list(RECEPTORS.items())
                              + list(HOUSEKEEPING.items())] if g not in vecs]
    print(f"matrix rows {n_rows}; genes found {len(vecs)}/{n_want}; missing {missing}")
    totals = load_total_per_library()
    logcpm = {}
    for g, v in vecs.items():
        logcpm[g] = [math.log2(v[i] / totals[i] * 1e6 + 1.0) for i in range(len(cols))]

    groups = {}
    for c in cols:
        groups.setdefault(group_of(c), []).append(c)
    print("group sizes:", {g: len(v) for g, v in sorted(groups.items())})

    panel = [g for g in ISG_PANEL.values() if g in logcpm]

    # ---- R1: time structure inside candidate IFN groups (digit <=4 vs >=5)
    print("\nR1 time structure (share of panel genes higher in digit>=5 block):")
    r1 = {}
    for gname, cs in sorted(groups.items()):
        low = [c for c in cs if digit_of(c) <= 4]
        high = [c for c in cs if digit_of(c) >= 5]
        if not low or not high:
            print(f"  {gname}: no split (n={len(cs)})")
            continue
        up = down = 0
        for g in panel:
            m_hi = st.fmean(logcpm[g][cols.index(c)] for c in high)
            m_lo = st.fmean(logcpm[g][cols.index(c)] for c in low)
            if m_hi > m_lo:
                up += 1
            else:
                down += 1
        share = max(up, down) / (up + down)
        r1[gname] = {"n_low": len(low), "n_high": len(high), "up": up, "down": down,
                     "share": round(share, 3), "structured": share >= 0.70}
        print(f"  {gname}: n={len(cs)} ({len(low)} low-digit / {len(high)} high-digit) "
              f"up={up} down={down} share={share:.2f} structured={share >= 0.70}")

    # ---- R2: panel level per group
    print("\nR2 mean log2 CPM over the ISG panel, per group:")
    levels = {}
    for gname, cs in sorted(groups.items()):
        vals = [st.fmean(logcpm[g][cols.index(c)] for c in cs) for g in panel]
        levels[gname] = round(st.fmean(vals), 3)
    for gname, v in sorted(levels.items(), key=lambda kv: kv[1]):
        print(f"  {gname}: {v:.3f}")

    out = {"n_matrix_rows": n_rows, "n_genes_requested": n_want,
           "genes_found": sorted(vecs), "missing": missing,
           "group_sizes": {g: len(v) for g, v in sorted(groups.items())},
           "R1_time_structure": r1, "R2_panel_level": levels,
           "panel_genes": panel}
    with open(os.path.join(HERE, "column_map.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)

    # ---- per-condition table (group x time block)
    allgenes = panel + [g for g in list(RECEPTORS.values()) + list(HOUSEKEEPING.values())
                        if g in logcpm]
    rows = []
    for gname, cs in sorted(groups.items()):
        for block, lo, hi in (("digit<=4", 1, 4), ("digit>=5", 5, 9), ("all", 0, 99)):
            sel = [c for c in cs if lo <= digit_of(c) <= hi]
            if not sel:
                continue
            rec = {"group": gname, "block": block, "n": len(sel)}
            for g in allgenes:
                rec[g] = round(st.fmean(logcpm[g][cols.index(c)] for c in sel), 3)
            rec["panel_mean"] = round(st.fmean(rec[g] for g in panel), 3)
            rows.append(rec)
    keys = ["group", "block", "n"] + allgenes + ["panel_mean"]
    with open(os.path.join(HERE, "module_levels.tsv"), "w", encoding="utf-8", newline="") as fh:
        import csv as _csv
        w = _csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nwrote column_map.json and module_levels.tsv ({len(rows)} condition rows)")


if __name__ == "__main__":
    main()
