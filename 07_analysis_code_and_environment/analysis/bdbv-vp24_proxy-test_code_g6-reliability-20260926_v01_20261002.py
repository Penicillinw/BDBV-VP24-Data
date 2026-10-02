"""G6 step 5 (v2): per-gene / manuscript-component tests, platform completeness,
and the stratified re-analysis that decides how GSE21158 may be read.

Two defects in v1 are fixed here and both are recorded in FINDINGS.md:

1. v1 matched genes on the GPL6098 `Symbol` column only. Four of the fourteen ISG
   module genes are not present under that symbol on this platform (ISG15, IFIT3,
   RSAD2, IFI6) and IFNLR1/GAPDH are absent outright. `module_mean` returned NaN for
   them and `spearman` then ranked NaN vectors into a spurious rho = +1.000, which
   v1 reported as "flat across most lines" - it was actually "not measured".
   v2 resolves genes by Symbol **or** by the platform's own `Synonym` column
   (GPL6098 lists ISG15 as `G1P2`, synonym `UCRP;IFI15;ISG15`), then states which
   genes remain absent from the array instead of reporting a number for them.

2. v1 reported only the pooled correlation. The induction phenotype is bimodal
   (five lines respond to IFN-alpha-2a, five do not), so the pooled rho is a
   stratification artefact. v2 reports the within-stratum values alongside it.

The residual "noise" question (regression to the mean) is kept, but only as a
footnote on the pooled statistic, which is no longer the headline.
"""

from __future__ import annotations

import gzip
import json
import math
import os
from collections import defaultdict

ROOT = r"G:\本迪布焦研究"
RAW = os.path.join(ROOT, "analysis", "g6_ifn_responsiveness_20260926", "raw")
OUT = os.path.join(ROOT, "analysis", "g6_ifn_responsiveness_20260926")

ISG = ["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3",
       "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6"]
RECEPTORS = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB"]
HK = ["ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "HPRT1"]
GENES = list(dict.fromkeys(ISG + RECEPTORS + HK))


def load():
    header = None
    symbol_of: dict[str, str] = {}
    rows: list[tuple[str, str, list[str]]] = []
    for line in open(os.path.join(RAW, "GPL6098_platform_table.txt"),
                     encoding="utf-8", errors="replace"):
        if line.startswith("ID\t"):
            header = line.rstrip("\n").split("\t")
            continue
        if header and line.strip() and not line.startswith("!"):
            p = line.rstrip("\n").split("\t")
            p += [""] * (len(header) - len(p))
            pid = p[0].strip()
            sym = p[7].strip().split("///")[0].strip()
            syn = [s.strip().upper() for s in p[13].split(";") if s.strip()]
            rows.append((pid, sym, syn))
            if sym and sym not in {"---", ""}:
                symbol_of[pid] = sym

    # resolve each requested gene by its own symbol first, then by the platform's
    # Synonym column (this is how ISG15 is recoverable: GPL6098 symbols it G1P2)
    #
    # A third path is required for genes that appear only under a legacy symbol whose
    # Synonym column does not carry the current name at all. GAPDH is the case here:
    # probes 1940184 and 1850181 are symbol GAPD (NM_002046.2, "glyceraldehyde-3-
    # phosphate dehydrogenase (GAPD)") with an EMPTY Synonym column, so symbol and
    # synonym search both miss it. Without this table the housekeeping module silently
    # drops GAPDH and the run reports "GAPDH absent from array", which is false.
    LEGACY_ALIASES = {"GAPDH": ("GAPD",), "ISG15": ("G1P2",), "IFNLR1": ("IL28RA",)}
    resolved: dict[str, list[str]] = {}
    for gene in GENES:
        up = gene.upper()
        hits = [pid for pid, sym, _ in rows if sym.upper() == up]
        if not hits:
            hits = [pid for pid, _, syn in rows if up in syn]
        if not hits:
            aliases = {a.upper() for a in LEGACY_ALIASES.get(gene, ())}
            hits = [pid for pid, sym, _ in rows if sym.upper() in aliases]
        if hits:
            resolved[gene] = hits

    probe2gene: dict[str, str] = {}
    for gene, probes in resolved.items():
        for pid in probes:
            probe2gene.setdefault(pid, gene)

    absent = [g for g in GENES if g not in resolved]

    titles, vecs = [], {}
    with gzip.open(os.path.join(RAW, "GSE21158_series_matrix.txt.gz"), "rt",
                   encoding="utf-8", errors="replace") as fh:
        in_t = False
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_t = True
                next(fh)
            elif in_t:
                if line.startswith("!series_matrix_table_end"):
                    break
                p = line.rstrip("\n").split("\t")
                if len(p) > 1:
                    try:
                        vecs[p[0].strip('"')] = [float(x) for x in p[1:]]
                    except ValueError:
                        pass
    return probe2gene, resolved, absent, titles, vecs


def design(titles):
    out = []
    for t in titles:
        low = t.lower()
        line = low.split()[0]
        cond = ("IFN24" if "24 hr interferon" in low
                else "IFN4" if "4 hr interferon" in low else "control")
        out.append((line, cond))
    return out


def spearman(a, b):
    def rank(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    da = math.sqrt(sum((v - ma) ** 2 for v in ra))
    db = math.sqrt(sum((v - mb) ** 2 for v in rb))
    return num / (da * db) if da and db else float("nan")


def main():
    probe2gene, resolved, absent, titles, vecs = load()
    des = design(titles)
    lines = sorted({d[0] for d in des})
    per: dict[str, dict[str, dict[str, list[float]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list)))
    for probe, vec in vecs.items():
        sym = probe2gene.get(probe)
        if not sym:
            continue
        for i, (line, cond) in enumerate(des):
            per[sym][line][cond].append(vec[i])

    print(f"genes resolved on GPL6098: {len(resolved)}/{len(GENES)}")
    for g in GENES:
        print(f"  {g:9s} " + (",".join(resolved[g]) if g in resolved else "ABSENT from array"))
    if absent:
        print(f"  -> absent, no number reported: {', '.join(absent)}")

    def module_mean(sym_list, line, cond, reps=None):
        vals = []
        for g in sym_list:
            v = per.get(g, {}).get(line, {}).get(cond, [])
            if not v:
                continue
            vals.append(sum(v) / len(v) if reps is None else sum(v[i] for i in reps) / len(reps))
        return sum(vals) / len(vals) if vals else float("nan")

    # reliability of the baseline score across the 3 control replicates
    isg_present = [g for g in ISG if g in resolved]
    rec_present = [g for g in RECEPTORS if g in resolved]
    hk_present = [g for g in HK if g in resolved]

    rep_scores = {r: [] for r in range(3)}
    for line in lines:
        for r in range(3):
            vals = []
            for g in isg_present:
                v = per.get(g, {}).get(line, {}).get("control", [])
                if len(v) > r:
                    vals.append(v[r])
            rep_scores[r].append(sum(vals) / len(vals) if vals else float("nan"))

    grand = sum(sum(v) for v in rep_scores.values()) / (3 * len(lines))
    ms_between = 3 * sum((sum(rep_scores[r][i] for r in range(3)) / 3 - grand) ** 2
                         for i in range(len(lines))) / (len(lines) - 1)
    ms_within = sum(sum((rep_scores[r][i] - sum(rep_scores[q][i] for q in range(3)) / 3) ** 2
                        for r in range(3)) for i in range(len(lines))) / (len(lines) * 2)
    reliability = (ms_between - ms_within) / (ms_between + 2 * ms_within) if ms_between > 0 else 0.0
    reliability = max(0.0, min(1.0, reliability))

    base = [module_mean(isg_present, l, "control") for l in lines]
    ind24 = [module_mean(isg_present, l, "IFN24") - module_mean(isg_present, l, "control") for l in lines]
    ind4 = [module_mean(isg_present, l, "IFN4") - module_mean(isg_present, l, "control") for l in lines]
    rec = [module_mean(rec_present, l, "control") for l in lines]
    hk = [module_mean(hk_present, l, "control") for l in lines]

    # the manuscript's exact component definitions
    base_2gene = [module_mean(["ISG15", "MX1"], l, "control") for l in lines]
    ind_2gene = [
        module_mean(["ISG15", "MX1"], l, "IFN24") - module_mean(["ISG15", "MX1"], l, "control")
        for l in lines
    ]
    base_typeI = [module_mean(["IFNAR1", "IFNAR2"], l, "control") for l in lines]
    base_typeIII = [module_mean(["IFNLR1", "IL10RB"], l, "control") for l in lines]
    delta_like = [(a + b + c) / 3 for a, b, c in zip(base_typeI, base_typeIII, base_2gene)]

    # ---- the decisive split: responders vs non-responders -------------------
    order = sorted(range(len(lines)), key=lambda i: ind24[i])
    responsive = [i for i in order if ind24[i] >= 0.3]
    unresponsive = [i for i in order if ind24[i] < 0.3]

    def sub(vals, idx):
        return [vals[i] for i in idx]

    def mean_of(vals, idx):
        xs = sub(vals, idx)
        return round(sum(xs) / len(xs), 3) if xs else float("nan")

    def rho_sub(base_vals, target_vals, idx):
        if len(idx) < 3:
            return None
        return round(spearman(sub(base_vals, idx), sub(target_vals, idx)), 3)

    stratified = {
        "rule": "responder = ISG module induction >= 0.3 log units at 24 h; the "
                "threshold is fixed on the response axis only, never on baseline",
        "pooled_spearman_baseline_vs_induction": round(spearman(base, ind24), 3),
        "responders": {
            "n": len(responsive),
            "lines": [lines[i] for i in responsive],
            "mean_baseline": mean_of(base, responsive),
            "mean_induction": mean_of(ind24, responsive),
            "spearman_baseline_vs_induction": rho_sub(base, ind24, responsive),
            "spearman_2gene_priming_vs_own_induction": rho_sub(base_2gene, ind_2gene, responsive),
            "spearman_typeI_vs_ISG_induction": rho_sub(base_typeI, ind24, responsive),
            "spearman_delta_like_vs_induction": rho_sub(delta_like, ind24, responsive),
        },
        "non_responders": {
            "n": len(unresponsive),
            "lines": [lines[i] for i in unresponsive],
            "mean_baseline": mean_of(base, unresponsive),
            "mean_induction": mean_of(ind24, unresponsive),
            "spearman_baseline_vs_induction": rho_sub(base, ind24, unresponsive),
        },
    }

    rho_obs = spearman(base, ind24)
    rho_dis = rho_obs / math.sqrt(reliability) if reliability > 0 else float("nan")
    out = {
        "n_cell_lines": len(lines),
        "platform_completeness": {
            "genes_requested": GENES,
            "genes_resolved": {g: resolved[g] for g in resolved},
            "genes_absent_from_array": absent,
            "ISG_module_used": isg_present,
            "receptor_module_used": rec_present,
            "housekeeping_module_used": hk_present,
            "note": "GPL6098 is an old array and three requested genes are reachable only "
                    "through a legacy name: ISG15 (probe 4590398, symbol G1P2, Synonym "
                    "UCRP;IFI15;ISG15) and IFNLR1 (probe 6520180, symbol IL28RA, Synonym "
                    "IFNLR;LICR2;IFNLR1;CRF2/12) resolve through the Synonym column, and "
                    "GAPDH (probes 1940184, 1850181, symbol GAPD, Synonym empty) resolves "
                    "only through the LEGACY_ALIASES table. Genes that remain unmatched "
                    "(IFIT3, RSAD2, IFI6) have no probe on this array and are reported as "
                    "absent rather than as a number.",
        },
        "manuscript_components_as_measurable_here": {
            "two_gene_priming_term": "ISG15 (probe 4590398, symbol G1P2) + MX1",
            "type_I_term": "IFNAR1 + IFNAR2 (both present)",
            "type_III_term": "IFNLR1 (probe 6520180, symbol IL28RA, Synonym "
                             "IFNLR;LICR2;IFNLR1;CRF2/12) + IL10RB - both present",
            "composite": "mean of the three terms above, i.e. six genes in total",
        },
        "stratified": stratified,
        "baseline_reliability_ICC": round(reliability, 3),
        "spearman_baseline_vs_induction24": round(rho_obs, 3),
        "disattenuated_spearman": round(rho_dis, 3),
        "spearman_baseline_vs_induction4h": round(spearman(base, ind4), 3),
        "spearman_receptor_vs_induction24": round(spearman(rec, ind24), 3),
        "spearman_housekeeping_vs_induction24": round(spearman(hk, ind24), 3),
        "manuscript_component_tests": {
            "baseline_ISG15_MX1_vs_own_induction": round(spearman(base_2gene, ind_2gene), 3),
            "baseline_IFNAR1_IFNAR2_vs_ISG_induction": round(spearman(base_typeI, ind24), 3),
            "baseline_IFNLR1_IL10RB_vs_ISG_induction": round(spearman(base_typeIII, ind24), 3),
            "delta_like_baseline_vs_ISG_induction": round(spearman(delta_like, ind24), 3),
        },
        "delta_like_values": [round(x, 3) for x in delta_like],
        "per_gene_spearman_baseline_vs_own_induction": {},
        "cell_lines": lines,
        "baseline_values": [round(x, 3) for x in base],
        "induction24_values": [round(x, 3) for x in ind24],
    }
    for g in ISG:
        if g not in resolved:
            out["per_gene_spearman_baseline_vs_own_induction"][g] = "not on array"
            continue
        b = [module_mean([g], l, "control") for l in lines]
        i24 = [module_mean([g], l, "IFN24") - module_mean([g], l, "control") for l in lines]
        try:
            out["per_gene_spearman_baseline_vs_own_induction"][g] = round(spearman(b, i24), 3)
        except Exception:  # noqa: BLE001
            continue

    with open(os.path.join(OUT, "g6_reliability.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
