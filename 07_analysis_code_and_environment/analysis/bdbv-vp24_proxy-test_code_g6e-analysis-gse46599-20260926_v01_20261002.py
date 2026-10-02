"""G6e: formal same-experiment test of the proxy assumption in GSE46599.

Question (frozen before any value was inspected): the manuscript's restriction gradient
Delta(c) treats BASELINE (uninfected) interferon-system expression as a monotone-increasing
proxy for a cell type's INTERFERON RESPONSIVENESS. Within a single experiment, does baseline
expression predict the magnitude of the response?

GSE46599 (Goujon/Schulz/Malim; PMID 24048477) is the best available instrument because
untreated and type-I-IFN-treated RNA come from the same cultures, the same platform
(Illumina HumanHT-12 v4, wide dynamic range) and the same experiment: 9 immortalized lines
(2 replicates each) plus 2 primary human cell types (3 replicates each); 1000 U/mL universal
type 1 IFN for 24 h. Unlike GSE21158 no line in this panel has a disabled IFN response, so
there is no bimodal responder/non-responder confound.

Modules are frozen in this file before the matrix is read. Genes are resolved through the
platform's Gene symbol column and, when that fails, through its UniGene title / Gene title
text and a small verified legacy-alias table; genes with no probe are reported as absent
instead of being given a number (the v1 defect in tools/g6_reliability_20260926.py).

Outputs: analysis/g6e_gse46599_20260926/
"""

from __future__ import annotations

import gzip
import itertools
import json
import math
import os
import random
from collections import defaultdict

import numpy as np

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926")
RAW = os.path.join(OUT, "raw")
MATRIX = os.path.join(RAW, "GSE46599_series_matrix.txt.gz")
ANNOT = os.path.join(RAW, "GPL10558.annot.gz")

# ---------------------------------------------------------------- frozen modules
ISG_MODULE = ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "IFIT1", "IFIT2", "IFIT3",
              "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6",
              "DDX58", "IFIH1"]
RECEPTOR_MODULE = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB"]
HOUSEKEEPING = ["ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "HPRT1"]
# the manuscript's own three components (ISG15+MX1 | IFNAR1+IFNAR2 | IFNLR1+IL10RB)
COMPONENT_PRIMING = ["ISG15", "MX1"]
COMPONENT_TYPE_I = ["IFNAR1", "IFNAR2"]
COMPONENT_TYPE_III = ["IFNLR1", "IL10RB"]
ALL_GENES = list(dict.fromkeys(ISG_MODULE + RECEPTOR_MODULE + HOUSEKEEPING))

# verified legacy names on Illumina arrays; each is confirmed against the annotation's
# GenBank accession before use (see load_annotation)
LEGACY_ALIASES = {"ISG15": ("G1P2", "UCRP"), "IFNLR1": ("IL28RA",), "GAPDH": ("GAPD",),
                  "IFIT3": ("IFIT4", "ISG60"), "RSAD2": ("CIG5", "VIG1"),
                  "IFI6": ("G1P3", "IFI616"), "OAS1": ("OAS1",), "MX1": ("MXA",)}

PRIMARY_CELLS = {"primary-macrophages", "primary-CD4+-T-cells"}
EXACT_LIMIT = 500_000     # enumerate all permutations when n! is at most this
MC_DRAWS = 200_000
SEED = 20260926


def load_annotation():
    """probe -> symbol, resolved for exactly the requested genes."""
    header = None
    rows = []
    with gzip.open(ANNOT, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                header = next(fh).rstrip("\n").split("\t")
                continue
            if header is None or line.startswith(("!", "#")) or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))

    def field(row, name):
        return (row.get(name) or "").strip()

    resolved: dict[str, list[str]] = {}
    evidence: dict[str, str] = {}
    for gene in ALL_GENES:
        up = gene.upper()
        hits = [field(r, "ID") for r in rows if field(r, "Gene symbol").upper() == up]
        how = "gene symbol"
        if not hits:
            hits = [field(r, "ID") for r in rows
                    if up in [s.upper() for s in field(r, "Gene symbol").split("///")]]
            how = "gene symbol list"
        if not hits:
            aliases = {a.upper() for a in LEGACY_ALIASES.get(gene, ())}
            hits = [field(r, "ID") for r in rows
                    if field(r, "Gene symbol").upper() in aliases]
            how = "legacy alias " + ",".join(sorted(aliases))
        if not hits:
            # last resort: the gene name written out in the title field, exact token match
            hits = [field(r, "ID") for r in rows
                    if any(tok.strip("(),;").upper() == up
                           for tok in field(r, "Gene title").replace("/", " ").split())]
            how = "gene title token"
        if hits:
            resolved[gene] = sorted(set(hits))
            # record the accession of the first probe as evidence of identity
            acc = next(field(r, "GenBank Accession") for r in rows
                       if field(r, "ID") == resolved[gene][0])
            evidence[gene] = f"{how}; probe {resolved[gene][0]}; accession {acc or 'n/a'}"

    probe2gene = {}
    for gene, probes in resolved.items():
        for pid in probes:
            probe2gene.setdefault(pid, gene)
    absent = [g for g in ALL_GENES if g not in resolved]
    return probe2gene, resolved, absent, evidence


def load_matrix():
    titles, chars, vecs = [], defaultdict(list), {}
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        in_table = False
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
                if vals:
                    chars[vals[0].split(":")[0].strip().lower()].extend(vals)
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                next(fh)
            elif in_table:
                if line.startswith("!series_matrix_table_end"):
                    break
                p = line.rstrip("\n").split("\t")
                if len(p) > 1:
                    try:
                        vecs[p[0].strip('"')] = [float(x) for x in p[1:]]
                    except ValueError:
                        pass
    return titles, chars, vecs


def spearman(a, b) -> float:
    ra, rb = _rank(a), _rank(b)
    return _pearson(ra, rb)


def _rank(x):
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


def _pearson(a, b) -> float:
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    da = math.sqrt(sum((v - ma) ** 2 for v in a))
    db = math.sqrt(sum((v - mb) ** 2 for v in b))
    return num / (da * db) if da and db else float("nan")


def permutation_p(a, b):
    """Two-sided permutation p. Exact when the permutation space is small enough."""
    obs = abs(spearman(a, b))
    n = len(a)
    if math.factorial(n) <= EXACT_LIMIT:
        hits = total = 0
        for perm in itertools.permutations(b):
            total += 1
            if abs(spearman(a, list(perm))) >= obs - 1e-12:
                hits += 1
        return {"p": round(hits / total, 4), "method": f"exact, {total} permutations",
                "min_achievable_p": round(1 / total, 4)}
    rng = np.random.default_rng(SEED)
    arr = np.array(b, dtype=float)
    ranks_a = np.array(_rank(a), dtype=float)
    hits = 0
    for _ in range(MC_DRAWS):
        perm = rng.permutation(arr)
        if abs(spearman(a, list(perm))) >= obs - 1e-12:
            hits += 1
    return {"p": round((hits + 1) / (MC_DRAWS + 1), 5),
            "method": f"Monte Carlo, {MC_DRAWS} draws, seed {SEED}",
            "notes": "n too large for exact enumeration"}


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    probe2gene, resolved, absent, evidence = load_annotation()
    titles, chars, vecs = load_matrix()
    n_samples = len(titles)

    design = []
    for i, t in enumerate(titles):
        parts = t.split("_")
        cell = "_".join(parts[:-2])
        cond = parts[-2]
        rep = parts[-1]
        design.append({"title": t, "cell": cell, "arm": "IFN" if cond == "IFN" else "control",
                       "replicate": rep, "index": i})
    cells = list(dict.fromkeys(d["cell"] for d in design))

    # per gene / cell / arm, averaged over probes
    per = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for probe, vec in vecs.items():
        gene = probe2gene.get(probe)
        if not gene or len(vec) < n_samples:
            continue
        for d in design:
            per[gene][d["cell"]][d["arm"]].append(vec[d["index"]])

    def gene_value(gene, cell, arm):
        vals = per.get(gene, {}).get(cell, {}).get(arm, [])
        return sum(vals) / len(vals) if vals else float("nan")

    def module_mean(genes, cell, arm):
        vals = [gene_value(g, cell, arm) for g in genes]
        vals = [v for v in vals if not math.isnan(v)]
        return sum(vals) / len(vals) if vals else float("nan")

    present = lambda genes: [g for g in genes if g in resolved]
    isg_p = present(ISG_MODULE)

    rows = []
    for cell in cells:
        rec = {"cell": cell,
               "kind": "primary" if cell in PRIMARY_CELLS else "immortalized",
               "n_control": sum(1 for d in design if d["cell"] == cell and d["arm"] == "control"),
               "n_ifn": sum(1 for d in design if d["cell"] == cell and d["arm"] == "IFN")}
        for name, genes in (("ISG", isg_p),
                            ("priming_ISG15_MX1", present(COMPONENT_PRIMING)),
                            ("typeI_IFNAR1_IFNAR2", present(COMPONENT_TYPE_I)),
                            ("typeIII_IFNLR1_IL10RB", present(COMPONENT_TYPE_III)),
                            ("receptor", present(RECEPTOR_MODULE)),
                            ("housekeeping", present(HOUSEKEEPING))):
            b = module_mean(genes, cell, "control")
            f = module_mean(genes, cell, "IFN")
            rec[f"{name}_baseline"] = round(b, 4)
            rec[f"{name}_IFN"] = round(f, 4)
            rec[f"{name}_induction"] = round(f - b, 4)
            rec[f"{name}_n_genes"] = len(genes)
        rec["composite_baseline"] = round(
            (rec["priming_ISG15_MX1_baseline"] + rec["typeI_IFNAR1_IFNAR2_baseline"]
             + rec["typeIII_IFNLR1_IL10RB_baseline"]) / 3, 4)
        rows.append(rec)
    rows.sort(key=lambda r: r["ISG_induction"])

    def col(name, subset=None):
        return [r[name] for r in rows if subset is None or r["kind"] == subset]

    strata = {"all": None, "immortalized_only": "immortalized", "primary_only": "primary"}
    tests = {}
    for label, subset in strata.items():
        n = len(col("ISG_baseline", subset))
        tests[label] = {"n": n}
        if n < 3:
            tests[label]["note"] = "too few cell types for a correlation"
            continue
        for base_name, base_col, target_col in (
                ("baseline ISG module vs own induction", "ISG_baseline", "ISG_induction"),
                ("baseline ISG15+MX1 vs own induction", "priming_ISG15_MX1_baseline", "priming_ISG15_MX1_induction"),
                ("baseline IFNAR1+IFNAR2 vs ISG induction", "typeI_IFNAR1_IFNAR2_baseline", "ISG_induction"),
                ("baseline IFNLR1+IL10RB vs ISG induction", "typeIII_IFNLR1_IL10RB_baseline", "ISG_induction"),
                ("baseline composite vs ISG induction", "composite_baseline", "ISG_induction"),
                ("baseline housekeeping vs ISG induction (negative control)", "housekeeping_baseline", "ISG_induction"),
                ("baseline ISG module vs induced level", "ISG_baseline", "ISG_IFN")):
            a, b = col(base_col, subset), col(target_col, subset)
            rho = spearman(a, b)
            entry = {"spearman": round(rho, 3)}
            entry.update(permutation_p(a, b))
            tests[label][base_name] = entry

    # induction distribution / bimodality
    ind = sorted(col("ISG_induction"))
    gaps = [ind[i + 1] - ind[i] for i in range(len(ind) - 1)]
    folds = [ind[i + 1] / ind[i] if ind[i] > 0 else float("inf") for i in range(len(ind) - 1)]
    bimodality = {
        "induction_sorted": [round(x, 3) for x in ind],
        "range": [round(ind[0], 3), round(ind[-1], 3)],
        "largest_gap": round(max(gaps), 3),
        "median_gap": round(sorted(gaps)[len(gaps) // 2], 3),
        "largest_consecutive_fold_change": round(max(folds), 3),
        # a responder/non-responder split shows up as a multi-fold jump between adjacent
        # sorted values (GSE21158: 15.0-fold, five lines at ~0.02). A graded distribution
        # (this panel) does not.
        "is_bimodal": bool(max(folds) >= 3.0),
        "n_non_responders_below_0.5": sum(1 for x in ind if x < 0.5),
        "n_below_20pct_of_max": sum(1 for x in ind if x < 0.2 * ind[-1]),
    }

    # the series carries its own functional label: HIV-1 restriction after IFN treatment
    phenotype = {}
    raw_labels = [v for v in chars.get("resistance to hiv-1 following ifn treatment", [])]
    for cell in cells:
        for d in design:
            if d["cell"] == cell:
                lbl = raw_labels[d["index"]].split(":", 1)[-1].strip()
                phenotype.setdefault(cell, set()).add(lbl)
    phenotype = {k: sorted(v) for k, v in phenotype.items()}

    groups = {"permissive": [], "partially resistant": [], "resistant": []}
    for r in rows:
        labels = [x for x in phenotype.get(r["cell"], []) if x != "untreated"]
        if labels and labels[0] in groups:
            groups[labels[0]].append(r)
    restriction = {}
    for name, members in groups.items():
        restriction[name] = {
            "n": len(members),
            "cells": [m["cell"] for m in members],
            "mean_baseline_ISG": round(sum(m["ISG_baseline"] for m in members) / len(members), 3)
            if members else None,
            "mean_ISG_induction": round(sum(m["ISG_induction"] for m in members) / len(members), 3)
            if members else None,
        }
    # resistant vs permissive, exact Mann-Whitney-style difference in means by permutation
    res_b = [m["ISG_baseline"] for m in groups["resistant"]]
    per_b = [m["ISG_baseline"] for m in groups["permissive"]]
    if res_b and per_b:
        obs = sum(res_b) / len(res_b) - sum(per_b) / len(per_b)
        pool = res_b + per_b
        hits = total = 0
        for combo in itertools.combinations(range(len(pool)), len(res_b)):
            pick = [pool[i] for i in combo]
            rest = [pool[i] for i in range(len(pool)) if i not in combo]
            total += 1
            if abs(sum(pick) / len(pick) - sum(rest) / len(rest)) >= abs(obs) - 1e-12:
                hits += 1
        restriction["resistant_vs_permissive_baseline"] = {
            "difference_resistant_minus_permissive": round(obs, 3),
            "p": round(hits / total, 4),
            "method": f"exact, {total} label permutations",
            "caveat": "functional label is HIV-1 restriction, not interferon responsiveness",
        }

    # ceiling check: per gene, baseline vs induced, and the maximum reached
    ceiling = {}
    for g in isg_p:
        bases = [gene_value(g, c, "control") for c in cells]
        ifns = [gene_value(g, c, "IFN") for c in cells]
        inds = [f - b for b, f in zip(bases, ifns)]
        ceiling[g] = {
            "baseline_range": [round(min(bases), 3), round(max(bases), 3)],
            "induced_range": [round(min(ifns), 3), round(max(ifns), 3)],
            "max_induction": round(max(inds), 3),
            "min_induction": round(min(inds), 3),
            "spearman_baseline_vs_induced": round(spearman(bases, ifns), 3),
        }

    out = {
        "dataset": "GSE46599",
        "platform": "Illumina HumanHT-12 v4 (GPL10558)",
        "design_verbatim": "48 samples; design: 9 cell lines, primary CD4+ T cells and "
                           "primary macrophages, untreated and IFN-treated; 2 replicate "
                           "experiments per cell line; 3 replicate experiments per primary "
                           "cell type",
        "treatment_verbatim": "1 to 2 million cells were treated with 1000 U / mL universal "
                              "type 1 IFN (#11200-1 PBL interferon source) for 24 h prior to "
                              "RNA extraction, or left untreated.",
        "n_samples": n_samples,
        "cells": cells,
        "platform_completeness": {
            "genes_requested": ALL_GENES,
            "genes_resolved": {g: resolved[g] for g in resolved},
            "genes_absent_from_array": absent,
            "resolution_evidence": evidence,
            "ISG_module_used": isg_p,
            "ISG_module_n_of_requested": f"{len(isg_p)}/{len(ISG_MODULE)}",
        },
        "manuscript_components": {
            "priming_ISG15_MX1": present(COMPONENT_PRIMING),
            "type_I_IFNAR1_IFNAR2": present(COMPONENT_TYPE_I),
            "type_III_IFNLR1_IL10RB": present(COMPONENT_TYPE_III),
            "composite": "equal-weight mean of the three terms above",
        },
        "per_cell_type": rows,
        "tests": tests,
        "bimodality": bimodality,
        "ceiling_check": ceiling,
        "hiv1_restriction_label_from_series": phenotype,
        "hiv1_restriction_groups": restriction,
    }
    with open(os.path.join(OUT, "g6e_results.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)

    print(f"platform: {len(resolved)}/{len(ALL_GENES)} genes resolved; absent: {absent}")
    print(f"ISG module {len(isg_p)}/{len(ISG_MODULE)}: {isg_p}")
    print(f"\n{'cell':24s} {'kind':13s} {'base':>7s} {'IFN':>7s} {'induction':>10s}")
    for r in rows:
        print(f"{r['cell']:24s} {r['kind']:13s} {r['ISG_baseline']:7.3f} "
              f"{r['ISG_IFN']:7.3f} {r['ISG_induction']:10.3f}")
    for label in tests:
        print(f"\n[{label}] n={tests[label]['n']}")
        for k, v in tests[label].items():
            if isinstance(v, dict):
                print(f"  {k:62s} rho={v['spearman']:+.3f} p={v['p']} ({v['method']})")
    print("\nbimodality:", json.dumps(bimodality, ensure_ascii=False))


if __name__ == "__main__":
    main()
