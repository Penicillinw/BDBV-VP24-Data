"""Third anchor, cross-sample scope: does the pre-treatment interferon-stimulated
LEVEL, rather than the induction it undergoes, track a functional viral outcome?

Dataset  GSE107549 (GPL10558, Illumina HumanHT-12 v4): PBMC RNA from 21 chronically
         HIV-infected, ART-suppressed subjects on pegylated IFN-alpha2a, visits v3
         (ART only, pre-immunotherapy), v5 and v8; response category per subject.

Why this is the right test for the manuscript
    The manuscript's pendant question is whether the cost of a defective VP24 scales
    with the absolute interferon-stimulated level a cell attains or with the fold
    induction it mounts.  PI ruling (2026-09-27): for the third anchor a cross-sample
    (subject/treatment) contrast is acceptable, not only a cross-cell-type one.  In
    this dataset the functional readout is virological control, and the published
    analysis attributes it to the pre-immunotherapy baseline state.

Criteria frozen before any value was read
    C1  Same platform annotation and same module definitions as the project's
        GSE46599 line (ISG17 / PRIMING / RECEPTOR / HK).
    C2  Module score per sample = mean over the module's genes of the mean of that
        gene's probes, on the log2 scale as deposited.  No cross-sample z-scoring,
        so that the score stays interpretable as a level.
    C3  Primary contrast: baseline (v3) module score, responders (RC) versus
        non-responders (RN).  Report AUC and an exact permutation p.  The rival
        quantity is the induction (v5-v3 and v8-v3), tested identically.  The
        comparison the manuscript cares about is AUC(level) versus AUC(fold).
    C4  Null distribution by exhaustive enumeration of the label assignments when
        the space allows, otherwise Monte Carlo with a fixed seed (reported).
    C5  If the baseline level separates the groups in the direction the source
        publication reports, it is recorded as positive evidence for the level
        reading; otherwise it is recorded as evidence against it, without spin.
    C6  Response category is taken only from the series metadata, never inferred.

Writes only into this directory.
"""

import csv
import gzip
import itertools
import json
import math
import os
import random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "third_anchor_gse107549_20260926")
MATRIX = os.path.join(OUT, "raw", "GSE107549_series_matrix.txt.gz")
ANNOT = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926", "raw", "GPL10558.annot.gz")

ISG17 = ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "IFIT1", "IFIT2", "IFIT3",
         "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6",
         "DDX58", "IFIH1"]
PRIMING = ["ISG15", "MX1"]
RECEPTOR = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB"]
HK = ["ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "HPRT1"]
MODULES = {"ISG17": ISG17, "PRIMING": PRIMING, "RECEPTOR": RECEPTOR, "HK": HK}

SEED = 20260927
MC_DRAWS = 200000


def read_annotation():
    """probe -> gene symbol (first symbol when several are listed)."""
    probe2gene = {}
    with gzip.open(ANNOT, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                header = fh.readline().rstrip("\n").split("\t")
                gi = header.index("Gene symbol")
                ii = header.index("ID")
                for row in fh:
                    if row.startswith("!platform_table_end"):
                        break
                    parts = row.rstrip("\n").split("\t")
                    if len(parts) <= max(gi, ii):
                        continue
                    sym = parts[gi].split("///")[0].strip()
                    if sym:
                        probe2gene[parts[ii]] = sym
                break
    return probe2gene


def read_matrix():
    """returns (sample ids in order, {probe: [values]}, meta dict of lists)."""
    meta = {}
    gsm, values = None, {}
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_"):
                key, rest = line.rstrip("\n").split("\t", 1)
                meta.setdefault(key, []).append([x.strip('"') for x in rest.split("\t")])
            if line.startswith("!series_matrix_table_begin"):
                header = [x.strip('"') for x in fh.readline().rstrip("\n").split("\t")]
                gsm = header[1:]
                for row in fh:
                    if row.startswith("!series_matrix_table_end"):
                        break
                    parts = row.rstrip("\n").split("\t")
                    probe = parts[0].strip('"')
                    try:
                        values[probe] = [float(x) for x in parts[1:]]
                    except ValueError:
                        continue
                break

    def field(prefix):
        for row in meta.get("!Sample_characteristics_ch1", []):
            if row and row[0].lower().startswith(prefix):
                return [x.split(":", 1)[1].strip() for x in row]
        return None

    return gsm, values, {
        "visit": field("visit"),
        "response": field("response"),
        "patient": field("patient id"),
    }


def module_scores(values, probe2gene, genes):
    """per-sample module score, log2 scale, no cross-sample standardisation."""
    per_gene = {}
    for gene in genes:
        probes = [p for p, s in probe2gene.items() if s == gene and p in values]
        if not probes:
            continue
        n = len(values[probes[0]])
        per_gene[gene] = [sum(values[p][i] for p in probes) / len(probes) for i in range(n)]
    if not per_gene:
        return None, []
    n = len(next(iter(per_gene.values())))
    score = [sum(per_gene[g][i] for g in per_gene) / len(per_gene) for i in range(n)]
    return score, sorted(per_gene)


def auc(cases, controls):
    """P(case > control) + 0.5 P(tie), rank based, no scipy."""
    wins = ties = 0
    for c in cases:
        for k in controls:
            if c > k:
                wins += 1
            elif c == k:
                ties += 1
    n = len(cases) * len(controls)
    return (wins + 0.5 * ties) / n if n else float("nan")


def permutation_p(pool, n_case, observed, stat):
    """exact when the label space is small, else Monte Carlo with a fixed seed."""
    n = len(pool)
    total = math.comb(n, n_case)
    if total <= 200000:
        ge = 0
        for idx in itertools.combinations(range(n), n_case):
            case = [pool[i] for i in idx]
            control = [pool[i] for i in range(n) if i not in idx]
            if stat(case, control) >= observed - 1e-12:
                ge += 1
        return ge / total, "exact", total
    rng = random.Random(SEED)
    ge = 0
    for _ in range(MC_DRAWS):
        idx = set(rng.sample(range(n), n_case))
        case = [pool[i] for i in sorted(idx)]
        control = [pool[i] for i in range(n) if i not in idx]
        if stat(case, control) >= observed - 1e-12:
            ge += 1
    return (ge + 1) / (MC_DRAWS + 1), "monte-carlo", MC_DRAWS


def main():
    probe2gene = read_annotation()
    gsm, values, meta = read_matrix()
    visit, response, patient = meta["visit"], meta["response"], meta["patient"]
    print("samples", len(gsm), "| probes", len(values), "| annotated probes", len(probe2gene))

    scores, genes_used = {}, {}
    for name, genes in MODULES.items():
        scores[name], genes_used[name] = module_scores(values, probe2gene, genes)
        print("  module %-9s genes resolved %2d/%2d" % (name, len(genes_used[name]), len(genes)))

    idx_by_subject = {}
    for i, (p, v) in enumerate(zip(patient, visit)):
        idx_by_subject.setdefault(p, {})[v] = (i, response[i])

    rows = []
    for subj, visits in sorted(idx_by_subject.items()):
        r = visits.get("v3", (None, None))[1]
        row = {"subject": subj, "response": r,
               "has_v3": "v3" in visits, "has_v5": "v5" in visits, "has_v8": "v8" in visits}
        for name in MODULES:
            if "v3" in visits and scores[name]:
                i = visits["v3"][0]
                row[name + "_v3"] = scores[name][i]
                if "v5" in visits:
                    row[name + "_v5"] = scores[name][visits["v5"][0]]
                if "v8" in visits:
                    row[name + "_v8"] = scores[name][visits["v8"][0]]
        rows.append(row)

    results = {}
    for name in ("ISG17", "PRIMING", "RECEPTOR", "HK"):
        rc = [r[name + "_v3"] for r in rows if r["response"] == "RC" and name + "_v3" in r]
        rn = [r[name + "_v3"] for r in rows if r["response"] == "RN" and name + "_v3" in r]
        if len(rc) < 3 or len(rn) < 3:
            continue
        a = auc(rc, rn)
        pool = rc + rn
        p, how, space = permutation_p(pool, len(rc), max(a, 1 - a), auc)
        p_one = permutation_p(pool, len(rc), a, auc)[0]
        folds_rc, folds_rn = [], []
        for r in rows:
            k5, k3 = name + "_v5", name + "_v3"
            if r["response"] == "RC" and k5 in r:
                folds_rc.append(r[k5] - r[k3])
            if r["response"] == "RN" and k5 in r:
                folds_rn.append(r[k5] - r[k3])
        a_fold = auc(folds_rc, folds_rn) if len(folds_rc) >= 3 and len(folds_rn) >= 3 else None
        results[name] = {
            "n_RC": len(rc), "n_RN": len(rn),
            "mean_RC_v3": sum(rc) / len(rc), "mean_RN_v3": sum(rn) / len(rn),
            "auc_level_v3": a, "p_perm_level": p, "p_perm_level_one_sided": p_one,
            "permutation": how, "permutation_space": space,
            "n_RC_fold": len(folds_rc), "n_RN_fold": len(folds_rn),
            "auc_fold_v5_v3": a_fold,
        }

    # Null distribution for the AUC: random gene modules of the same size, drawn
    # from the genes this array carries, scored exactly like ISG17.
    rng = random.Random(SEED)
    genes_available = sorted({g for g in probe2gene.values() if any(
        p in values for p, s in probe2gene.items() if s == g)})
    rc_idx = [i for i, r in enumerate(rows)
              if r["response"] == "RC" and "ISG17_v3" in r]
    rn_idx = [i for i, r in enumerate(rows)
              if r["response"] == "RN" and "ISG17_v3" in r]
    null = []
    for _ in range(2000):
        genes = rng.sample(genes_available, len(ISG17))
        score, used = module_scores(values, probe2gene, genes)
        if score is None or len(used) < len(ISG17):
            continue
        case = [score[i] for i in rc_idx]
        control = [score[i] for i in rn_idx]
        null.append(auc(case, control))
    null.sort()
    if null and "ISG17" in results:
        obs = results["ISG17"]["auc_level_v3"]
        pct = 100.0 * sum(1 for v in null if v < obs) / len(null)
        results["ISG17"]["null_median_auc"] = null[len(null) // 2]
        results["ISG17"]["null_p95_auc"] = null[int(0.95 * (len(null) - 1))]
        results["ISG17"]["null_max_auc"] = null[-1]
        results["ISG17"]["percentile_of_observed"] = pct
        results["null_random_modules"] = {"n_draws": len(null), "seed": SEED,
                                          "module_size": len(ISG17)}

    with open(os.path.join(OUT, "subject_scores.tsv"), "w", encoding="utf-8", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(OUT, "third_anchor_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"modules": results, "genes_used": genes_used,
                   "n_samples": len(gsm), "n_subjects": len(rows)}, fh,
                  ensure_ascii=False, indent=1)

    for name, rec in results.items():
        if name == "null_random_modules":
            continue
        print("\n%s  RC=%d RN=%d" % (name, rec["n_RC"], rec["n_RN"]))
        print("  baseline v3 mean: RC %.4f  RN %.4f" % (rec["mean_RC_v3"], rec["mean_RN_v3"]))
        print("  AUC(level) = %.4f   p = %.4g (%s, space %d)"
              % (rec["auc_level_v3"], rec["p_perm_level"], rec["permutation"], rec["permutation_space"]))
        print("  AUC(fold v5-v3) = %s" % ("%.4f" % rec["auc_fold_v5_v3"]
                                           if rec["auc_fold_v5_v3"] is not None else "n/a"))


if __name__ == "__main__":
    main()
