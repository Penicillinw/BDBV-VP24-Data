"""G6c replication arm: same frozen question as GSE327707, in an independent
cohort, platform and laboratory (GSE306664, Allen Institute, 10x fixed RNA
profiling).

FROZEN BEFORE RUNNING
  * Modules identical to the other G6 arms (14-gene ISG, 4-gene receptor,
    6-gene housekeeping, ISG15+MX1 priming term, delta_like composite).
  * Cell QC: keep cells with >= 200 UMIs and >= 100 genes. A sample needs
    >= 50 kept cells to enter the analysis; exclusions are reported.
  * Normalisation: log1p(counts per 10,000) per cell, then mean over cells
    (pseudobulk per sample).
  * Primary endpoint: within donor x cell type, Spearman between BASELINE ISG
    module score and the IFN-alpha induction of the same module (n = 5 donors per
    cell type), then pooled with ranks taken within cell type.
  * Secondary: IFN-beta (type I), IFN-lambda1 (type III) and IFN-gamma (type II),
    plus the cross-cell-type version (mean baseline vs mean induction).
  * Controls: housekeeping module, partial correlation controlling it, disjoint
    gene halves, leave-one-donor-out, and the random 14-gene null.

Step 1 (this script, `--extract-pseudobulk`) writes pseudobulk.tsv + sample_qc.tsv.
Step 2 (same script, default) reads those and prints/writes the statistics.
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
import re
import sys

import numpy as np

ROOT = r"G:\本迪布焦研究"
BASE = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926")
RAW = os.path.join(BASE, "raw")
H5 = os.path.join(RAW, "gse306664_h5")
OUT = os.path.join(BASE, "gse306664")

ISG = ["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3",
       "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6"]
RECEPTORS = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB"]
HK = ["ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "HPRT1"]
GENES = list(dict.fromkeys(ISG + RECEPTORS + HK))

CELL_TYPES = ["Bcell", "Monocyte", "NK", "Tcell"]
TREATMENT_OF = {"IFNa": "IFNa", "IFNb": "IFNb", "IFNg": "IFNg", "IFN-L1": "IFN-L1",
                "none": "none", "Fresh": "Fresh",
                "culture_no_stim": "culture_no_stim", "culture_IFNa": "culture_IFNa"}


def load_design() -> dict[str, dict[str, str]]:
    text = open(os.path.join(RAW, "GSE306664_brief.txt"), encoding="utf-8", errors="replace").read()
    design: dict[str, dict[str, str]] = {}
    for block in re.split(r"(?m)^\^SAMPLE = ", text)[1:]:
        gsm = block.split("\n", 1)[0].strip()
        title = re.search(r"(?m)^!Sample_title = (.*)$", block)
        if not title:
            continue
        tokens = title.group(1).split()
        donor = next((t.split("-")[-1] for t in tokens if t.endswith("BW")), "")
        cell_type = next((t for t in tokens if t in {"Bcell", "Monocyte", "NK", "Tcell", "PBMC"}), "")
        treatment = next((t for t in tokens if t in TREATMENT_OF), "")
        design[gsm] = {"donor": donor, "cell_type": cell_type, "treatment": treatment,
                       "title": title.group(1)}
    return design


def extract_pseudobulk() -> None:
    import h5py
    import scipy.sparse as sp

    os.makedirs(OUT, exist_ok=True)
    design = load_design()
    files = sorted(os.listdir(H5))
    genes: list[str] | None = None
    columns: list[str] = []
    matrix: list[np.ndarray] = []
    qc: list[dict] = []
    for i, name in enumerate(files, 1):
        gsm = name.split("_", 1)[0]
        meta = design.get(gsm)
        if meta is None:
            print(f"  !! no design entry for {name}")
            continue
        with h5py.File(os.path.join(H5, name), "r") as fh:
            m = fh["matrix"]
            feat = [x.decode() for x in m["features/name"][:]]
            data = m["data"][:].astype(np.float64)
            indices = m["indices"][:]
            indptr = m["indptr"][:]
            shape = tuple(int(x) for x in m["shape"][:])
            n_genes = int(shape[0])
            cell_counts = np.bincount(
                np.repeat(np.arange(shape[1]), np.diff(indptr)),
                weights=data,
                minlength=shape[1],
            )
            n_genes_per_cell = np.repeat(np.arange(shape[1]), np.diff(indptr))
            genes_per_cell = np.bincount(n_genes_per_cell, minlength=shape[1])
            keep = (cell_counts >= 200) & (genes_per_cell >= 100)
            cells_kept = int(keep.sum())
            if genes is None:
                genes = feat
            assert genes == feat, f"feature order differs in {name}"
            # column-normalise to 10k, log1p, then average over kept cells
            X = sp.csc_matrix((data, indices, indptr), shape=(n_genes, shape[1]))
            keep_idx = np.where(keep)[0]
            pb = np.zeros(n_genes)
            if len(keep_idx):
                sub = X[:, keep_idx].tocsc()
                totals = np.asarray(sub.sum(axis=0)).ravel()
                inv = np.where(totals > 0, 1e4 / np.maximum(totals, 1), 0.0)
                norm = sub.multiply(inv[None, :]).tocsc()
                norm.data = np.log1p(norm.data)
                pb = np.asarray(norm.mean(axis=1)).ravel()
        columns.append(gsm)
        matrix.append(pb)
        qc.append({
            "gsm": gsm, "file": name, "donor": meta["donor"], "cell_type": meta["cell_type"],
            "treatment": meta["treatment"], "cells_total": shape[1], "cells_kept": cells_kept,
            "median_umi": round(float(np.median(cell_counts[keep])) if cells_kept else 0.0, 1),
        })
        if i % 10 == 0:
            print(f"  pseudobulk {i}/{len(files)}", flush=True)
    with open(os.path.join(OUT, "sample_qc.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(qc[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(qc)
    with open(os.path.join(OUT, "pseudobulk.tsv"), "w", encoding="utf-8", newline="") as fh:
        fh.write("gene\t" + "\t".join(columns) + "\n")
        for gi, gene in enumerate(genes or []):
            fh.write(gene + "\t" + "\t".join(f"{matrix[c][gi]:.5f}" for c in range(len(columns))) + "\n")
    print(f"pseudobulk written: {len(genes or [])} genes x {len(columns)} samples")


def spearman(a, b) -> float:
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


def pearson(a, b) -> float:
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    da = math.sqrt(sum((v - ma) ** 2 for v in a))
    db = math.sqrt(sum((v - mb) ** 2 for v in b))
    return num / (da * db) if da and db else float("nan")


def partial(a, b, z) -> float:
    ra, rb, rz = rank_within(a), rank_within(b), rank_within(z)
    rab, raz, rbz = pearson(ra, rb), pearson(ra, rz), pearson(rb, rz)
    denom = math.sqrt((1 - raz ** 2) * (1 - rbz ** 2))
    return (rab - raz * rbz) / denom if denom else float("nan")


def rank_within(v: list[float]) -> list[float]:
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    if "--extract-pseudobulk" in sys.argv:
        extract_pseudobulk()
        return

    genes: list[str] = []
    expr: dict[str, dict[str, float]] = {}
    with open(os.path.join(OUT, "pseudobulk.tsv"), encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)[1:]
        for row in reader:
            genes.append(row[0])
            expr[row[0]] = {header[i]: float(row[i + 1]) for i in range(len(header))}
    qc = list(csv.DictReader(open(os.path.join(OUT, "sample_qc.tsv"), encoding="utf-8"), delimiter="\t"))
    usable = [r for r in qc if int(r["cells_kept"]) >= 50]
    print(f"samples: {len(qc)}  usable (>=50 cells): {len(usable)}")
    for ct in CELL_TYPES:
        donors = sorted({r["donor"] for r in usable if r["cell_type"] == ct})
        print(f"  {ct}: donors {donors}")
    absent = [g for g in GENES if g not in expr]
    if absent:
        print(f"  module genes absent from probe set: {absent}")

    def score(genes_, gsm) -> float:
        vals = [expr[g][gsm] for g in genes_ if g in expr]
        return sum(vals) / len(vals) if vals else float("nan")

    for r in usable:
        r["isg"] = score(ISG, r["gsm"])
        r["isg_2gene"] = score(["ISG15", "MX1"], r["gsm"])
        r["typeI"] = score(["IFNAR1", "IFNAR2"], r["gsm"])
        r["typeIII"] = score(["IFNLR1", "IL10RB"], r["gsm"])
        r["hk"] = score(HK, r["gsm"])
        r["delta_like"] = (r["typeI"] + r["typeIII"] + r["isg_2gene"]) / 3

    def get(cell_type, treatment, donor=None):
        out = {}
        for r in usable:
            if r["cell_type"] == cell_type and r["treatment"] == treatment:
                if donor is None or r["donor"] == donor:
                    out[r["donor"]] = r
        return out

    # ---- which treatments exist per cell type -------------------------------
    treatment_map = {}
    for ct in CELL_TYPES:
        treatment_map[ct] = sorted({r["treatment"] for r in usable if r["cell_type"] == ct})

    # ---- primary: IFN-alpha -------------------------------------------------
    per_cell = {}
    for ct in CELL_TYPES:
        base, ifn = get(ct, "none"), get(ct, "IFNa")
        donors = sorted(set(base) & set(ifn))
        if len(donors) < 3:
            per_cell[ct] = {"n_donors": len(donors), "note": "too few paired donors"}
            continue
        b = [base[d]["isg"] for d in donors]
        ind = [ifn[d]["isg"] - base[d]["isg"] for d in donors]
        per_cell[ct] = {
            "donors": donors,
            "baseline_isg": [round(v, 3) for v in b],
            "induction_isg": [round(v, 3) for v in ind],
            "rho_baseline_isg_vs_induction": round(spearman(b, ind), 3),
            "rho_baseline_2gene_vs_own_induction": round(
                spearman([base[d]["isg_2gene"] for d in donors],
                         [ifn[d]["isg_2gene"] - base[d]["isg_2gene"] for d in donors]), 3),
            "rho_baseline_typeI_vs_isg_induction": round(
                spearman([base[d]["typeI"] for d in donors], ind), 3),
            "rho_delta_like_vs_isg_induction": round(
                spearman([base[d]["delta_like"] for d in donors], ind), 3),
            "rho_housekeeping_vs_isg_induction": round(
                spearman([base[d]["hk"] for d in donors], ind), 3),
        }

    def pooled(treatment: str, cell_types: list[str], exclude: set[str] = frozenset(),
               half: str | None = None, control: str = "none"):
        xa: list[float] = []
        ya: list[float] = []
        hk: list[float] = []
        for ct in cell_types:
            base, trt = get(ct, control), get(ct, treatment)
            donors = [d for d in sorted(set(base) & set(trt)) if d not in exclude]
            if half == "predictor":
                xb = [score(ISG[0::2], base[d]["gsm"]) for d in donors]
                yb = [score(ISG[1::2], trt[d]["gsm"]) - score(ISG[1::2], base[d]["gsm"]) for d in donors]
            elif half == "response":
                xb = [score(ISG[1::2], base[d]["gsm"]) for d in donors]
                yb = [score(ISG[0::2], trt[d]["gsm"]) - score(ISG[0::2], base[d]["gsm"]) for d in donors]
            else:
                xb = [base[d]["isg"] for d in donors]
                yb = [trt[d]["isg"] - base[d]["isg"] for d in donors]
            xa += rank_within(xb)
            ya += rank_within(yb)
            hk += [base[d]["hk"] for d in donors]
        return xa, ya, hk

    primary = {}
    for treatment in ("IFNa", "IFNb", "IFNg", "IFN-L1"):
        x, y, hk = pooled(treatment, CELL_TYPES)
        if len(x) < 6:
            primary[treatment] = {"n_units": len(x), "note": "insufficient paired samples"}
            continue
        rho = spearman(x, y)
        rng = random.Random(20260926)
        draws = 20000
        hits = 0
        for _ in range(draws):
            rng.shuffle(y)
            if abs(spearman(x, y)) >= abs(rho) - 1e-12:
                hits += 1
        primary[treatment] = {
            "n_units": len(x),
            "spearman_baseline_rank_vs_induction_rank": round(rho, 3),
            "p_permutation": round((hits + 1) / (draws + 1), 4),
            "partial_controlling_housekeeping": round(partial(x, y, hk), 3),
        }

    # cross-cell-type view, IFN-alpha
    cross = []
    for ct in CELL_TYPES:
        base, ifn = get(ct, "none"), get(ct, "IFNa")
        donors = sorted(set(base) & set(ifn))
        if not donors:
            continue
        cross.append({
            "cell_type": ct,
            "mean_baseline_isg": round(sum(base[d]["isg"] for d in donors) / len(donors), 3),
            "mean_induction_isg": round(
                sum(ifn[d]["isg"] - base[d]["isg"] for d in donors) / len(donors), 3),
        })
    cross_rho = (round(spearman([c["mean_baseline_isg"] for c in cross],
                                [c["mean_induction_isg"] for c in cross]), 3) if len(cross) >= 3 else None)

    # leave-one-donor-out on the primary endpoint
    all_donors = sorted({r["donor"] for r in usable})
    loo = {}
    for donor in all_donors:
        x, y, _ = pooled("IFNa", CELL_TYPES, exclude={donor})
        loo[donor] = round(spearman(x, y), 3) if len(x) >= 6 else None

    # disjoint halves
    xa, ya, _ = pooled("IFNa", CELL_TYPES, half="predictor")
    xb, yb, _ = pooled("IFNa", CELL_TYPES, half="response")
    halves = {
        "predictor_7_vs_response_other_7": round(spearman(xa, ya), 3) if len(xa) >= 6 else None,
        "predictor_other_7_vs_response_7": round(spearman(xb, yb), 3) if len(xb) >= 6 else None,
    }

    # random gene-set null on the same pooled construction
    background = [g for g in expr if sum(expr[g].values()) / len(expr[g]) >= 1.0]
    rng = random.Random(20260926)
    null_rhos = []
    for _ in range(500):
        take = rng.sample(background, 14)
        xs: list[float] = []
        ys: list[float] = []
        for ct in CELL_TYPES:
            base, ifn = get(ct, "none"), get(ct, "IFNa")
            donors = sorted(set(base) & set(ifn))
            xs += rank_within([sum(expr[g][base[d]["gsm"]] for g in take) / 14 for d in donors])
            ys += rank_within([sum(expr[g][ifn[d]["gsm"]] - expr[g][base[d]["gsm"]] for g in take) / 14
                               for d in donors])
        if len(xs) >= 6:
            null_rhos.append(abs(spearman(xs, ys)))
    null_rhos.sort()
    obs = abs(primary.get("IFNa", {}).get("spearman_baseline_rank_vs_induction_rank", float("nan")))
    null_summary = {
        "draws": len(null_rhos),
        "median_abs_rho": round(null_rhos[len(null_rhos) // 2], 3) if null_rhos else None,
        "p95_abs_rho": round(null_rhos[int(0.95 * len(null_rhos))], 3) if null_rhos else None,
        "observed_percentile": round(100 * sum(1 for v in null_rhos if v <= obs) / len(null_rhos), 1)
        if null_rhos else None,
    }

    # ceiling check
    ceiling = {}
    for ct in CELL_TYPES:
        base, ifn = get(ct, "none"), get(ct, "IFNa")
        donors = sorted(set(base) & set(ifn))
        if not donors:
            continue
        b = [base[d]["isg"] for d in donors]
        p = [ifn[d]["isg"] for d in donors]
        ceiling[ct] = {
            "baseline_range": round(max(b) - min(b), 3),
            "post_ifn_range": round(max(p) - min(p), 3),
            "spearman_baseline_vs_post_level": round(spearman(b, p), 3),
        }

    result = {
        "dataset": "GSE306664",
        "design": "Allen Institute HIRIS-Atlas: healthy donors, FACS-enriched B / Monocyte / NK / T "
                  "cells, 21 h stimulation with IFN-alpha, IFN-beta, IFN-gamma, IFN-lambda1 "
                  "plus untreated controls; 10x fixed RNA profiling",
        "samples_total": len(qc),
        "samples_usable": len(usable),
        "module_genes": {"ISG": ISG, "receptors": RECEPTORS, "housekeeping": HK},
        "genes_absent": absent,
        "treatments_per_cell_type": treatment_map,
        "per_cell_type_ifna": per_cell,
        "primary_pooled": primary,
        "cross_cell_type_ifna": {"rows": cross, "spearman_mean_baseline_vs_mean_induction": cross_rho},
        "leave_one_donor_out_ifna": loo,
        "disjoint_gene_halves_ifna": halves,
        "random_gene_set_null_ifna": null_summary,
        "ceiling_check_ifna": ceiling,
    }
    json.dump(result, open(os.path.join(OUT, "g6c_gse306664_results.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(json.dumps({k: result[k] for k in
                      ("per_cell_type_ifna", "primary_pooled", "cross_cell_type_ifna",
                       "leave_one_donor_out_ifna", "disjoint_gene_halves_ifna",
                       "random_gene_set_null_ifna", "ceiling_check_ifna")},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
