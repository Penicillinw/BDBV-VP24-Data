"""G9: is the manuscript's Δ ranking an artefact of one atlas and one gene set?

Part A — perturbation of the frozen HPA score: leave-one-gene-out, an extended ISG
module, PC1 weighting, rank-based scoring, and no-log transformation; report rank
correlation against the frozen Δ and the cell types whose rank is unstable.

Part B — independent atlas. The workspace already caches Tabula Sapiens organs,
HLiCA liver, HBCA microglia and the Krasnow lung atlas as h5ad. For matched cell
types, compute the same three components from pseudobulk CPM and compare the
ordering with HPA.

Outputs: analysis/g9_delta_robustness_20260926/
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections import defaultdict

import anndata as ad
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = os.path.join(ROOT, "data", "scRNAseq_raw")
WIDE = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925", "hpa_ifn_landscape_wide.tsv")
GRADIENT = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925", "restriction_gradient.tsv")
OUT = os.path.join(ROOT, "analysis", "g9_delta_robustness_20260926")

FROZEN = {
    "IFN_I_capacity": ["IFNAR1", "IFNAR2"],
    "IFN_III_capacity": ["IFNLR1", "IL10RB"],
    "ISG_priming": ["ISG15", "MX1"],
}
EXTENDED_ISG = ["ISG15", "MX1", "MX2", "OAS1", "BST2", "IFITM1", "IFITM2"]
SIX = [g for v in FROZEN.values() for g in v]

ATLASES = {
    "TS_small_intestine.h5ad": ("cell_type", "total_counts"),
    "TS_large_intestine.h5ad": ("cell_type", "total_counts"),
    "TS_kidney.h5ad": ("cell_type", "total_counts"),
    "TS_testis.h5ad": ("cell_type", "total_counts"),
    "TS_bone_marrow.h5ad": ("cell_type", "total_counts"),
    "TS_Vasculature_cellxgene_source.h5ad": ("cell_type", "total_counts"),
    "HLiCA_hepatocyte.h5ad": ("cell_type", "nCount_RNA"),
    "HBCA_microglia.h5ad": ("cell_type", "total_UMIs"),
    "Krasnow_lung_10X.h5ad": ("cell_type", "nUMI"),
}

# HPA cell type -> list of (atlas file, atlas label). Lineage-level and deliberately lossy.
MAP: dict[str, list[tuple[str, str]]] = {
    "enterocytes": [
        ("TS_small_intestine.h5ad", "enterocyte of epithelium proper of duodenum"),
        ("TS_small_intestine.h5ad", "enterocyte of epithelium proper of ileum"),
        ("TS_small_intestine.h5ad", "enterocyte of epithelium proper of small intestine"),
    ],
    "colonocytes": [("TS_large_intestine.h5ad", "enterocyte of epithelium of large intestine")],
    "goblet cells": [
        ("TS_small_intestine.h5ad", "small intestine goblet cell"),
        ("TS_large_intestine.h5ad", "goblet cell"),
    ],
    "hepatocytes": [("HLiCA_hepatocyte.h5ad", "hepatocyte")],
    "microglia": [("HBCA_microglia.h5ad", "microglial cell")],
    "macrophages": [
        ("Krasnow_lung_10X.h5ad", "alveolar macrophage"),
        ("TS_Vasculature_cellxgene_source.h5ad", "macrophage"),
    ],
    "monocytes": [
        ("TS_bone_marrow.h5ad", "monocyte"),
        ("Krasnow_lung_10X.h5ad", "classical monocyte"),
    ],
    "neutrophils": [("TS_bone_marrow.h5ad", "neutrophil")],
    "vascular endothelial cells": [
        ("TS_Vasculature_cellxgene_source.h5ad", "endothelial cell"),
        ("Krasnow_lung_10X.h5ad", "capillary endothelial cell"),
    ],
    "fibroblasts": [
        ("TS_Vasculature_cellxgene_source.h5ad", "fibroblast"),
        ("Krasnow_lung_10X.h5ad", "fibroblast"),
    ],
    "vascular smooth muscle cells": [
        ("TS_Vasculature_cellxgene_source.h5ad", "smooth muscle cell")
    ],
    "pericytes": [("TS_Vasculature_cellxgene_source.h5ad", "pericyte")],
    "t-cells": [
        ("TS_bone_marrow.h5ad", "CD4-positive, alpha-beta T cell"),
        ("TS_bone_marrow.h5ad", "CD8-positive, alpha-beta T cell"),
    ],
    "b-cells": [("TS_bone_marrow.h5ad", "B cell")],
    "plasma cells": [("TS_bone_marrow.h5ad", "plasma cell")],
    "nk-cells": [("TS_bone_marrow.h5ad", "natural killer cell")],
    "alveolar cells type 2": [("Krasnow_lung_10X.h5ad", "pulmonary alveolar type 2 cell")],
    "peritubular myoid cells": [("TS_testis.h5ad", "peritubular myoid cell")],
    "late spermatids": [("TS_testis.h5ad", "spermatid")],
    "late primary spermatocytes": [("TS_testis.h5ad", "spermatocyte")],
    "proximal tubule cells": [("TS_kidney.h5ad", "kidney epithelial cell")],
}


def spearman(a: list[float], b: list[float]) -> float:
    def rank(x: list[float]) -> list[float]:
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


def component_mean(series_by_gene: dict[str, dict[str, float]], genes: list[str], ct: str) -> float:
    vals = [series_by_gene[g].get(ct, 0.0) for g in genes if g in series_by_gene]
    return sum(vals) / len(vals) if vals else 0.0


def delta_from_components(series_by_gene: dict[str, dict[str, float]], groups: dict[str, list[str]], cts: list[str]) -> dict[str, float]:
    out = {}
    for ct in cts:
        parts = [component_mean(series_by_gene, g, ct) for g in groups.values()]
        out[ct] = sum(parts) / len(parts)
    return out


def part_a() -> dict[str, object]:
    wide = list(csv.DictReader(open(WIDE, encoding="utf-8"), delimiter="\t"))
    cts = [r["cell_type"] for r in wide]
    raw = {
        g: {r["cell_type"]: float(r[g]) for r in wide if r.get(g) not in (None, "")}
        for g in SIX + EXTENDED_ISG
        if all(g in r for r in wide)
    }
    log = {g: {ct: math.log10(raw[g][ct] + 1.0) for ct in cts} for g in raw}

    frozen = delta_from_components(log, FROZEN, cts)

    variants: dict[str, dict[str, float]] = {}
    for drop in SIX:
        groups = {k: [g for g in v if g != drop] for k, v in FROZEN.items()}
        groups = {k: v for k, v in groups.items() if v}
        variants[f"drop_{drop}"] = delta_from_components(log, groups, cts)
    variants["extended_ISG"] = delta_from_components(
        log,
        {"IFN_I_capacity": FROZEN["IFN_I_capacity"], "IFN_III_capacity": FROZEN["IFN_III_capacity"], "ISG_priming": EXTENDED_ISG},
        cts,
    )
    variants["no_log"] = delta_from_components(raw, FROZEN, cts)
    # rank-mean of the three components
    comps = {name: {ct: component_mean(log, genes, ct) for ct in cts} for name, genes in FROZEN.items()}
    ranked = {}
    for name, vals in comps.items():
        order = sorted(cts, key=lambda c: vals[c])
        ranked[name] = {ct: float(i + 1) for i, ct in enumerate(order)}
    variants["rank_mean"] = {ct: sum(ranked[n][ct] for n in ranked) / len(ranked) for ct in cts}
    # PC1 over the six genes
    mat = np.array([[log[g][ct] for g in SIX] for ct in cts], dtype=float)
    mat = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-9)
    u, s, vt = np.linalg.svd(mat, full_matrices=False)
    pc1 = u[:, 0] * s[0]
    if pc1.sum() < 0:
        pc1 = -pc1
    variants["pc1"] = {ct: float(pc1[i]) for i, ct in enumerate(cts)}

    base = [frozen[ct] for ct in cts]
    results = {}
    for name, vals in variants.items():
        rho = spearman(base, [vals[ct] for ct in cts])
        order_base = {ct: i for i, ct in enumerate(sorted(cts, key=lambda c: -frozen[c]))}
        order_var = {ct: i for i, ct in enumerate(sorted(cts, key=lambda c: -vals[c]))}
        shift = {ct: abs(order_base[ct] - order_var[ct]) for ct in cts}
        results[name] = {
            "spearman_vs_frozen": round(rho, 3),
            "max_rank_shift": max(shift.values()),
            "n_cells_shift_gt_20": sum(1 for v in shift.values() if v > 20),
            "largest_shifts": sorted(shift.items(), key=lambda kv: -kv[1])[:8],
        }

    # extremes of the frozen score under the most divergent variant
    worst = min(results, key=lambda k: results[k]["spearman_vs_frozen"])
    top_frozen = [ct for ct in sorted(cts, key=lambda c: -frozen[c])[:15]]
    top_worst = [ct for ct in sorted(cts, key=lambda c: -variants[worst][c])[:15]]

    out_dir_rows = [
        {"cell_type": ct, "delta_frozen": frozen[ct], **{f"delta_{k}": v[ct] for k, v in variants.items()}}
        for ct in cts
    ]
    with open(os.path.join(OUT, "hpa_delta_variants.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out_dir_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(out_dir_rows)

    return {
        "n_cell_types": len(cts),
        "variants": results,
        "most_divergent_variant": worst,
        "top15_frozen": top_frozen,
        "top15_under_most_divergent_variant": top_worst,
        "top15_overlap": len(set(top_frozen) & set(top_worst)),
    }


def part_b() -> dict[str, object]:
    genes_needed = sorted(set(SIX) | {"KPNA1", "KPNA5", "KPNA6"})
    per_atlas: dict[str, dict[str, dict[str, float]]] = {}
    for fname, (label_col, total_col) in ATLASES.items():
        path = os.path.join(SC, fname)
        if not os.path.exists(path):
            continue
        handle = ad.read_h5ad(path, backed="r")
        names = handle.var["feature_name"].astype(str).to_numpy()
        idx = [i for i, g in enumerate(names) if g in genes_needed]
        present = [str(names[i]) for i in idx]
        if not idx:
            handle.file.close()
            print(f"  {fname}: no target genes found, skipped")
            continue
        sub = handle[:, idx].to_memory()
        handle.file.close()
        X = np.asarray(sub.X.todense() if hasattr(sub.X, "todense") else sub.X, dtype=float)
        obs = sub.obs
        if total_col not in obs.columns:
            continue
        totals = obs[total_col].to_numpy(dtype=float)
        labels = obs[label_col].astype(str).to_numpy()
        acc: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
        counts: dict[str, int] = defaultdict(int)
        for j, g in enumerate(present):
            col = X[:, j]
            for lab in np.unique(labels):
                mask = labels == lab
                acc[lab][g] = float(col[mask].sum())
                counts[lab] = int(mask.sum())
        per_atlas[fname] = {
            lab: {
                **{
                    g: math.log10(acc[lab].get(g, 0.0) / max(1.0, float(totals[labels == lab].sum())) * 1e6 + 1.0)
                    for g in present
                },
                "n_cells": counts[lab],
            }
            for lab in acc
        }
        print(f"  {fname}: {len(acc)} cell types, genes present {len(present)}/{len(genes_needed)}")

    hpa = {r["cell_type"]: r for r in csv.DictReader(open(WIDE, encoding="utf-8"), delimiter="\t")}
    pairs = []
    for hpa_ct, targets in MAP.items():
        if hpa_ct not in hpa:
            continue
        hpa_vals = {
            g: math.log10(float(hpa[hpa_ct][g]) + 1.0) if hpa[hpa_ct].get(g) not in (None, "") else 0.0
            for g in SIX
        }
        hpa_delta = sum(
            sum(hpa_vals[g] for g in grp) / len(grp) for grp in FROZEN.values()
        ) / len(FROZEN)
        for fname, lab in targets:
            rec = per_atlas.get(fname, {}).get(lab)
            if not rec:
                continue
            at_delta = sum(
                sum(rec.get(g, 0.0) for g in grp) / len(grp) for grp in FROZEN.values()
            ) / len(FROZEN)
            pairs.append(
                {
                    "hpa_cell_type": hpa_ct,
                    "atlas": fname,
                    "atlas_label": lab,
                    "n_cells": rec["n_cells"],
                    "hpa_delta_log10": round(hpa_delta, 4),
                    "atlas_delta_log10": round(at_delta, 4),
                }
            )

    with open(os.path.join(OUT, "atlas_matched_pairs.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(pairs[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(pairs)

    rho_absolute = spearman(
        [p["hpa_delta_log10"] for p in pairs],
        [p["atlas_delta_log10"] for p in pairs],
    )

    # Cross-atlas absolute CPM is not comparable: each cached atlas is a separate
    # dataset with its own library baseline. Standardise WITHIN each atlas (and
    # within HPA) before comparing orderings, then correlate per atlas as well as
    # pooled.
    def zscore(vals: dict[str, float]) -> dict[str, float]:
        xs = list(vals.values())
        mu = sum(xs) / len(xs)
        sd = math.sqrt(sum((v - mu) ** 2 for v in xs) / len(xs)) or 1.0
        return {k: (v - mu) / sd for k, v in vals.items()}

    hpa_delta_raw = {
        r["cell_type"]: sum(
            sum(math.log10(float(r[g]) + 1.0) for g in grp) / len(grp) for grp in FROZEN.values()
        ) / len(FROZEN)
        for r in csv.DictReader(open(WIDE, encoding="utf-8"), delimiter="\t")
    }
    hpa_z = zscore(hpa_delta_raw)
    atlas_z: dict[str, dict[str, float]] = {}
    for fname, cells in per_atlas.items():
        usable = {
            lab: sum(
                sum(rec.get(g, 0.0) for g in grp) / len(grp) for grp in FROZEN.values()
            ) / len(FROZEN)
            for lab, rec in cells.items()
            if rec.get("n_cells", 0) >= 50
        }
        if len(usable) >= 5:
            atlas_z[fname] = zscore(usable)

    standardised = []
    for p in pairs:
        zc = atlas_z.get(p["atlas"], {}).get(p["atlas_label"])
        if zc is None:
            continue
        standardised.append({**p, "hpa_delta_z": round(hpa_z[p["hpa_cell_type"]], 4), "atlas_delta_z": round(zc, 4)})
    with open(os.path.join(OUT, "atlas_matched_pairs_standardised.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(standardised[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(standardised)

    rho_standardised = spearman(
        [p["hpa_delta_z"] for p in standardised],
        [p["atlas_delta_z"] for p in standardised],
    )
    per_atlas_rho = {}
    for fname in sorted({p["atlas"] for p in standardised}):
        sub = [p for p in standardised if p["atlas"] == fname]
        if len(sub) >= 4:
            per_atlas_rho[fname] = {
                "n_pairs": len(sub),
                "spearman": round(
                    spearman([p["hpa_delta_z"] for p in sub], [p["atlas_delta_z"] for p in sub]), 3
                ),
            }

    return {
        "n_matched_pairs": len(pairs),
        "spearman_absolute_cpm": round(rho_absolute, 3),
        "n_standardised_pairs": len(standardised),
        "spearman_within_dataset_standardised": round(rho_standardised, 3),
        "per_atlas_spearman": per_atlas_rho,
        "pairs": pairs,
        "atlases_used": sorted({p["atlas"] for p in pairs}),
    }


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    print("Part A: HPA score perturbations")
    a = part_a()
    print(json.dumps(a, ensure_ascii=False, indent=1)[:2000])
    print("\nPart B: independent atlas")
    b = part_b()
    print(json.dumps({k: v for k, v in b.items() if k != "pairs"}, ensure_ascii=False, indent=1))
    with open(os.path.join(OUT, "robustness_summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"part_a_hpa_perturbations": a, "part_b_independent_atlas": b}, fh, ensure_ascii=False, indent=1)
    print("\nwrote", os.path.join(OUT, "robustness_summary.json"))


if __name__ == "__main__":
    main()
