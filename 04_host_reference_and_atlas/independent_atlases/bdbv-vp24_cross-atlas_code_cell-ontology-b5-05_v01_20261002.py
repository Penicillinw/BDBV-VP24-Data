"""B5 step 5 - cross-atlas agreement on harmonised Cell Ontology labels, and
label coverage for the compartments the manuscript predicts.

Routes compared side by side (frozen before the numbers were seen):

  B5-hierarchy   atlas CL terms that are the HPA term or a descendant of it
                 (primary route of B5)
  B5-exact       HPA CL term == atlas CL term, no hierarchy expansion
  g9-replay      the earlier hand-written 21-name / 9-atlas map, reproduced
                 under three value conventions to separate the effect of the
                 label set from the effect of the value convention

Every correlation carries n, a label-permutation p, the smallest |rho|
detectable at 80% power (alpha = 0.05) and the achieved power at the observed
rho.  Because a pooled Spearman across datasets depends on how each dataset is
standardised, two standardisation universes are reported for both routes:

  all_labels   each dataset z-scored over its own full label set
               (HPA over all 154 types; each atlas over its own >=floor labels)
  matched_only each dataset z-scored over the matched labels only
               (removes between-atlas scale differences)

Outputs: agreement.tsv, label_coverage.tsv, raw/b5_matched_pairs.tsv,
         raw/g9_replay_pairs.tsv, raw/b5_summary.json
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
from collections import defaultdict

from b5_common import (
    MODULE,
    Ontology,
    RAW,
    ROOT,
    achieved_power,
    atlas_inventory,
    hpa_delta,
    load_hpa,
    min_detectable_rho,
    spearman,
    spearman_p_value,
    zscore,
)

OUT = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
N_PERM = 20000
SEED = 20260926
random.seed(SEED)

# The compartments the manuscript's prediction names, each anchored on one
# root CL term (frozen from the CL release kept in raw/cl.obo).  "airway" uses
# the lower-respiratory-tract epithelial root because that is the level at
# which the only airway-containing atlas (Krasnow lung) annotates cells; the
# covered terms are listed explicitly in label_coverage.tsv.
COMPARTMENTS = {
    "urothelial": "CL:0000731",
    "airway / lower respiratory epithelium": "CL:0002632",
    "myeloid": "CL:0000763",
    "retinal pigment epithelium": "CL:0002586",
    "photoreceptors": "CL:0000210",
    "Sertoli": "CL:0000216",
    "Leydig": "CL:0000178",
    "choroid plexus": "CL:0000706",
}

# Hand-written map of the earlier g9 probe, reproduced verbatim.
G9_MAP = {
    "enterocytes": [("TS_small_intestine.h5ad", "enterocyte of epithelium proper of duodenum"),
                    ("TS_small_intestine.h5ad", "enterocyte of epithelium proper of ileum"),
                    ("TS_small_intestine.h5ad", "enterocyte of epithelium proper of small intestine")],
    "colonocytes": [("TS_large_intestine.h5ad", "enterocyte of epithelium of large intestine")],
    "goblet cells": [("TS_small_intestine.h5ad", "small intestine goblet cell"),
                     ("TS_large_intestine.h5ad", "goblet cell")],
    "hepatocytes": [("HLiCA_hepatocyte.h5ad", "hepatocyte")],
    "microglia": [("HBCA_microglia.h5ad", "microglial cell")],
    "macrophages": [("Krasnow_lung_10X.h5ad", "alveolar macrophage"),
                    ("TS_Vasculature_cellxgene_source.h5ad", "macrophage")],
    "monocytes": [("TS_bone_marrow.h5ad", "monocyte"),
                  ("Krasnow_lung_10X.h5ad", "classical monocyte")],
    "neutrophils": [("TS_bone_marrow.h5ad", "neutrophil")],
    "vascular endothelial cells": [("TS_Vasculature_cellxgene_source.h5ad", "endothelial cell"),
                                   ("Krasnow_lung_10X.h5ad", "capillary endothelial cell")],
    "fibroblasts": [("TS_Vasculature_cellxgene_source.h5ad", "fibroblast"),
                    ("Krasnow_lung_10X.h5ad", "fibroblast")],
    "vascular smooth muscle cells": [("TS_Vasculature_cellxgene_source.h5ad", "smooth muscle cell")],
    "pericytes": [("TS_Vasculature_cellxgene_source.h5ad", "pericyte")],
    "t-cells": [("TS_bone_marrow.h5ad", "CD4-positive, alpha-beta T cell"),
                ("TS_bone_marrow.h5ad", "CD8-positive, alpha-beta T cell")],
    "b-cells": [("TS_bone_marrow.h5ad", "B cell")],
    "plasma cells": [("TS_bone_marrow.h5ad", "plasma cell")],
    "nk-cells": [("TS_bone_marrow.h5ad", "natural killer cell")],
    "alveolar cells type 2": [("Krasnow_lung_10X.h5ad", "pulmonary alveolar type 2 cell")],
    "peritubular myoid cells": [("TS_testis.h5ad", "peritubular myoid cell")],
    "late spermatids": [("TS_testis.h5ad", "spermatid")],
    "late primary spermatocytes": [("TS_testis.h5ad", "spermatocyte")],
    "proximal tubule cells": [("TS_kidney.h5ad", "kidney epithelial cell")],
}

GENES = [g for v in MODULE.values() for g in v]


def load_pseudobulk() -> dict[tuple[str, str], dict[str, float]]:
    rows = {}
    for r in csv.DictReader(open(os.path.join(RAW, "atlas_pseudobulk.tsv"), encoding="utf-8"), delimiter="\t"):
        rows[(r["atlas"], r["cl_id"])] = {
            "n_cells": int(r["n_cells"]),
            "lib_sum": float(r["lib_sum"]),
            "delta": float(r["delta"]),
            "delta_lognorm": float(r["delta_lognorm_route"]),
            **{f"cpm_{g}": float(r[f"cpm_{g}"]) for g in GENES},
        }
    return rows


def delta_from_cpm(cpm: dict[str, float]) -> float:
    comps = []
    for genes in MODULE.values():
        comps.append(sum(math.log10(cpm.get(g, 0.0) + 1.0) for g in genes) / len(genes))
    return sum(comps) / len(comps)


def pooled_cpm(members: list[dict[str, float]]) -> dict[str, float]:
    """Library-weighted pool of already-CPM values (no extra 1e6 factor)."""
    lib = sum(m["lib_sum"] for m in members)
    if lib <= 0:
        return {g: 0.0 for g in GENES}
    return {g: sum(m[f"cpm_{g}"] * m["lib_sum"] for m in members) / lib for g in GENES}


def perm_p(x: list[float], y: list[float], obs: float, n_perm: int) -> float:
    y = list(y)
    ge = 0
    for _ in range(n_perm):
        random.shuffle(y)
        if abs(spearman(x, y)) >= abs(obs) - 1e-12:
            ge += 1
    return (1 + ge) / (1 + n_perm)


def build_pairs(route, floor, onto, pb, testable, atlas_labels, hpa_delta_all):
    pairs = []
    for m in testable:
        cl, hpa_ct = m["cl_id"], m["hpa_cell_type"]
        for atlas, terms in atlas_labels.items():
            if route == "B5-exact":
                rec = pb.get((atlas, cl))
                matched = [cl] if rec and rec["n_cells"] >= floor else []
            else:
                matched = [
                    t for t in terms
                    if (atlas, t) in pb and cl in onto.ancestors(t) and pb[(atlas, t)]["n_cells"] >= floor
                ]
            if not matched:
                continue
            members = [pb[(atlas, t)] for t in matched]
            pairs.append(
                {
                    "route": route,
                    "hpa_cell_type": hpa_ct,
                    "cl_id": cl,
                    "atlas": atlas,
                    "atlas_terms": ";".join(matched),
                    "n_atlas_terms": len(matched),
                    "n_cells": sum(m["n_cells"] for m in members),
                    "hpa_delta": hpa_delta_all[hpa_ct],
                    "atlas_delta": delta_from_cpm(pooled_cpm(members)),
                }
            )
    return pairs


def main() -> None:
    onto = Ontology(os.path.join(RAW, "cl.obo"))
    hpa = load_hpa()
    pb = load_pseudobulk()
    inv = atlas_inventory()

    mapping = list(csv.DictReader(open(os.path.join(OUT, "label_mapping.tsv"), encoding="utf-8"), delimiter="\t"))
    testable = [m for m in mapping if m["testable"] == "yes"]

    hpa_delta_all = {ct: hpa_delta(row) for ct, row in hpa.items()}
    hpa_z_all = zscore(hpa_delta_all)

    rows = []
    pair_rows: list[dict[str, object]] = []

    for floor in (50, 200, 500):
        atlas_labels = defaultdict(set)
        for r in inv:
            if int(r["n_cells"]) >= floor:
                atlas_labels[r["atlas"]].add(r["cl_id"])
        atlas_stats = {}
        for a, labs in atlas_labels.items():
            vals = [pb[(a, t)]["delta"] for t in labs if (a, t) in pb]
            if len(vals) >= 2:
                mu = sum(vals) / len(vals)
                sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals)) or 1.0
                atlas_stats[a] = (mu, sd)

        for route in ("B5-hierarchy", "B5-exact"):
            pairs = build_pairs(route, floor, onto, pb, testable, atlas_labels, hpa_delta_all)
            if not pairs:
                continue
            # matched-only standardisation
            z_h = zscore({p["hpa_cell_type"]: p["hpa_delta"] for p in pairs})
            z_a = zscore({f"{p['cl_id']}|{p['hpa_cell_type']}": p["atlas_delta"] for p in pairs})
            for p in pairs:
                p["x_matched_only"] = z_h[p["hpa_cell_type"]]
                p["y_matched_only"] = z_a[f"{p['cl_id']}|{p['hpa_cell_type']}"]
                p["x_all_labels"] = hpa_z_all[p["hpa_cell_type"]]
                mu, sd = atlas_stats[p["atlas"]]
                p["y_all_labels"] = (p["atlas_delta"] - mu) / sd

            for universe in ("all_labels", "matched_only"):
                xs = [p[f"x_{universe}"] for p in pairs]
                ys = [p[f"y_{universe}"] for p in pairs]
                n = len(xs)
                rho = spearman(xs, ys)
                tag = f"{route}|floor{floor}|{universe}"
                rows.append(
                    {
                        "analysis": tag,
                        "route": route,
                        "scope": "pooled",
                        "atlas": "",
                        "cell_floor": floor,
                        "z_universe": universe,
                        "n_pairs": n,
                        "spearman": round(rho, 4),
                        "p_perm": round(perm_p(xs, ys, rho, 4000 if n > 60 else N_PERM), 4),
                        "p_perm_n": 4000 if n > 60 else N_PERM,
                        "p_t_approx": round(spearman_p_value(rho, n), 4),
                        "min_detectable_rho_power80": round(min_detectable_rho(n), 3),
                        "achieved_power_at_observed": round(achieved_power(rho, n), 3),
                    }
                )

            if floor == 50:
                # per-atlas
                for atlas in sorted({p["atlas"] for p in pairs}):
                    sub = [p for p in pairs if p["atlas"] == atlas]
                    xs = [p["x_all_labels"] for p in sub]
                    ys = [p["y_all_labels"] for p in sub]
                    rho = spearman(xs, ys) if len(sub) >= 2 else float("nan")
                    rows.append(
                        {
                            "analysis": f"{route}|floor50|per_atlas",
                            "route": route,
                            "scope": "per_atlas",
                            "atlas": atlas,
                            "cell_floor": 50,
                            "z_universe": "all_labels",
                            "n_pairs": len(sub),
                            "spearman": round(rho, 4) if math.isfinite(rho) else "",
                            "p_perm": round(perm_p(xs, ys, rho, 4000), 4) if len(sub) >= 4 else "",
                            "p_perm_n": 4000 if len(sub) >= 4 else "",
                            "p_t_approx": round(spearman_p_value(rho, len(sub)), 4) if len(sub) >= 3 else "",
                            "min_detectable_rho_power80": round(min_detectable_rho(len(sub)), 3),
                            "achieved_power_at_observed": round(achieved_power(rho, len(sub)), 3) if len(sub) >= 4 else "",
                        }
                    )
                # mean of per-atlas rho (each atlas counted once)
                per = [r for r in rows if r["analysis"] == f"{route}|floor50|per_atlas" and r["n_pairs"] >= 5 and r["spearman"] != ""]
                if per:
                    mean_rho = sum(float(r["spearman"]) for r in per) / len(per)
                    rows.append(
                        {
                            "analysis": f"{route}|floor50|mean_of_per_atlas",
                            "route": route,
                            "scope": "mean_of_per_atlas",
                            "atlas": f"atlases with n>=5: {len(per)}",
                            "cell_floor": 50,
                            "z_universe": "all_labels",
                            "n_pairs": sum(int(r["n_pairs"]) for r in per),
                            "spearman": round(mean_rho, 4),
                            "p_perm": "",
                            "p_perm_n": "",
                            "p_t_approx": "",
                            "min_detectable_rho_power80": "",
                            "achieved_power_at_observed": "",
                        }
                    )
                # leave-one-atlas-out
                for atlas in sorted({p["atlas"] for p in pairs}):
                    sub = [p for p in pairs if p["atlas"] != atlas]
                    if len(sub) < 10:
                        continue
                    rho = spearman([p["x_all_labels"] for p in sub], [p["y_all_labels"] for p in sub])
                    rows.append(
                        {
                            "analysis": f"{route}|floor50|leave_one_atlas_out",
                            "route": route,
                            "scope": "leave_one_atlas_out",
                            "atlas": f"without {atlas}",
                            "cell_floor": 50,
                            "z_universe": "all_labels",
                            "n_pairs": len(sub),
                            "spearman": round(rho, 4),
                            "p_perm": "",
                            "p_perm_n": "",
                            "p_t_approx": "",
                            "min_detectable_rho_power80": "",
                            "achieved_power_at_observed": "",
                        }
                    )
                # de-duplication: a generic HPA term and a specific one can both
                # resolve to the same atlas label set, which pseudo-replicates.
                # Collapse to one row per (atlas, atlas_terms) and average the
                # HPA-side z over that group.
                groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
                for p in pairs:
                    groups[(p["atlas"], p["atlas_terms"])].append(p)
                gx = [sum(q["x_all_labels"] for q in g) / len(g) for g in groups.values()]
                gy = [g[0]["y_all_labels"] for g in groups.values()]
                if len(gx) >= 5:
                    rho = spearman(gx, gy)
                    rows.append(
                        {
                            "analysis": f"{route}|floor50|deduplicated_atlas_labelsets",
                            "route": route,
                            "scope": "deduplicated",
                            "atlas": f"{len(gx)} distinct atlas label sets (from {len(pairs)} pairs)",
                            "cell_floor": 50,
                            "z_universe": "all_labels",
                            "n_pairs": len(gx),
                            "spearman": round(rho, 4),
                            "p_perm": round(perm_p(gx, gy, rho, 4000), 4),
                            "p_perm_n": 4000,
                            "p_t_approx": round(spearman_p_value(rho, len(gx)), 4),
                            "min_detectable_rho_power80": round(min_detectable_rho(len(gx)), 3),
                            "achieved_power_at_observed": round(achieved_power(rho, len(gx)), 3),
                        }
                    )
                if route == "B5-hierarchy":
                    pair_rows = pairs

    # ---- g9 replay ------------------------------------------------------
    g9_pairs = []
    name_to_cl = {}
    for r in inv:
        name_to_cl.setdefault((r["atlas"], r["cl_label"]), r["cl_id"])
    for hpa_ct, targets in G9_MAP.items():
        for atlas, atlas_label in targets:
            cl = name_to_cl.get((atlas, atlas_label))
            if cl is None or (atlas, cl) not in pb:
                continue
            g9_pairs.append(
                {
                    "hpa_cell_type": hpa_ct,
                    "atlas": atlas,
                    "atlas_label": atlas_label,
                    "cl_id": cl,
                    "n_cells": pb[(atlas, cl)]["n_cells"],
                    "hpa_delta": hpa_delta_all[hpa_ct],
                    "atlas_delta_rawcpm_route": pb[(atlas, cl)]["delta"],
                    "atlas_delta_lognorm_route": pb[(atlas, cl)]["delta_lognorm"],
                }
            )

    # per-atlas reference distributions for the g9 convention: each atlas
    # z-scored over its own >=50-cell label set, HPA z-scored over all 154 types
    atlas_z_stats = {}
    for atlas in {r["atlas"] for r in inv}:
        vals = [
            pb[(r["atlas"], r["cl_id"])]["delta_lognorm"]
            for r in inv
            if r["atlas"] == atlas and int(r["n_cells"]) >= 50 and (r["atlas"], r["cl_id"]) in pb
        ]
        if len(vals) >= 2:
            mu = sum(vals) / len(vals)
            sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals)) or 1.0
            atlas_z_stats[atlas] = (mu, sd)

    def add_g9_row(tag, values_key, mode):
        if mode == "absolute":
            xs = [p["hpa_delta"] for p in g9_pairs]
            ys = [p[values_key] for p in g9_pairs]
        elif mode == "global_z":
            z_h = zscore({p["hpa_cell_type"]: p["hpa_delta"] for p in g9_pairs})
            z_a = zscore({f"{p['cl_id']}|{p['atlas']}": p[values_key] for p in g9_pairs})
            xs = [z_h[p["hpa_cell_type"]] for p in g9_pairs]
            ys = [z_a[f"{p['cl_id']}|{p['atlas']}"] for p in g9_pairs]
        else:  # per_atlas_z - the convention the earlier probe actually used
            mu_h = sum(hpa_delta_all.values()) / len(hpa_delta_all)
            sd_h = math.sqrt(sum((v - mu_h) ** 2 for v in hpa_delta_all.values()) / len(hpa_delta_all)) or 1.0
            xs, ys = [], []
            for p in g9_pairs:
                mu_a, sd_a = atlas_z_stats[p["atlas"]]
                xs.append((p["hpa_delta"] - mu_h) / sd_h)
                ys.append((p[values_key] - mu_a) / sd_a)
        rho = spearman(xs, ys)
        rows.append(
            {
                "analysis": tag,
                "route": "g9-replay",
                "scope": "pooled",
                "atlas": "",
                "cell_floor": 50,
                "z_universe": {"absolute": "absolute", "global_z": "matched_only", "per_atlas_z": "all_labels"}[mode],
                "n_pairs": len(xs),
                "spearman": round(rho, 4),
                "p_perm": round(perm_p(xs, ys, rho, N_PERM), 4),
                "p_perm_n": N_PERM,
                "p_t_approx": round(spearman_p_value(rho, len(xs)), 4),
                "min_detectable_rho_power80": round(min_detectable_rho(len(xs)), 3),
                "achieved_power_at_observed": round(achieved_power(rho, len(xs)), 3),
            }
        )

    add_g9_row("g9-replay|handmap|raw-cpm-counts", "atlas_delta_rawcpm_route", "absolute")
    add_g9_row("g9-replay|handmap|lognorm-values", "atlas_delta_lognorm_route", "absolute")
    add_g9_row("g9-replay|handmap|raw-cpm-counts|per-atlas-z", "atlas_delta_rawcpm_route", "per_atlas_z")
    add_g9_row("g9-replay|handmap|lognorm-values|per-atlas-z (reproduces 0.256)", "atlas_delta_lognorm_route", "per_atlas_z")

    with open(os.path.join(OUT, "agreement.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(RAW, "b5_matched_pairs.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pair_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(pair_rows)

    with open(os.path.join(RAW, "g9_replay_pairs.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(g9_pairs[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(g9_pairs)

    # ---- compartment coverage -------------------------------------------
    cov_rows = []
    for comp, root in COMPARTMENTS.items():
        for atlas in sorted({r["atlas"] for r in inv}):
            terms = [r for r in inv if r["atlas"] == atlas and root in onto.ancestors(r["cl_id"])]
            cells = sum(int(r["n_cells"]) for r in terms)
            big = [r for r in terms if int(r["n_cells"]) >= 50]
            cov_rows.append(
                {
                    "compartment": comp,
                    "root_cl_id": root,
                    "root_cl_label": onto.name(root),
                    "atlas": atlas,
                    "represented": "yes" if terms else "no",
                    "total_cells": cells,
                    "n_labels": len(terms),
                    "n_labels_ge_50_cells": len(big),
                    "terms": ";".join(
                        f"{r['cl_id']}({r['cl_label']},{r['n_cells']})"
                        for r in sorted(terms, key=lambda r: -int(r["n_cells"]))[:8]
                    ),
                }
            )
    with open(os.path.join(OUT, "label_coverage.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cov_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(cov_rows)

    summary = {
        "n_hpa_types": len(mapping),
        "n_testable_hpa_types": len(testable),
        "n_pairs_hierarchy_floor50": len(pair_rows),
        "n_atlases_with_pairs": len({p["atlas"] for p in pair_rows}),
        "n_g9_replay_pairs": len(g9_pairs),
        "compartment_coverage": {
            comp: {
                "atlases_with_any_cells": sum(1 for r in cov_rows if r["compartment"] == comp and r["represented"] == "yes"),
                "atlases_with_ge50_cells": sum(1 for r in cov_rows if r["compartment"] == comp and int(r["n_labels_ge_50_cells"]) > 0),
                "total_cells_all_atlases": sum(r["total_cells"] for r in cov_rows if r["compartment"] == comp),
            }
            for comp in COMPARTMENTS
        },
    }
    with open(os.path.join(RAW, "b5_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    print(json.dumps(summary, ensure_ascii=False, indent=1))
    print()
    for r in rows:
        if r["scope"] in ("pooled", "mean_of_per_atlas"):
            print(f'{r["analysis"]:62s} n={r["n_pairs"]:4d} rho={r["spearman"]}')


if __name__ == "__main__":
    main()
