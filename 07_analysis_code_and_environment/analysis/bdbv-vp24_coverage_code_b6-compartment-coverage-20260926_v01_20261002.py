"""B6: coverage audit of the compartments the manuscript's partition names.

The score exists only for the 154 cell types in the HPA single-cell-type table. This
script asks, for every compartment the manuscript's prediction names (and for the
compartments its motivation cites), whether a matching cell type is present, and it
records the gaps explicitly so the manuscript can state them.
"""

from __future__ import annotations

import csv
import json
import os

ROOT = r"G:\本迪布焦研究"
GRADIENT = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                        "restriction_gradient.tsv")
OUT = os.path.join(ROOT, "analysis", "b6_coverage_20260926")

# compartments named in the prediction, with the HPA cell types that stand in for them.
# membership lists are written out explicitly (no substring guessing) and frozen here.
COMPARTMENTS = {
    "urothelium": ["urothelial cells"],
    "enteric epithelium": ["enterocytes", "colonocytes", "goblet cells", "paneth cells",
                           "tuft cells", "enteric stem cells",
                           "enteric transient amplifying cells",
                           "foveolar cells", "mucous neck cells", "parietal cells",
                           "gastric chief cells", "gastric progenitor cells"],
    "airway epithelium": ["respiratory basal cells", "respiratory ciliated cells",
                          "respiratory secretory cells", "respiratory deuterosomal cells",
                          "respiratory ionocytes", "alveolar cells type 1",
                          "alveolar cells type 2", "transitional alveolar cells"],
    "ocular surface epithelium": ["ocular epithelial cells", "conjunctival goblet cells"],
    "hepatobiliary": ["hepatocytes", "kupffer cells", "cholangiocytes",
                      "hepatic stellate cells"],
    "myeloid": ["monocytes", "macrophages", "kupffer cells", "microglia", "pdcs", "cdc",
                "neutrophils", "mast cells", "hofbauer cells", "monocyte progenitors",
                "neutrophil progenitors"],
    "retinal pigment epithelium": ["retinal pigment epithelial cells"],
    "photoreceptors": ["rod photoreceptor cells", "cone photoreceptor cells"],
    "retinal neurons": ["retinal amacrine cells", "retinal bipolar cells",
                        "retinal horizontal cells", "retinal ganglion cells"],
    "Muller glia": ["müller glia", "bergmann glia"],
    "seminiferous epithelium": ["sertoli cells", "undifferentiated spermatogonia",
                                "differentiating spermatogonia", "early primary spermatocytes",
                                "late primary spermatocytes", "early spermatids",
                                "late spermatids"],
    "interstitial cells of the testis": ["leydig cells", "peritubular myoid cells"],
    "epididymis": ["epididymal basal cells", "epididymal clear cells",
                   "epididymal principal cells",
                   "epididymal efferent duct absorptive cells",
                   "epididymal efferent duct ciliated cells"],
    "ventricular system": ["choroid plexus epithelial cells", "ependymal cells"],
    "kidney tubule": ["proximal tubule cells", "loop of henle epithelial cells",
                      "distal convoluted tubule cells", "podocytes",
                      "renal collecting duct principal cells",
                      "renal collecting duct intercalated cells",
                      "renal connecting tubule cells"],
    "endothelium": ["vascular endothelial cells", "lymphatic endothelial cells"],
    "plasma cells": ["plasma cells"],
    "trophoblast": ["extravillous trophoblasts", "cytotrophoblasts",
                    "syncytiotrophoblasts", "migrating cytotrophoblasts"],
    # compartments the motivation cites but that the reference may not cover
    "spleen": ["splenic red pulp macrophages", "splenic white pulp", "marginal zone B cells"],
    "thymus": ["thymocytes", "medullary thymic epithelial cells"],
    "adrenal gland": ["adrenal cortex cells", "adrenal medulla cells"],
    "pituitary": ["somatotrophs", "lactotrophs", "corticotrophs", "gonadotrophs",
                  "thyrotrophs", "pituicytes/fscs", "pituitary stem cells"],
    "skeletal muscle": ["myonuclei", "myosatellite cells",
                        "fibro-adipogenic progenitors"],
    "cardiac muscle": ["cardiomyocytes", "epicardial cells"],
    "adipose": ["adipocytes"],
    "skin": ["basal keratinocytes", "suprabasal keratinocytes", "melanocytes",
             "fibroblasts", "schwann cells"],
}


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    rows = list(csv.DictReader(open(GRADIENT, encoding="utf-8"), delimiter="\t"))
    present = {r["cell_type"]: int(r["rank"]) for r in rows}
    delta = {r["cell_type"]: float(r["delta_restriction"]) for r in rows}

    out = []
    for comp, members in COMPARTMENTS.items():
        found = [m for m in members if m in present]
        missing = [m for m in members if m not in present]
        ranks = sorted(present[m] for m in found)
        dvals = [delta[m] for m in found]
        out.append({
            "compartment": comp,
            "n_members_listed": len(members),
            "n_present_in_hpa": len(found),
            "n_missing": len(missing),
            "missing_members": ";".join(missing),
            "coverage": round(len(found) / len(members), 2) if members else 0.0,
            "best_rank": ranks[0] if ranks else None,
            "worst_rank": ranks[-1] if ranks else None,
            "median_delta": round(sorted(dvals)[len(dvals) // 2], 3) if dvals else None,
            "present_members": ";".join(found),
        })

    out.sort(key=lambda r: (r["n_present_in_hpa"] == 0, r["coverage"]))
    with open(os.path.join(OUT, "compartment_coverage.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    absent = [r["compartment"] for r in out if r["n_present_in_hpa"] == 0]
    partial = [r["compartment"] for r in out if 0 < r["coverage"] < 1]
    print(f"HPA cell types available: {len(rows)}")
    print(f"compartments audited: {len(out)}; fully absent: {len(absent)}; "
          f"partial: {len(partial)}")
    print("absent:", absent)
    print()
    print(f"{'compartment':34s} {'present/listed':>15s} {'best':>5s} {'worst':>6s} "
          f"{'medDelta':>9s}")
    for r in out:
        print(f"{r['compartment']:34s} "
              f"{str(r['n_present_in_hpa']) + '/' + str(r['n_members_listed']):>15s} "
              f"{str(r['best_rank']):>5s} {str(r['worst_rank']):>6s} "
              f"{str(r['median_delta']):>9s}")
    json.dump({"n_hpa_cell_types": len(rows), "compartments": out,
               "absent_compartments": absent, "partial_compartments": partial},
              open(os.path.join(OUT, "coverage_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("\nwrote", os.path.join(OUT, "compartment_coverage.tsv"))


if __name__ == "__main__":
    main()
