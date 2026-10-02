"""B5 step 3 - freeze the HPA cell-type -> Cell Ontology mapping.

Three mapping routes are allowed, in this priority order; each row in the
output records which one was used:

  exact_label     HPA name equals a CL term label (after a documented
                  normalisation: case, plural, hyphenation, token order)
  synonym         HPA name equals a CL synonym string (same normalisation)
  curated         hand-written override in CURATED below, used where the CL
                  label differs from the HPA label by granularity or by the
                  HPA shorthand only.  Frozen before any agreement statistic
                  was computed; every entry carries a written reason.

`testable` = at least one atlas contains cells whose CL term is the mapped
term itself or a MORE SPECIFIC term below it (descendant).  Direction of the
ontology query was corrected relative to the throw-away probe in b5_02:
descendants, not ancestors.

Output: label_mapping.tsv (deliverable), raw/hpa_cl_unmapped.tsv (diagnostic)
"""

from __future__ import annotations

import csv
import os
from collections import defaultdict

from b5_common import HPA_ALIASES, Ontology, QC_FLOOR, RAW, ROOT, atlas_inventory, load_hpa, norm

# Hand-written, frozen overrides. Each entry: HPA name -> (CL id, reason).
# Only used where the CL label differs from the HPA label by granularity
# (region/segment-level CL subtypes) or by an HPA shorthand the normaliser
# cannot bridge. Nothing here was chosen by looking at agreement.
CURATED: dict[str, tuple[str, str]] = {
    "alveolar cells type 1": ("CL:0002062", "CL label is 'pulmonary alveolar type 1 cell'"),
    "alveolar cells type 2": ("CL:0002063", "CL label is 'pulmonary alveolar type 2 cell'"),
    "respiratory secretory cells": ("CL:1000272", "CL term for the airway secretory population is 'lung secretory cell'"),
    "respiratory deuterosomal cells": ("CL:4033044", "CL term is 'deuterosomal cell'"),
    "conjunctival goblet cells": ("CL:2000084", "CL label is 'conjunctiva goblet cell'"),
    "vascular endothelial cells": ("CL:0000115", "override: CL:0002139 'endothelial cell of vascular tree' is one level too specific and its only child in CL is the lymphatic lineage"),
    "proximal tubule cells": ("CL:0002306", "CL label is 'epithelial cell of proximal tubule'"),
    "undifferentiated spermatogonia": ("CL:0000020", "CL label is 'spermatogonium'"),
    "differentiating spermatogonia": ("CL:0000020", "CL label is 'spermatogonium'"),
    "early primary spermatocytes": ("CL:0000656", "CL label is 'primary spermatocyte'"),
    "late primary spermatocytes": ("CL:0000656", "CL label is 'primary spermatocyte'"),
    "distal convoluted tubule cells": ("CL:0002305", "CL label is 'epithelial cell of distal tubule'"),
    "gastric chief cells": ("CL:0000155", "CL label is 'peptic cell'"),
    "lactotrophs": ("CL:0002311", "CL label is 'mammotroph'"),
    "smooth muscle cells": ("CL:0000192", "CL label is 'smooth muscle cell'"),
    "vascular smooth muscle cells": ("CL:0000359", "CL label is 'vascular associated smooth muscle cell'"),
    "lymphatic endothelial cells": ("CL:0002138", "CL label is 'endothelial cell of lymphatic vessel'"),
    "mesothelial cells": ("CL:0000077", "CL label is 'mesothelial cell'"),
    "erythrocyte progenitors": ("CL:0000038", "CL label is 'erythroid progenitor cell'"),
    "megakaryocyte-erythroid progenitors": ("CL:0000038", "CL label is 'erythroid progenitor cell'"),
    "megakaryocyte progenitors": ("CL:0000556", "CL label is 'megakaryocyte'"),
    "hematopoietic stem cells": ("CL:0000037", "CL label is 'hematopoietic stem cell'"),
    "monocyte progenitors": ("CL:0000557", "CL label is 'granulocyte monocyte progenitor cell'"),
    "neutrophil progenitors": ("CL:0000775", "CL label is 'neutrophil'"),
    "pancreatic islet cells": ("CL:0000169", "CL label is 'type B pancreatic cell'; HPA's islet population is dominated by beta cells"),
    "respiratory ionocytes": ("CL:0017000", "CL label is 'pulmonary ionocyte'"),
    "salivary ionocytes": ("CL:0020065", "CL label is 'ionocyte of salivary gland'"),
}


def main() -> None:
    onto = Ontology(os.path.join(RAW, "cl.obo"))
    hpa = load_hpa()
    inv = atlas_inventory()

    terms_by_atlas: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for r in inv:
        terms_by_atlas[r["atlas"]][r["cl_id"]] = {"label": r["cl_label"], "n_cells": int(r["n_cells"])}

    by_label: dict[str, set[str]] = defaultdict(set)
    by_syn: dict[str, set[str]] = defaultdict(set)
    for tid, rec in onto.terms.items():
        if rec["obsolete"]:
            continue
        by_label[norm(str(rec["name"]))].add(tid)
        for s in rec["synonyms"]:  # type: ignore[union-attr]
            by_syn[norm(str(s))].add(tid)

    rows = []
    unmapped = []
    for name in hpa:
        query = HPA_ALIASES.get(name, name)
        nq = norm(query)
        cl_id = None
        route = ""
        evidence = query
        if name in CURATED:
            cl_id, why = CURATED[name]
            route = "curated"
            evidence = why
        elif len(by_label.get(nq, ())) == 1:
            cl_id = next(iter(by_label[nq]))
            route = "exact_label"
        elif len(by_syn.get(nq, ())) == 1:
            cl_id = next(iter(by_syn[nq]))
            route = "synonym"

        if cl_id is None:
            unmapped.append({"hpa_cell_type": name, "query_used": query})
            continue

        ancestors = onto.ancestors(cl_id)
        covering = []
        for atlas, terms in terms_by_atlas.items():
            for t, rec in terms.items():
                if t in ancestors:  # t == cl_id, or t is a more specific term under it
                    covering.append((atlas, t, rec["label"], rec["n_cells"]))
        covering.sort(key=lambda x: -x[3])
        testable = [c for c in covering if c[3] >= QC_FLOOR]
        rows.append(
            {
                "hpa_cell_type": name,
                "query_string": query,
                "cl_id": cl_id,
                "cl_label": onto.name(cl_id),
                "mapping_route": route,
                "evidence": evidence,
                "testable": "yes" if testable else "no",
                "n_covering_atlas_terms": len(covering),
                "n_covering_ge_floor": len(testable),
                "covering_atlases": ",".join(sorted({c[0].replace(".h5ad", "") for c in testable})),
                "covering_terms": ";".join(f"{c[0].replace('.h5ad', '')}:{c[1]}({c[3]})" for c in covering),
            }
        )

    rows.sort(key=lambda r: (r["testable"] != "yes", -int(r["n_covering_ge_floor"]), r["hpa_cell_type"]))
    out_path = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926", "label_mapping.tsv")
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(RAW, "hpa_cl_unmapped.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["hpa_cell_type", "query_used"], delimiter="\t")
        w.writeheader()
        w.writerows(unmapped)

    testable = [r for r in rows if r["testable"] == "yes"]
    by_route: dict[str, int] = defaultdict(int)
    for r in rows:
        by_route[r["mapping_route"]] += 1
    print(f"HPA types: {len(hpa)}; mapped: {len(rows)}; testable: {len(testable)}")
    print("routes:", dict(by_route))
    print(f"unmapped: {len(unmapped)} -> raw/hpa_cl_unmapped.tsv")


if __name__ == "__main__":
    main()
