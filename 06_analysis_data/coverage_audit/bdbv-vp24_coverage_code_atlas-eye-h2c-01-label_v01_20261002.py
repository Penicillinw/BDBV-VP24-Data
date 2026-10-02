"""H2c step 1 - diagnosis of why the photoreceptor labels failed to connect.

Question (from the 2026-09-26 night handoff, section 14.1): the low-end
compartment test in ``analysis/h2_lowend_atlas_hunt_20260926`` left
**photoreceptor cells** unclosed with the note "label not connected", and left
the **choroid plexus** arm as "methodologically uninformative". This script
establishes, with evidence, *why*, before any re-run is attempted.

Frozen diagnostic criteria (written down before any agreement statistic is
re-computed):

  D1  For every retina-related HPA reference cell type, record whether B5's
      frozen mapping (``b5_harmonised_atlas_20260926/label_mapping.tsv``)
      contains a row.  Absent row => the failure is on the *reference* side.
  D2  For every CL id carried by the new atlases, compute the reflexive
      is_a+part_of ancestor closure and list the HPA reference types it reaches
      under the *frozen* mapping.  Empty => the atlas term is unreachable.
  D3  Re-derive B5's own proposal table (``raw/hpa_cl_candidates.tsv``) for the
      unconnected types and report the ``match_type`` tier recorded there.
      B5 accepts only ``curated`` / ``exact_label`` / a *unique* ``synonym``
      match; ``token_overlap`` was recorded but never accepted automatically
      (``b5_02_mapping.py`` docstring, ``b5_03_curate_mapping.py``).
  D4  For the choroid-plexus arm, test the *mechanical* explanation directly:
      if an atlas contributes exactly one CL term, within-atlas z-scoring
      forces its z to 0, so the arm cannot carry information no matter how many
      cells it contains.  Also test whether any *other* atlas carries the same
      CL term (the only route by which the arm could be recovered without new
      data).

Repair is deliberately NOT performed here: it is tested as a separate,
explicitly labelled sensitivity arm in ``h2c_02_agreement_retina.py``, with the
frozen analysis re-run unchanged alongside it.

Writes only inside ``analysis/h2c_photoreceptor_20260926/``.
Outputs: label_diagnosis.tsv, diagnosis_summary.json
"""

from __future__ import annotations

import csv
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
H2 = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926")
REF = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                   "restriction_gradient.tsv")

sys.path.insert(0, B5)  # reuse the frozen helpers verbatim
from b5_common import Ontology, norm  # noqa: E402

CAND = os.path.join(B5, "raw", "hpa_cl_candidates.tsv")
MAP = os.path.join(B5, "label_mapping.tsv")
ATLAS_LABELS = os.path.join(H2, "raw", "new_atlas_labels.tsv")

RETINA_KINDS = [
    "rod photoreceptor cells",
    "cone photoreceptor cells",
    "retinal amacrine cells",
    "retinal bipolar cells",
    "retinal horizontal cells",
    "m\u00fcller glia",
    "retinal pigment epithelial cells",
    "retinal ganglion cells",
]

# Proposed repair, restricted to the retina family.  This is a *curated* route in
# exactly the sense B5 accepts ("curated ... every entry carries a written
# reason"), but it is MY addition, not part of the frozen B5 artefact, and it is
# therefore only ever used in a clearly-labelled sensitivity arm (h2c_02).
# Cards: HPA label -> (CL id, written reason).  Chosen from the retina family
# only, before any agreement statistic was re-computed, and each is a
# one-to-one choice (CL carries exactly one such term).
PROPOSED = {
    "rod photoreceptor cells": ("CL:0000604", "CL's rod photoreceptor term is 'retinal rod cell' (the only rod term in CL); HPA's label is its plural, region-free form"),
    "cone photoreceptor cells": ("CL:0000573", "CL's cone photoreceptor term is 'retinal cone cell' (the only cone term in CL); HPA's label is its plural, region-free form"),
    "retinal amacrine cells": ("CL:0000561", "CL label is 'amacrine cell'; in CL it sits under retinal neuron"),
    "retinal bipolar cells": ("CL:0000748", "CL label is 'retinal bipolar neuron'"),
    "retinal horizontal cells": ("CL:0000745", "CL label is 'retina horizontal cell' (spelling differs by one letter)"),
    "m\u00fcller glia": ("CL:0000636", "CL label is 'Mueller cell' (ASCII transliteration of M\u00fcller cell); the retinal glia"),
}


def load_tsv(path: str) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    onto = Ontology(os.path.join(B5, "raw", "cl.obo"))
    ref = {r["cell_type"]: r for r in load_tsv(REF)}
    mapping = load_tsv(MAP)
    cand = load_tsv(CAND)
    atlas_labels = load_tsv(ATLAS_LABELS)

    # --- reference-side CL ids actually carried by the frozen mapping --------
    term_to_hpa: dict[str, list[str]] = defaultdict(list)
    mapped_hpa: set[str] = set()
    for r in mapping:
        term_to_hpa[r["cl_id"]].append(r["hpa_cell_type"])
        mapped_hpa.add(r["hpa_cell_type"])

    # --- atlas terms per dataset -------------------------------------------
    by_atlas: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in atlas_labels:
        if r["cl_id"].startswith("CL:"):
            by_atlas[r["atlas"]].append(r)

    # --- D3: candidate tiers recorded by B5 for the retina kinds ------------
    tiers: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    best: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for r in cand:
        name = r["hpa_cell_type"]
        if name not in RETINA_KINDS:
            continue
        tiers[name][r["match_type"]] += 1
        best[name].append((int(r["n_atlas_terms_covered"]), r["cl_id"], r["cl_label"]))

    rows: list[dict[str, object]] = []
    for name in RETINA_KINDS:
        rec = ref.get(name)
        in_map = name in mapped_hpa
        cl_id = next((r["cl_id"] for r in mapping if r["hpa_cell_type"] == name), "")
        tiers_here = dict(tiers.get(name, {}))
        # best token_overlap proposals, ranked by how many atlas terms they cover
        props = sorted(best.get(name, []), key=lambda t: (-t[0], t[1]))[:4]
        rows.append({
            "hpa_cell_type": name,
            "in_reference_table": bool(rec),
            "reference_delta": round(float(rec["delta_restriction"]), 4) if rec else "",
            "reference_rank": int(rec["rank"]) if rec else "",
            "in_b5_frozen_mapping": in_map,
            "b5_cl_id": cl_id,
            "b5_cl_label": onto.name(cl_id) if cl_id else "",
            "b5_norm_query": norm(name),
            "b5_cl_norm_label": norm(onto.name(cl_id)) if cl_id else "",
            "candidate_tiers": ";".join(f"{k}={v}" for k, v in sorted(tiers_here.items())),
            "top_proposals": " | ".join(f"{c}({lbl},covers {n})" for n, c, lbl in props),
            "proposed_cl_id": PROPOSED.get(name, ("", ""))[0],
            "proposed_cl_label": onto.name(PROPOSED[name][0]) if name in PROPOSED else "",
            "proposed_reason": PROPOSED.get(name, ("", ""))[1],
            "proposed_parent_chain": "",
            "proposed_covers": 0,
            "diagnosis": "",
            "reaches_atlas_terms_now": 0,
            "reaches_atlas_terms_now_detail": "",
        })

    # --- D2: which atlas terms are reachable from the frozen mapping --------
    def edges_of(cl: str) -> str:
        """Full direct is_a and part_of edges of a CL term (evidence, not a chain)."""
        rec = onto.terms.get(cl, {})
        isa = [f"{p} {onto.name(p)}" for p in rec.get("is_a", [])]  # type: ignore[union-attr]
        po = [f"{p} {onto.name(p)}" for p in rec.get("parents_alt", [])]  # type: ignore[union-attr]
        return ("is_a: " + ("; ".join(isa) if isa else "none") +
                " || part_of: " + ("; ".join(po) if po else "none"))

    reach_cache: dict[str, set[str]] = {}
    for atlas, recs in by_atlas.items():
        for r in recs:
            cl = r["cl_id"]
            if cl not in reach_cache:
                reach_cache[cl] = {h for t in onto.ancestors(cl)
                                   for h in term_to_hpa.get(t, [])}
    for row in rows:
        name = str(row["hpa_cell_type"])
        hits = []
        for atlas, recs in by_atlas.items():
            for r in recs:
                if name in reach_cache[r["cl_id"]]:
                    hits.append(f"{atlas.split('_')[0]}:{r['cl_id']}"
                                f"({r['cl_label']},{r['n_cells']} cells)")
        row["reaches_atlas_terms_now"] = len(hits)
        row["reaches_atlas_terms_now_detail"] = "; ".join(hits)
        prop = str(row["proposed_cl_id"])
        if prop:
            row["proposed_parent_chain"] = edges_of(prop)
            covered = [f"{a.split('_')[0]}:{r['cl_label']}({r['n_cells']})"
                       for a, recs in by_atlas.items() for r in recs
                       if r["cl_id"] == prop and int(r["n_cells"]) >= 50]
            row["proposed_covers"] = len(covered)
            row["top_proposals"] = ("PROPOSED " + prop + " covers " +
                                    ("; ".join(covered) if covered else "nothing >= 50 cells"))
        if not row["in_reference_table"]:
            row["diagnosis"] = "not an HPA reference cell type"
        elif row["in_b5_frozen_mapping"]:
            row["diagnosis"] = ("mapped by B5; connection failure if any is "
                                "on the atlas side (see reaches_atlas_terms_now)")
        else:
            row["diagnosis"] = (
                "REFERENCE-SIDE GAP: B5 froze no CL id for this HPA label; "
                "every atlas term below is therefore unreachable "
                f"(tiers proposed by B5: {row['candidate_tiers']})")

    with open(os.path.join(HERE, "label_diagnosis.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # --- D4: choroid plexus -------------------------------------------------
    cp_term = "CL:0000706"
    cp_atlases = {a: [r for r in recs if r["cl_id"] == cp_term]
                  for a, recs in by_atlas.items()}
    cp_atlases = {a: v for a, v in cp_atlases.items() if v}
    cp_rec = ref.get("choroid plexus epithelial cells", {})
    cp_summary = {
        "reference_type": "choroid plexus epithelial cells",
        "reference_delta": round(float(cp_rec["delta_restriction"]), 4) if cp_rec else None,
        "mapped_cl": next((r["cl_id"] for r in mapping
                           if r["hpa_cell_type"] == "choroid plexus epithelial cells"), None),
        "atlases_carrying_CL:0000706": {
            a: {"n_cl_terms_in_atlas": len(by_atlas[a]),
                "n_cells_here": int(v[0]["n_cells"]),
                "counts_layer": v[0]["source"]}
            for a, v in cp_atlases.items()},
        "reading": ("within-atlas z of a term that is the atlas's ONLY term is "
                    "identically 0 -> variance-based arm carries no information; "
                    "recovery requires a second atlas carrying the same term"),
        "second_atlas_found": len(cp_atlases) > 1,
    }

    # --- D2b: cross-check the local 19-atlas inventory ----------------------
    local_terms = set()
    inv = os.path.join(B5, "raw", "atlas_celltype_inventory.tsv")
    if os.path.exists(inv):
        for r in load_tsv(inv):
            local_terms.add(r["cl_id"])
    cp_summary["CL:0000706_in_local_19_atlases"] = cp_term in local_terms

    # --- per-atlas label counts, to show granularity of the retina atlases --
    atlas_summary = {
        a: {"n_cl_terms": len(recs),
            "n_terms_ge_50_cells": sum(1 for r in recs if int(r["n_cells"]) >= 50),
            "terms": [f"{r['cl_id']}({r['cl_label']})" for r in recs]}
        for a, recs in by_atlas.items()}

    out = {
        "rule": "diagnosis only; no agreement statistic recomputed here",
        "reference_table_n": len(ref),
        "b5_frozen_mapping_n": len(mapped_hpa),
        "reference_types_without_cl_id": len(ref) - len(set(ref) & mapped_hpa),
        "retina_rows": rows,
        "choroid_plexus": cp_summary,
        "new_atlas_label_summary": atlas_summary,
    }
    with open(os.path.join(HERE, "diagnosis_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)

    print(json.dumps({"reference_types_without_cl_id":
                      out["reference_types_without_cl_id"],
                      "choroid_plexus": cp_summary}, ensure_ascii=False, indent=1))
    print("\nretina rows")
    for r in rows:
        print(f"  {r['hpa_cell_type']:32s} map={str(r['in_b5_frozen_mapping']):5s} "
              f"reaches={r['reaches_atlas_terms_now']:2d}  {r['candidate_tiers']}")


if __name__ == "__main__":
    main()
