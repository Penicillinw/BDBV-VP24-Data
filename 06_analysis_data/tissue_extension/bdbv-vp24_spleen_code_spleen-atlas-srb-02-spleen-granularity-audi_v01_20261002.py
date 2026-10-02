"""SRB step 2 - what spleen cell-type resolution is actually public? (criteria R3/R4)

The coverage audit (tools/b6_compartment_coverage_20260926.py) defines the spleen
compartment as exactly three members, and none of them exists in the HPA 154-type
reference:

    splenic red pulp macrophages   (CL:0000874)
    splenic white pulp             (no single CL cell type; white pulp is a tissue region)
    marginal zone B cells          (CL:0000845)

Before "rebuilding the reference score" for the spleen it is therefore necessary
to know whether any public single-cell dataset annotates those populations at all.
This script answers that from the CELLxGENE Discover dataset index, and records the
finest resolution that *is* available for spleen tissue.

Frozen criteria:
  R3a  the dataset index is the live CELLxGENE Discover public index, cached to
       raw/cellxgene_datasets_20260926.json, with the retrieval date recorded;
  R3b  a "spleen dataset" is one whose tissue ontology list contains UBERON:0002106;
  R3c  a spleen-specific population counts as available only if some dataset
       annotates its exact Cell Ontology identifier (no name fuzziness, no
       parent-term substitution);
  R4   the finest available granularity is reported as the full union of cell-type
       labels over all spleen datasets, with per-label dataset counts.

Writes raw/cellxgene_datasets_20260926.json, raw/spleen_dataset_audit.tsv,
raw/spleen_celltype_labels.tsv and raw/spleen_granularity.json.
"""

from __future__ import annotations

import csv
import json
import os
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926")
RAW = os.path.join(OUT, "raw")
CACHE = os.path.join(RAW, "cellxgene_datasets_20260926.json")

API = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
UA = {"User-Agent": "Mozilla/5.0 (research; BDBV-VP24 spleen-reference-rebuild)"}

SPLEEN_TISSUE = "UBERON:0002106"
TARGETS = {
    "marginal zone B cell of spleen": "CL:0000845",
    "splenic red pulp macrophage": "CL:0000874",
    "endothelial cell of venous sinus of red pulp of spleen": "CL:1000397",
    "endothelial cell of venous sinus of spleen": "CL:0002651",
}


def fetch() -> list[dict]:
    os.makedirs(RAW, exist_ok=True)
    if os.path.exists(CACHE):
        return json.load(open(CACHE, encoding="utf-8"))
    req = urllib.request.Request(API, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as r:
        data = json.load(r)
    json.dump(data, open(CACHE, "w", encoding="utf-8"))
    return data


def main() -> None:
    data = fetch()
    print(f"CELLxGENE public datasets indexed: {len(data)}")

    # R3c - availability of the three named spleen populations
    have = {name: [] for name in TARGETS}
    for d in data:
        ids = {c.get("ontology_term_id") for c in (d.get("cell_type") or [])}
        for name, cl in TARGETS.items():
            if cl in ids:
                have[name].append(d["dataset_id"])

    # R3b/R4 - spleen datasets and their label resolution
    spleen = []
    label_datasets: dict[tuple[str, str], set[str]] = {}
    for d in data:
        tis = {t.get("ontology_term_id") for t in (d.get("tissue") or [])}
        if SPLEEN_TISSUE not in tis:
            continue
        spleen.append(d)
        for c in (d.get("cell_type") or []):
            label_datasets.setdefault((c.get("label", ""), c.get("ontology_term_id", "")), set()).add(d["dataset_id"])

    rows = []
    for d in sorted(spleen, key=lambda x: -(x.get("cell_count") or 0)):
        rows.append({
            "dataset_id": d["dataset_id"],
            "title": (d.get("title") or "")[:110],
            "cell_count": d.get("cell_count"),
            "n_cell_types": len(d.get("cell_type") or []),
            "assay": ";".join(a.get("label", "") for a in (d.get("assay") or [])),
            "donor_id": ";".join(d.get("donor_id") or []),
            "development_stage": ";".join(x.get("label", "") for x in (d.get("development_stage") or [])),
            "collection_id": d.get("collection_id", ""),
        })
    with open(os.path.join(RAW, "spleen_dataset_audit.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    lab_rows = [{"cell_type": k[0], "ontology_term_id": k[1], "n_datasets": len(v),
                 "datasets": ";".join(sorted(v))}
                for k, v in sorted(label_datasets.items(), key=lambda kv: -len(kv[1]))]
    with open(os.path.join(RAW, "spleen_celltype_labels.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(lab_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(lab_rows)

    verdict = {name: {"cl_id": cl, "n_datasets": len(have[name]), "dataset_ids": have[name][:10]}
               for name, cl in TARGETS.items()}
    summary = {
        "retrieved": time.strftime("%Y-%m-%d"),
        "n_public_datasets_indexed": len(data),
        "n_spleen_datasets": len(spleen),
        "n_distinct_spleen_cell_type_labels": len(label_datasets),
        "spleen_specific_population_availability": verdict,
        "criterion": ("R3c: available only if a public dataset annotates the exact CL id; "
                      "name similarity and parent terms are not accepted"),
    }
    json.dump(summary, open(os.path.join(RAW, "spleen_granularity.json"), "w", encoding="utf-8"), indent=1)

    print(f"spleen datasets: {len(spleen)}; distinct spleen cell-type labels: {len(label_datasets)}")
    for name, v in verdict.items():
        print(f"  {name:58s} {v['cl_id']:12s} n_datasets = {v['n_datasets']}")
    print("wrote raw/spleen_dataset_audit.tsv, raw/spleen_celltype_labels.tsv, raw/spleen_granularity.json")


if __name__ == "__main__":
    main()
