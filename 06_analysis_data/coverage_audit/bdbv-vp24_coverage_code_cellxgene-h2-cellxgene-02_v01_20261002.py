"""H2 line: scan CellxGene Discover for human single-cell atlases covering the
low-end compartments (urothelium, RPE, photoreceptors, Sertoli, choroid plexus,
spleen).

Writes raw API responses plus a filtered TSV. Read-only with respect to the
rest of the repository: everything lands in this directory.
"""
from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
RAW.mkdir(parents=True, exist_ok=True)

API = "https://api.cellxgene.cziscience.com/curation/v1"
UA = {"User-Agent": "Mozilla/5.0 (research; BDBV-VP24 atlascan)"}

# tissue labels as they appear in CellxGene metadata (lower-case match)
TARGETS = {
    "urothelium": ["bladder", "ureter", "renal pelvis", "urothel"],
    "rpe": ["retinal pigment", "rpe", "eye", "retina"],
    "photoreceptor": ["retina", "eye", "photoreceptor"],
    "sertoli": ["testis", "testicle", "seminiferous"],
    "choroid_plexus": ["choroid plexus", "brain", "ventricle"],
    "spleen": ["spleen", "splenic"],
}


def get(url: str, timeout: int = 90):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def main() -> None:
    collections = get(f"{API}/collections?visibility=PUBLIC")
    (RAW / "cellxgene_collections.json").write_text(
        json.dumps(collections, indent=1), encoding="utf-8"
    )
    print(f"collections: {len(collections)}")

    rows = []
    hits = {k: [] for k in TARGETS}
    for ci, coll in enumerate(collections):
        cid = coll.get("collection_id")
        cname = coll.get("name", "")
        doi = coll.get("doi", "")
        for ds in coll.get("datasets", []) or []:
            tissues = ds.get("tissue") or []
            taxa = ds.get("organism") or []
            if not any("Homo sapiens" in str(t) or t == "Homo sapiens" for t in taxa):
                continue
            blob = " ".join(str(t) for t in tissues).lower()
            matched = [k for k, kws in TARGETS.items() if any(kw in blob for kw in kws)]
            if not matched:
                continue
            row = {
                "collection_id": cid,
                "collection_name": cname,
                "collection_doi": doi,
                "dataset_id": ds.get("dataset_id"),
                "dataset_name": ds.get("title") or ds.get("name"),
                "cell_count": ds.get("cell_count"),
                "tissues": ";".join(str(t) for t in tissues),
                "assay": ";".join(str(a.get("label", a)) if isinstance(a, dict) else str(a)
                                  for a in (ds.get("assay") or [])),
                "disease": ";".join(str(d.get("label", d)) if isinstance(d, dict) else str(d)
                                    for d in (ds.get("disease") or [])),
                "matched_targets": ";".join(matched),
                "explorer_url": ds.get("explorer_url"),
            }
            rows.append(row)
            for k in matched:
                hits[k].append(row)
        if ci % 50 == 0:
            time.sleep(0.2)

    (RAW / "cellxgene_target_datasets.json").write_text(
        json.dumps(rows, indent=1), encoding="utf-8"
    )
    cols = [
        "matched_targets", "cell_count", "tissues", "dataset_name",
        "collection_name", "collection_doi", "dataset_id", "collection_id",
        "assay", "disease", "explorer_url",
    ]
    with (HERE / "cellxgene_target_datasets.tsv").open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["matched_targets"], -(x["cell_count"] or 0))):
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")

    summary = {k: len(v) for k, v in hits.items()}
    print("matched datasets per target:", json.dumps(summary, indent=1))
    (HERE / "cellxgene_summary.json").write_text(
        json.dumps({"n_collections": len(collections),
                    "n_target_datasets": len(rows),
                    "per_target": summary}, indent=1),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
