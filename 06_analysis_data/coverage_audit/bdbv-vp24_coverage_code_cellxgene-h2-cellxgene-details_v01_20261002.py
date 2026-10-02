"""H2 line: fetch CellxGene collection details for collections whose tissue
labels match the six low-end compartments, then test each dataset's cell_type
vocabulary for the specific cell types of interest.

Read-only outside this directory.
"""
from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
RAW.mkdir(parents=True, exist_ok=True)
API = "https://api.cellxgene.cziscience.com/curation/v1"
UA = {"User-Agent": "Mozilla/5.0 (research; BDBV-VP24 atlascan)"}

TISSUE_RULES = {
    "urothelium": r"bladder|ureter|renal pelvis|urothel",
    "eye": r"retina|fovea|macula|\beye\b|iris|cornea|choriocapillaris|optic",
    "testis": r"testis|testicle|seminiferous|epididymis",
    "choroid_plexus": r"choroid plexus",
    "spleen": r"spleen|splenic",
}

# cell_type label patterns that count as the compartment of interest
CELLTYPE_RULES = {
    "urothelium": re.compile(r"urothel|umbrella cell|bladder epithel", re.I),
    "RPE": re.compile(r"retinal pigment|(^|\b)rpe(\b|$)", re.I),
    "photoreceptor": re.compile(r"photoreceptor|\brod\b|\bcone\b|cone photoreceptor", re.I),
    "Sertoli": re.compile(r"sertoli", re.I),
    "choroid_plexus_epithelium": re.compile(r"choroid plexus", re.I),
    "spleen_pulp": re.compile(r"splen|red pulp|white pulp|marginal zone", re.I),
}


def get(url: str, timeout: int = 90):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def main() -> None:
    colls = json.loads((RAW / "cellxgene_collections.json").read_text(encoding="utf-8"))

    want: dict[str, set[str]] = {}
    for c in colls:
        for ds in c.get("datasets") or []:
            taxa = [t.get("label") for t in (ds.get("organism") or [])]
            if "Homo sapiens" not in taxa:
                continue
            labels = " ; ".join(
                str(t.get("label", "")) for t in (ds.get("tissue") or [])
            )
            for tgt, pat in TISSUE_RULES.items():
                if re.search(pat, labels, re.I):
                    want.setdefault(c.get("collection_id"), set()).add(tgt)

    print(f"collections with matching tissue labels: {len(want)}")

    detail_dir = RAW / "cellxgene_collection_details"
    detail_dir.mkdir(exist_ok=True)

    n_fetched = 0
    for cid in sorted(want):
        fp = detail_dir / f"{cid}.json"
        if fp.exists():
            continue
        try:
            j = get(f"{API}/collections/{cid}")
            fp.write_text(json.dumps(j, indent=1), encoding="utf-8")
            n_fetched += 1
        except urllib.error.HTTPError as exc:
            (detail_dir / f"{cid}.error.txt").write_text(str(exc), encoding="utf-8")
        except Exception as exc:  # noqa: BLE001
            (detail_dir / f"{cid}.error.txt").write_text(repr(exc), encoding="utf-8")
        time.sleep(0.25)
    print(f"fetched {n_fetched} new collection details")

    rows = []
    for cid in sorted(want):
        fp = detail_dir / f"{cid}.json"
        if not fp.exists():
            continue
        j = json.loads(fp.read_text(encoding="utf-8"))
        for ds in j.get("datasets") or []:
            taxa = [t.get("label") for t in (ds.get("organism") or [])]
            if "Homo sapiens" not in taxa:
                continue
            if ds.get("tombstone"):
                continue
            cts = [c.get("label", "") for c in (ds.get("cell_type") or [])]
            hit_types = sorted(
                {k for k, pat in CELLTYPE_RULES.items() if any(pat.search(c) for c in cts)}
            )
            tissue_labels = ";".join(str(t.get("label", "")) for t in (ds.get("tissue") or []))
            tissue_hits = sorted(
                {k for k, pat in TISSUE_RULES.items() if re.search(pat, tissue_labels, re.I)}
            )
            assets = ds.get("assets") or []
            h5ad = [a for a in assets if (a.get("filetype") or "").upper() == "H5AD"]
            rds = [a for a in assets if (a.get("filetype") or "").upper() == "RDS"]
            rows.append({
                "collection_id": cid,
                "collection_name": j.get("name"),
                "collection_doi": j.get("doi"),
                "dataset_id": ds.get("dataset_id"),
                "dataset_title": ds.get("title"),
                "cell_count": ds.get("cell_count"),
                "primary_cell_count": ds.get("primary_cell_count"),
                "is_primary_data": ds.get("is_primary_data"),
                "tissue": tissue_labels,
                "tissue_hits": ";".join(tissue_hits),
                "celltype_hits": ";".join(hit_types),
                "n_cell_types": len(cts),
                "cell_type_labels": "|".join(cts),
                "assay": ";".join(a.get("label", "") for a in (ds.get("assay") or [])),
                "suspension_type": ";".join(ds.get("suspension_type") or []),
                "disease": ";".join(d.get("label", "") for d in (ds.get("disease") or [])),
                "h5ad_url": h5ad[0]["url"] if h5ad else "",
                "h5ad_size": h5ad[0].get("filesize") if h5ad else "",
                "rds_url": rds[0]["url"] if rds else "",
                "explorer_url": ds.get("explorer_url"),
            })

    cols = list(rows[0].keys()) if rows else []
    with (HERE / "cellxgene_dataset_details.tsv").open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["celltype_hits"] == "", -(x["cell_count"] or 0))):
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")

    withhit = [r for r in rows if r["celltype_hits"]]
    print(f"human datasets scanned: {len(rows)}; with target cell types: {len(withhit)}")
    per = {}
    for r in withhit:
        for k in r["celltype_hits"].split(";"):
            per.setdefault(k, []).append((r["dataset_id"], r["cell_count"], r["collection_name"]))
    for k, v in sorted(per.items()):
        print(f"  {k}: {len(v)} datasets")
    (HERE / "cellxgene_celltype_hits.json").write_text(
        json.dumps({k: v for k, v in per.items()}, indent=1), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
