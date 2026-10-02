"""H2 line: GEO (E-utilities) and Europe PMC queries for human single-cell
atlases covering the six low-end compartments. Raw responses are stored.
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
RAW.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (research; BDBV-VP24 lowend-atlas-scan)"}

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

GEO_QUERIES = [
    ('bladder_singlecell', 'human bladder AND (single cell OR single-cell OR scRNA) AND (urothelium OR urothelial)'),
    ('bladder_atlas', 'bladder urothelium atlas AND "Homo sapiens"[Organism]'),
    ('ureter_singlecell', 'ureter OR "renal pelvis" AND (single cell RNA) AND "Homo sapiens"[Organism]'),
    ('rpe_singlecell', '"retinal pigment epithelium" AND (single cell OR scRNA OR "single-cell RNA")'),
    ('retina_atlas', 'human retina single cell atlas AND "Homo sapiens"[Organism]'),
    ('photoreceptor_singlecell', 'photoreceptor AND (single cell RNA) AND "Homo sapiens"[Organism]'),
    ('testis_singlecell', 'human testis single cell AND (Sertoli OR spermatogonia)'),
    ('sertoli_singlecell', 'Sertoli cell AND (single cell RNA) AND "Homo sapiens"[Organism]'),
    ('choroid_plexus_singlecell', '"choroid plexus" AND (single cell OR single-nucleus) AND human'),
    ('spleen_singlecell', 'human spleen AND (single cell OR single-nucleus) AND atlas'),
    ('tabula_sapiens', '"Tabula Sapiens"'),
    ('pan_organ_atlas', 'human cell atlas AND (multi-organ OR pan-tissue OR 24 organs) AND single cell'),
]

EPMC_QUERIES = [
    ('epmc_bladder_atlas', '(TITLE_ABS:"bladder" AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"single cell"))'),
    ('epmc_rpe_atlas', '(TITLE_ABS:"retinal pigment epithelium" AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"single cell" OR TITLE_ABS:"atlas"))'),
    ('epmc_retina_atlas', '(TITLE_ABS:"human retina" AND TITLE_ABS:"atlas")'),
    ('epmc_testis_atlas', '(TITLE_ABS:"human testis" AND TITLE_ABS:"single-cell")'),
    ('epmc_sertoli', '(TITLE_ABS:"Sertoli" AND TITLE_ABS:"human")'),
    ('epmc_choroid_plexus', '(TITLE_ABS:"choroid plexus" AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"single cell" OR TITLE_ABS:"transcriptomic"))'),
    ('epmc_spleen_atlas', '(TITLE_ABS:"human spleen" AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"atlas"))'),
    ('epmc_pan_organ', '(TITLE_ABS:"Tabula Sapiens" OR TITLE_ABS:"human cell atlas")'),
]


def fetch(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def geo_search(term: str, retmax: int = 20):
    q = urllib.parse.urlencode({
        "db": "gds", "term": term, "retmax": retmax, "retmode": "json",
        "sort": "relevance",
    })
    js = json.loads(fetch(f"{EUTILS}/esearch.fcgi?{q}"))
    ids = js["esearchresult"].get("idlist", [])
    out = {"query": term, "count": int(js["esearchresult"].get("count", 0)), "entries": []}
    if not ids:
        return out
    q2 = urllib.parse.urlencode({"db": "gds", "id": ",".join(ids), "retmode": "json"})
    sj = json.loads(fetch(f"{EUTILS}/esummary.fcgi?{q2}"))
    for uid in sj.get("result", {}).get("uids", []):
        rec = sj["result"][uid]
        out["entries"].append({
            "uid": uid,
            "accession": rec.get("accession"),
            "title": rec.get("title"),
            "n_samples": rec.get("n_samples"),
            "gdstype": rec.get("gdstype"),
            "taxon": rec.get("taxon"),
            "summary": (rec.get("summary") or "")[:400],
        })
    return out


def epmc_search(query: str, page_size: int = 15):
    q = urllib.parse.urlencode({
        "query": query, "format": "json", "pageSize": page_size, "resultType": "core",
        "sort": "CITED desc",
    })
    js = json.loads(fetch(f"{EPMC}?{q}"))
    res = js.get("resultList", {}).get("result", [])
    return {
        "query": query,
        "hitCount": js.get("hitCount"),
        "top": [
            {
                "pmid": r.get("pmid"), "pmcid": r.get("pmcid"),
                "doi": r.get("doi"), "year": r.get("pubYear"),
                "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title"),
                "title": r.get("title"),
                "is_preprint": r.get("source") in {"PPR"},
                "source": r.get("source"),
            }
            for r in res
        ],
    }


def main() -> None:
    out = {"geo": {}, "epmc": {}}
    for name, term in GEO_QUERIES:
        try:
            out["geo"][name] = geo_search(term)
            print(f"GEO {name}: count={out['geo'][name]['count']}")
        except Exception as exc:  # noqa: BLE001
            out["geo"][name] = {"query": term, "error": repr(exc)}
            print(f"GEO {name}: ERROR {exc}")
        time.sleep(0.5)
    for name, q in EPMC_QUERIES:
        try:
            out["epmc"][name] = epmc_search(q)
            print(f"EPMC {name}: hits={out['epmc'][name]['hitCount']}")
        except Exception as exc:  # noqa: BLE001
            out["epmc"][name] = {"query": q, "error": repr(exc)}
            print(f"EPMC {name}: ERROR {exc}")
        time.sleep(0.5)
    (RAW / "geo_epmc_scan.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
