"""Exhaustive evidence trail for one question:

  Does ANY public dataset contain BDBV and EBOV (or other orthoebolaviruses)
  side by side in the same cell type / same experimental system AND profile the
  HOST transcriptome?

Every query, hit count and per-record verdict is written to
analysis/t4_host_tx_hunt_20260925/ so that a negative answer is auditable.
"""

import json
import os
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "t4_host_tx_hunt_20260925")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def fetch(url, tries=3, timeout=90):
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2)
    raise RuntimeError(f"failed: {url} ({last})")


def esearch(db, term, retmax=40):
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmax": retmax, "retmode": "json"})
    return json.loads(fetch(url))["esearchresult"]


def esummary(db, ids):
    if not ids:
        return {}
    url = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": db, "id": ",".join(ids), "retmode": "json"})
    return json.loads(fetch(url))["result"]


def epmc(query, page_size=25):
    url = EPMC + "?" + urllib.parse.urlencode(
        {"query": query, "resultType": "core", "format": "json",
         "pageSize": page_size})
    return json.loads(fetch(url))


GEO_QUERIES = [
    ("geo_bundibugyo", "Bundibugyo"),
    ("geo_bundibugyo_rnaseq", "Bundibugyo AND (RNA-seq OR transcriptome OR expression profiling)"),
    ("geo_bundibugyo_celltype", "Bundibugyo AND (Huh7 OR hepatocyte OR macrophage OR monocyte OR dendritic OR organoid)"),
    ("geo_bundibugyo_human", "Bundibugyo AND Homo sapiens[Organism]"),
    ("geo_ebov_and_bdbv", "Ebola AND Bundibugyo"),
    ("geo_filovirus_comparative", "(Ebola OR Bundibugyo) AND (comparative OR dual OR four viruses)"),
    ("sra_bundibugyo", "Bundibugyo"),
    ("bioproject_bundibugyo", "Bundibugyo"),
]

EPMC_QUERIES = [
    ("epmc_bdbv_transcriptome", '(TITLE_ABS:"Bundibugyo" AND (TITLE_ABS:"transcriptome" OR TITLE_ABS:"RNA-seq" OR TITLE_ABS:"gene expression"))'),
    ("epmc_bdbv_vs_ebov", '(TITLE_ABS:"Bundibugyo" AND TITLE_ABS:"Ebola")'),
    ("epmc_bdbv_macrophage", '(TITLE_ABS:"Bundibugyo" AND (TITLE_ABS:"macrophage" OR TITLE_ABS:"monocyte" OR TITLE_ABS:"dendritic"))'),
]


def main():
    os.makedirs(OUT, exist_ok=True)
    report = {"geo": {}, "epmc": {}, "notes": []}

    for tag, term in GEO_QUERIES:
        db = "sra" if tag.startswith("sra") else "bioproject" if tag.startswith("bioproject") else "gds"
        try:
            res = esearch(db, term, retmax=30)
        except Exception as exc:  # noqa: BLE001
            report["geo"][tag] = {"db": db, "term": term, "error": str(exc)}
            continue
        ids = res.get("idlist", [])
        entries = []
        if ids:
            summ = esummary(db, ids)
            for uid in summ.get("uids", []):
                rec = summ[uid]
                entries.append({
                    "uid": uid,
                    "accession": rec.get("accession") or rec.get("project_acc"),
                    "title": (rec.get("title") or "")[:180],
                    "n_samples": rec.get("n_samples"),
                    "gdstype": rec.get("gdstype"),
                    "taxon": rec.get("taxon") or rec.get("organism"),
                    "summary": (rec.get("summary") or "")[:400],
                })
        report["geo"][tag] = {"db": db, "term": term,
                              "count": int(res.get("count", 0)), "entries": entries}
        print(f"[{db}] {tag}: {res.get('count')} hits")

    for tag, query in EPMC_QUERIES:
        try:
            payload = epmc(query)
        except Exception as exc:  # noqa: BLE001
            report["epmc"][tag] = {"query": query, "error": str(exc)}
            continue
        results = payload.get("resultList", {}).get("result", [])
        report["epmc"][tag] = {
            "query": query,
            "hitCount": payload.get("hitCount"),
            "top": [{"pmid": r.get("pmid"), "year": r.get("pubYear"),
                     "journal": (r.get("journalInfo", {}).get("journal", {}) or {}).get("title"),
                     "title": (r.get("title") or "")[:180]} for r in results],
        }
        print(f"[epmc] {tag}: {payload.get('hitCount')} hits")

    json.dump(report, open(os.path.join(OUT, "search_evidence.json"), "w",
                           encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\nwritten", os.path.join(OUT, "search_evidence.json"))


if __name__ == "__main__":
    main()
