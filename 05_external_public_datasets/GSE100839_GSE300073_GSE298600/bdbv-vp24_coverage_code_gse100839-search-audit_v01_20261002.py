"""Audited public-resource search for the remaining gaps of the BDBV VP24 manuscript.

Gaps searched (each must end in a verdict: found / not found / found-but-unusable):
  D1  decisive experiment: same host cell type, two orthoebolaviruses, host transcriptome
  D2  low-end compartments: cell-type-resolved interferon data for urothelium, retinal
      pigment epithelium, photoreceptors, Sertoli, Leydig, choroid plexus, spleen
  D3  a third functional-anchor panel (viral readout + interferon + several cell types)
  D4  a primary-epithelium IFN-I/IFN-III panel with a released matrix (GSE342661 held)
  D5  paired mRNA-protein resources with wider protein coverage

Every query, its date, its hit count and the identifiers returned are written to
search_audit.tsv / search_audit.json, so the negative result is checkable rather than
asserted. Network calls are read-only; nothing is written outside this directory.
"""

from __future__ import annotations

import csv
import datetime
import json
import os
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = datetime.date.today().isoformat()
UA = {"User-Agent": "Mozilla/5.0 (research; manuscript gap audit)"}

GEO = [("D1", "gds", "Bundibugyo"),
       ("D1", "gds", "Bundibugyo AND Homo sapiens[Organism]"),
       ("D1", "gds", "ebolavirus AND transcriptome"),
       ("D1", "gds", "filovirus AND (RNA-seq OR \"single cell\")"),
       ("D1", "gds", "filovirus AND organoid"),
       ("D1", "sra", "Bundibugyo"),
       ("D1", "bioproject", "Bundibugyo"),
       ("D2", "gds", "spleen AND (single cell OR snRNA) AND Homo sapiens[Organism]"),
       ("D2", "gds", "(retina OR \"retinal pigment\") AND interferon AND Homo sapiens[Organism]"),
       ("D2", "gds", "(urothelium OR bladder) AND interferon AND Homo sapiens[Organism]"),
       ("D2", "gds", "(\"choroid plexus\" OR Sertoli OR Leydig) AND interferon AND Homo sapiens[Organism]"),
       ("D2", "gds", "\"ARPE-19\" AND (interferon OR IFN)"),
       ("D3", "gds", "(interferon OR IFN) AND (pseudotype OR \"viral load\" OR replication) AND (RNA-seq OR transcriptome)"),
       ("D4", "gds", "(\"airway epithelium\" OR \"nasal epithelium\" OR \"bronchial epithelium\") AND interferon AND (lambda OR \"type III\")"),
       ("D5", "gds", "proteome AND transcriptome AND (\"multiple tissues\" OR \"tissue atlas\") AND Homo sapiens[Organism]")]

EPMC = [("D1", 'Bundibugyo AND (transcriptom* OR "RNA-seq" OR sequencing)'),
        ("D1", '(Bundibugyo OR filovirus) AND (organoid OR "gut organoid") AND infection'),
        ("D1", 'filovirus AND "single-cell" AND (host response OR transcriptom*)'),
        ("D2", '("retinal pigment epithelium" OR ARPE-19) AND interferon AND transcriptome'),
        ("D2", 'urothelium AND interferon AND (transcriptome OR "RNA-seq")'),
        ("D3", 'interferon AND (filovirus OR Ebola) AND (pretreatment OR "dose-dependent") AND transcriptome'),
        ("D4", '("nasal epithelium" OR "airway epithelium") AND ("IFN-lambda" OR "type III interferon") AND transcriptome'),
        ("D5", 'paired proteome transcriptome atlas human tissues coverage')]


def get_json(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.load(fh)


def entrez(db, term, retmax=10):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=" + db
           + "&retmode=json&retmax=" + str(retmax) + "&term=" + urllib.parse.quote(term))
    d = get_json(url)
    res = d.get("esearchresult", {})
    return int(res.get("count", 0)), res.get("idlist", [])


def epmc(term, page=8):
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
           + urllib.parse.quote(term) + "&format=json&pageSize=" + str(page)
           + "&resultType=lite")
    d = get_json(url)
    hits = d.get("hitCount", 0)
    out = [{"id": r.get("pmid") or r.get("id"), "year": r.get("pubYear"),
            "journal": r.get("journalTitle", ""), "title": r.get("title", "")}
           for r in d.get("resultList", {}).get("result", [])]
    return hits, out


def main():
    rows = []
    detail = {}
    for gap, db, term in GEO:
        try:
            n, ids = entrez(db, term)
            rows.append({"gap": gap, "source": f"NCBI {db}", "query": term,
                         "date": DATE, "hits": n, "ids": ";".join(ids)})
            detail[f"{db}:{term}"] = ids
        except Exception as exc:                                  # network failure
            rows.append({"gap": gap, "source": f"NCBI {db}", "query": term,
                         "date": DATE, "hits": "ERROR", "ids": str(exc)[:120]})
    for gap, term in EPMC:
        try:
            n, recs = epmc(term)
            rows.append({"gap": gap, "source": "Europe PMC", "query": term,
                         "date": DATE, "hits": n,
                         "ids": ";".join(str(r["id"]) for r in recs)})
            detail[f"epmc:{term}"] = recs
        except Exception as exc:
            rows.append({"gap": gap, "source": "Europe PMC", "query": term,
                         "date": DATE, "hits": "ERROR", "ids": str(exc)[:120]})
    with open(os.path.join(HERE, "search_audit.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gap", "source", "query", "date", "hits", "ids"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "search_audit.json"), "w", encoding="utf-8") as fh:
        json.dump({"date": DATE, "rows": rows, "detail": detail}, fh,
                  ensure_ascii=False, indent=1)
    for r in rows:
        print(f"{r['gap']:3s} {r['source']:12s} hits={str(r['hits']):>6s}  {r['query'][:78]}")


if __name__ == "__main__":
    main()
