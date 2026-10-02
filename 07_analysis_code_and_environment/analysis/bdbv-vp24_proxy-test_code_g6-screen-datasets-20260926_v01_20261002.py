"""G6 step 1: screen public data for multi-cell-type / multi-cell-line interferon
stimulation experiments with matched controls.

Writes analysis/g6_ifn_responsiveness_20260926/datasets_screened.tsv and the raw
GEO responses.
"""

from __future__ import annotations

import csv
import json
import os
import time
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "g6_ifn_responsiveness_20260926")
RAW = os.path.join(OUT, "raw")

QUERIES = [
    '"interferon alpha"[Title] AND (cell line[Title] OR cell lines[Title]) '
    "AND Homo sapiens[Organism]",
    '"interferon beta"[Title] AND (cell line[Title] OR cell lines[Title]) '
    "AND Homo sapiens[Organism]",
    '"interferon lambda"[Title] OR "interferon-lambda"[Title] AND Homo sapiens[Organism]',
    '(IFN-alpha[Title] OR IFN-alpha[All Fields]) AND stimulation[All Fields] '
    "AND Homo sapiens[Organism] AND expression profiling by high throughput sequencing[DataSet Type]",
    'interferon[Title] AND (epithelial[Title] OR fibroblast[Title] OR endothelial[Title] '
    "OR hepatocyte[Title] OR monocyte[Title] OR macrophage[Title]) AND Homo sapiens[Organism]",
    '"interferon-stimulated genes"[Title] AND (panel[Title] OR atlas[Title] OR landscape[Title])',
    '"type I interferon"[Title] AND response[Title] AND human[Title]',
    'interferon[Title] AND (cell lines[Title] OR cell-line[Title]) AND treatment[All Fields]',
]


def get(url: str, tries: int = 4) -> str:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=120) as resp:  # noqa: S310 (fixed https host)
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - transient
            last = exc
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"request failed: {url}\n{last}")


def esearch(term: str, retmax: int = 40) -> dict:
    return json.loads(
        get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gds&retmax={retmax}&retmode=json&term={quote(term)}"
        )
    )["esearchresult"]


def esummary(ids: list[str]) -> dict:
    return json.loads(
        get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=gds&retmode=json&id={','.join(ids)}"
        )
    )["result"]


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    seen: dict[str, dict] = {}
    for term in QUERIES:
        try:
            res = esearch(term)
        except SystemExit as exc:
            print("query failed:", term, exc)
            continue
        ids = res.get("idlist", [])
        print(f"count={res.get('count')} ids={len(ids)} :: {term[:90]}")
        with open(os.path.join(RAW, f"esearch_{abs(hash(term)) % 10**8}.json"), "w", encoding="utf-8") as fh:
            json.dump(res, fh, ensure_ascii=False, indent=1)
        if not ids:
            continue
        summ = esummary(ids)
        with open(os.path.join(RAW, f"esummary_{abs(hash(term)) % 10**8}.json"), "w", encoding="utf-8") as fh:
            json.dump(summ, fh, ensure_ascii=False, indent=1)
        for uid in summ.get("uids", []):
            rec = summ[uid]
            acc = rec.get("accession", "")
            if not acc:
                continue
            if acc not in seen:
                seen[acc] = {
                    "accession": acc,
                    "title": rec.get("title", ""),
                    "n_samples": rec.get("n_samples", ""),
                    "gdstype": rec.get("gdstype", ""),
                    "taxon": rec.get("taxon", ""),
                    "summary": (rec.get("summary", "") or "")[:400],
                    "matched_queries": term[:60],
                }
        time.sleep(0.5)

    rows = sorted(seen.values(), key=lambda r: -int(r["n_samples"] or 0))
    with open(os.path.join(OUT, "datasets_screened.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n{len(rows)} unique series; top 25 by sample count:")
    for r in rows[:25]:
        print(f"  {r['accession']:12s} n={str(r['n_samples']):>4s}  {r['title'][:95]}")


if __name__ == "__main__":
    main()
