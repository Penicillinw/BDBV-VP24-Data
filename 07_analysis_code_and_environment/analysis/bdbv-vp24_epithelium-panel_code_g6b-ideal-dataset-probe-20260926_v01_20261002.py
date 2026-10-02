"""Targeted probe for the ideal G6 dataset: several NON-transformed human cell types or
organoids, interferon stimulation, matched controls, in one experiment."""

from __future__ import annotations

import json
import os
import time
from urllib.parse import quote
from urllib.request import Request, urlopen

OUT = r"G:\本迪布焦研究\analysis\g6b_primary_epithelial_20260926"
RAW = os.path.join(OUT, "raw")

QUERIES = [
    '(interferon OR IFN) AND organoid AND Homo sapiens[Organism] AND expression profiling[DataSet Type]',
    '(interferon OR IFN) AND (primary cells[Title] OR primary human cells[Title]) AND Homo sapiens[Organism]',
    '"interferon lambda" OR "IFN-lambda" AND (intestinal[Title] OR airway[Title] OR epithelial[Title])',
    'interferon AND (enteroid OR colonoid OR organoids[Title]) AND treated[All Fields]',
    '(interferon OR IFN) AND (keratinocytes[Title] AND fibroblasts[Title]) AND Homo sapiens[Organism]',
    '(interferon OR IFN) AND (airway[Title] AND epithelial[Title]) AND (ALI[All Fields] OR differentiated[All Fields])',
    '(interferon OR IFN) AND (monocyte-derived[Title] OR macrophages[Title] OR dendritic cells[Title]) AND Homo sapiens[Organism]',
    'interferon AND stimulation AND (single cell[All Fields] OR scRNA[All Fields]) AND Homo sapiens[Organism]',
]


def get(url: str, tries: int = 4) -> str:
    last = None
    for a in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=120) as r:  # noqa: S310
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (a + 1))
    raise SystemExit(str(last))


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    seen: dict[str, dict] = {}
    for term in QUERIES:
        res = json.loads(
            get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                f"?db=gds&retmax=30&retmode=json&term={quote(term)}")
        )["esearchresult"]
        ids = res.get("idlist", [])
        if not ids:
            print(f"0 hits :: {term[:70]}")
            continue
        summ = json.loads(
            get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=gds&retmode=json&id={','.join(ids)}")
        )["result"]
        for uid in summ.get("uids", []):
            r = summ[uid]
            acc = r.get("accession", "")
            if acc.startswith("GSE") and acc not in seen:
                seen[acc] = {
                    "accession": acc,
                    "n_samples": r.get("n_samples", ""),
                    "title": r.get("title", ""),
                    "gdstype": r.get("gdstype", ""),
                    "summary": (r.get("summary", "") or "")[:300],
                }
        print(f"{res.get('count')} hits :: {term[:70]}")
        time.sleep(0.5)

    rows = sorted(seen.values(), key=lambda x: -int(x["n_samples"] or 0))
    with open(os.path.join(OUT, "ideal_probe.tsv"), "w", encoding="utf-8", newline="") as fh:
        import csv

        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} unique series")
    for r in rows[:20]:
        print(f"  {r['accession']:12s} n={str(r['n_samples']):>4s}  {r['title'][:92]}")


if __name__ == "__main__":
    main()
