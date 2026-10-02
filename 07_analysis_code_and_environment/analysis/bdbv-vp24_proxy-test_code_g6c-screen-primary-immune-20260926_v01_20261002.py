"""G6c step 1: screen public data for interferon stimulation of PRIMARY human
immune cells (monocyte / macrophage / dendritic cell / microglia) with matched
untreated controls in the same experiment.

Scope note: the main G6 line used GSE21158 (10 cancer cell lines). G6c asks the
same question - does baseline receptor/ISG expression proxy the magnitude of the
interferon response - in the cell types the manuscript actually ranks.

Outputs:
  analysis/g6c_primary_immune_20260926/datasets_screened.tsv
  analysis/g6c_primary_immune_20260926/raw/<tag>_esearch.json
  analysis/g6c_primary_immune_20260926/raw/<tag>_esummary.json
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926")
RAW = os.path.join(OUT, "raw")

# tag -> GEO (gds) query string
QUERIES = {
    "q01_mono_mac_ifn": '(interferon[Title] OR IFN[Title]) AND (monocyte*[Title] OR '
                        "macrophage*[Title] OR PBMC[Title]) AND Homo sapiens[Organism]",
    "q02_dc_microglia_ifn": '(interferon[Title] OR IFN[Title]) AND (dendritic[Title] OR '
                           "microglia*[Title]) AND Homo sapiens[Organism]",
    "q03_ifn_primary_all": 'interferon[Title] AND (primary[Title] OR monocyte-derived[Title]) '
                           "AND Homo sapiens[Organism]",
    "q04_ifnalpha_stim_seq": '("interferon alpha"[All Fields] OR "interferon-alpha"[All Fields] '
                             'OR IFNA[All Fields]) AND stimulation[All Fields] '
                             "AND Homo sapiens[Organism] "
                             "AND expression profiling by high throughput sequencing[DataSet Type]",
    "q05_macrophage_polarization": 'macrophage[Title] AND (polarization[Title] OR activation[Title]) '
                                   "AND interferon[All Fields] AND Homo sapiens[Organism]",
    "q06_sc_ifn_pbmc": '("single cell"[All Fields] OR scRNA[All Fields]) AND interferon[All Fields] '
                       "AND (PBMC[All Fields] OR monocyte[All Fields] OR macrophage[All Fields]) "
                       "AND Homo sapiens[Organism]",
    "q07_isg_primary": '"interferon-stimulated genes"[All Fields] AND (monocyte*[Title] OR '
                       "macrophage*[Title] OR dendritic[Title]) AND Homo sapiens[Organism]",
    "q08_ifn_dose_time": 'interferon[Title] AND (dose[Title] OR time course[Title] OR kinetics[Title]) '
                         "AND Homo sapiens[Organism]",
    "q09_whole_blood_ifn": 'interferon[Title] AND (blood[Title] OR leukocyte*[Title]) '
                           "AND Homo sapiens[Organism]",
    "q10_monocyte_donors_ifn": '(monocyte*[Title] OR macrophage*[Title]) AND interferon[All Fields] '
                               "AND donor*[All Fields] AND Homo sapiens[Organism] "
                               "AND expression profiling by high throughput sequencing[DataSet Type]",
}

RETMAX = 50


def get(url: str, tries: int = 4) -> str:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=120) as resp:  # noqa: S310 (fixed https host)
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - transient network failures
            last = exc
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"request failed: {url}\n{last}")


def esearch(term: str) -> dict:
    payload = get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=gds&retmax={RETMAX}&retmode=json&term={quote(term)}"
    )
    return json.loads(payload)["esearchresult"]


def esummary(ids: list[str]) -> dict:
    payload = get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=gds&retmode=json&id={','.join(ids)}"
    )
    return json.loads(payload)["result"]


def main() -> None:
    try:  # console may be a legacy code page; never let a printed title kill the screen
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(RAW, exist_ok=True)
    rows: dict[str, dict] = {}
    for tag, term in QUERIES.items():
        result = esearch(term)
        ids = result.get("idlist", [])
        with open(os.path.join(RAW, f"{tag}_esearch.json"), "w", encoding="utf-8") as fh:
            json.dump({"term": term, "result": result}, fh, ensure_ascii=False, indent=1)
        if not ids:
            print(f"{tag}: 0 hits")
            continue
        summ = esummary(ids)
        with open(os.path.join(RAW, f"{tag}_esummary.json"), "w", encoding="utf-8") as fh:
            json.dump(summ, fh, ensure_ascii=False, indent=1)
        for gid in ids:
            rec = summ.get(gid, {})
            if not rec:
                continue
            key = rec.get("accession", gid)
            row = rows.setdefault(
                key,
                {
                    "accession": key,
                    "gds_id": gid,
                    "title": (rec.get("title") or "").replace("\t", " "),
                    "taxon": rec.get("taxon", ""),
                    "gdstype": rec.get("gdstype", ""),
                    "n_samples": rec.get("n_samples", ""),
                    "pdat": rec.get("pdat", ""),
                    "queries": [],
                    "summary": (rec.get("summary") or "").replace("\t", " ")[:1200],
                },
            )
            row["queries"].append(tag)
            print(f"{tag}: {key} | {rec.get('taxon','')} | {rec.get('n_samples','')} samples | {rec.get('title','')[:90]}")
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "datasets_screened.tsv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "accession",
                "gds_id",
                "title",
                "taxon",
                "gdstype",
                "n_samples",
                "pdat",
                "queries",
                "summary",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in sorted(rows.values(), key=lambda r: r["accession"]):
            row = dict(row)
            row["queries"] = ",".join(row["queries"])
            writer.writerow(row)
    print(f"\n{len(rows)} unique series -> {path}")


if __name__ == "__main__":
    main()
