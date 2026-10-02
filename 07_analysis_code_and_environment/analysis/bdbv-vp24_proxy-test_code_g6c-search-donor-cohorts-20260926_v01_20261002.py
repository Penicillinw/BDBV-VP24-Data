"""G6c step 3: second, donor-cohort-oriented GEO screen.

Round 1 (tools/g6c_screen_primary_immune_20260926.py) found one design that can
test the G6 question in primary immune cells (GSE306664). A second, independent
design would let the answer be cross-checked, so this round targets datasets with
many HUMAN DONORS rather than many cell lines.
"""

from __future__ import annotations

import json
import os
import sys
import time
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926", "raw")

QUERIES = {
    "c01_ifna_healthy_donors": 'interferon alpha[Title] AND (healthy donor[Title] OR '
                               "healthy donors[Title] OR healthy volunteer[Title]) "
                               "AND Homo sapiens[Organism]",
    "c02_ifna_monocyte_donors": '"interferon alpha"[All Fields] AND monocyte[All Fields] '
                                "AND donor[Title] AND Homo sapiens[Organism] "
                                "AND expression profiling by high throughput sequencing[DataSet Type]",
    "c03_pbmc_ifn_response": 'interferon[Title] AND (PBMC[Title] OR PBMCs[Title]) AND '
                             "Homo sapiens[Organism]",
    "c04_interindividual_ifn": 'interferon[Title] AND (inter-individual[All Fields] OR '
                               "variability[Title] OR variability[All Fields]) "
                               "AND Homo sapiens[Organism]",
    "c05_blood_exvivo_ifn": 'interferon[Title] AND (whole blood[Title] OR leukocytes[Title] OR '
                            "blood[Title]) AND (ex vivo[All Fields] OR stimulation[Title]) "
                            "AND Homo sapiens[Organism]",
    "c06_ifnb_macrophage_donors": '(interferon beta[Title] OR IFN-beta[Title]) AND '
                                  "(macrophage*[Title] OR monocyte*[Title]) AND Homo sapiens[Organism]",
    "c07_isg_baseline_induction": '(baseline[Title] OR constitutive[Title]) AND '
                                  "(interferon[Title] OR ISG[Title]) AND Homo sapiens[Organism]",
    "c08_lupus_ifn_monocyte": 'interferon[Title] AND monocyte*[Title] AND (lupus[Title] OR SLE[Title]) '
                              "AND Homo sapiens[Organism]",
}


def get(url: str, tries: int = 4) -> str:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=120) as resp:  # noqa: S310
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    print(f"FAILED {url}: {last}")
    return ""


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(RAW, exist_ok=True)
    for tag, term in QUERIES.items():
        payload = get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gds&retmax=40&retmode=json&term={quote(term)}"
        )
        if not payload:
            continue
        result = json.loads(payload)["esearchresult"]
        ids = result.get("idlist", [])
        with open(os.path.join(RAW, f"{tag}_esearch.json"), "w", encoding="utf-8") as fh:
            json.dump({"term": term, "result": result}, fh, ensure_ascii=False, indent=1)
        if not ids:
            print(f"== {tag}: 0")
            continue
        summ = json.loads(
            get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=gds&retmode=json&id={','.join(ids)}"
            )
        )["result"]
        with open(os.path.join(RAW, f"{tag}_esummary.json"), "w", encoding="utf-8") as fh:
            json.dump(summ, fh, ensure_ascii=False, indent=1)
        print(f"== {tag}: {len(ids)}")
        for gid in ids:
            rec = summ.get(gid, {})
            if not rec.get("taxon", "").startswith("Homo"):
                continue
            print(
                f"   {rec.get('accession',''):>10} n={rec.get('n_samples',''):>4} "
                f"{rec.get('gdstype','')[:34]:<34} {rec.get('title','')[:88]}"
            )
        time.sleep(0.4)


if __name__ == "__main__":
    main()
