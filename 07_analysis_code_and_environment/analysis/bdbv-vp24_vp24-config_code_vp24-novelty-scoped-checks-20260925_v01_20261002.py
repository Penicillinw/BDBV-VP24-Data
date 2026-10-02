"""Decisive scoped checks: does a specific paper's indexed text contain the string?

Europe PMC indexes the full text of PMC author manuscripts even when the article
is not open access, so `EXT_ID:<pmid> AND "<string>"` answers directly whether a
paper wrote a residue out.  Results are saved as JSON next to the other raw
responses.
"""

import json
import os
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "vp24_novelty_lit_20260925")
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = {"User-Agent": "Mozilla/5.0 (research novelty audit; contact: local)"}

PMIDS = ("27974555", "37243162", "35069519", "37647113")
STRINGS = ("P83S", "S83", "Q135", "N135", "H140", "R140", "V141A", "A141",
           "four positions", "four residues", "native", "naturally",
           "conserved", "signature", "BDBV residues", "exclusive")


def hit_count(query):
    url = f"{BASE}?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "pageSize": 1}
    )
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=90) as resp:
        return json.loads(resp.read().decode("utf-8", "replace"))["hitCount"]


def main():
    table = {}
    for pmid in PMIDS:
        table[pmid] = {}
        for needle in STRINGS:
            query = f'EXT_ID:{pmid} AND "{needle}"'
            try:
                count = hit_count(query)
            except Exception as exc:  # noqa: BLE001
                count = f"error: {exc}"
            table[pmid][needle] = count
            print(f"{pmid:>9s}  {needle:<14s} -> {count}")
            time.sleep(0.4)
        print()
    with open(os.path.join(OUT, "raw", "scoped_string_checks.json"), "w",
              encoding="utf-8") as fh:
        json.dump(table, fh, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
