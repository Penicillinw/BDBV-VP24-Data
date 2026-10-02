"""Main-agent check: release status and contents of GEO GSE342661.

Read-only.  Raw responses are dumped next to this file under raw/.  Criterion
frozen before looking at the answer:

    * USABLE only if a sample-level expression matrix is downloadable without
      contacting the authors AND the sample metadata maps every column back to
      a treatment condition;
    * otherwise USABLE-AFTER-REQUANT if raw reads exist, else UNAVAILABLE.
"""

from __future__ import annotations

import os
import re
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
os.makedirs(RAW, exist_ok=True)
ACC = "GSE342661"


def get(url):
    req = urllib.request.Request(url, headers={
        "User-Agent": "Mozilla/5.0 (research metadata check)",
        "Accept": "text/html,application/xml,text/plain,*/*"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return fh.read()


def save(name, data):
    path = os.path.join(RAW, name)
    with open(path, "wb") as fh:
        fh.write(data)
    return path


def main():
    text = ""
    try:
        data = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
                   + ACC + "&targ=self&form=text&view=brief")
        save(ACC + "_soft_brief.txt", data)
        text = data.decode("utf-8", "replace")
        print("=== GEO SOFT brief ===")
        for ln in text.splitlines():
            if re.match(r"^\s*!Series_(title|status|submission_date|"
                        r"last_update_date|pubmed_id|summary|overall_design|"
                        r"type|platform_id|relation|sample_id)", ln):
                print(ln[:300])
    except Exception as exc:                     # noqa: BLE001
        print("!! accession request failed:", exc)

    print("sample_id lines:",
          len([l for l in text.splitlines() if l.startswith("!Series_sample_id")]))

    try:
        data = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + ACC)
        save(ACC + "_acc_page.html", data)
        html = data.decode("utf-8", "replace")
        print("\n=== supplementary links ===")
        found = sorted(set(re.findall(
            r'href="([^"]*(?:suppl|GSE342661)[^"]*)"', html)))
        for f in found[:40]:
            print(" ", f)
        if not found:
            print("  (none found)")
        for key in ("not yet available", "private", "embargo", "held"):
            m = re.search(key, html, re.I)
            if m:
                s = max(0, m.start() - 120)
                print("  [%s] ...%s..." % (key, html[s:m.start() + 160]))
                break
    except Exception as exc:                     # noqa: BLE001
        print("!! accession page request failed:", exc)

    try:
        es = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                 "?db=gds&term=" + ACC + "&retmode=json").decode()
        save(ACC + "_esearch.json", es.encode())
        print("\n=== esearch ===")
        print(es.strip()[:400])
    except Exception as exc:                     # noqa: BLE001
        print("!! esearch failed:", exc)


if __name__ == "__main__":
    main()
