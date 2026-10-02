"""Direct test of the material question: immortalized lines versus primary human cells
in the same interferon-stimulation experiment (GEO GSE46599).

Design from the GEO record: 48 samples; 9 cell lines, primary CD4+ T cells and primary
macrophages; untreated and IFN-treated; 2 replicate experiments per cell line and 3 per
primary cell type. If the transformed lines are bimodal (some failing to respond at all)
while the primary cells respond, the cancer-cell-line material is demonstrably unsuitable.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import re
from collections import defaultdict
from urllib.request import Request, urlopen

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926")
RAW = os.path.join(OUT, "raw")
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46599"


def get(url: str, dest: str | None = None, tries: int = 4) -> bytes:
    if dest and os.path.exists(dest) and os.path.getsize(dest) > 0:
        return open(dest, "rb").read()
    last = None
    for a in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=600) as r:  # noqa: S310
                data = r.read()
            if dest:
                with open(dest, "wb") as fh:
                    fh.write(data)
            return data
        except Exception as exc:  # noqa: BLE001
            last = exc
            import time

            time.sleep(3 * (a + 1))
    raise SystemExit(f"fetch failed {url}: {last}")


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    # sample design first
    txt = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
              "?acc=GSE46599&targ=self&form=text&view=quick").decode("utf-8", "replace")
    with open(os.path.join(RAW, "GSE46599_design.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt)
    for key in ["title", "overall_design", "type", "platform_id", "summary"]:
        m = re.search(rf"!Series_{key} = (.+)", txt)
        if m:
            print(f"{key}: {m.group(1)[:300]}")
    gsms = re.findall(r"=GSM\d+", txt)
    print("GSM entries in design text:", len(gsms))

    # supplementary file listing
    supp = get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46599/suppl/").decode(
        "utf-8", "replace")
    files = re.findall(r'href="([^"?/][^"]*)"', supp)
    print("supplementary files:", files)

    # sample titles from the series matrix header
    try:
        raw = get(f"{BASE}/matrix/GSE46599_series_matrix.txt.gz",
                  os.path.join(RAW, "GSE46599_series_matrix.txt.gz"))
        with gzip.open(os.path.join(RAW, "GSE46599_series_matrix.txt.gz"), "rt",
                       encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("!Sample_title"):
                    titles = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
                    print("n titles:", len(titles))
                    for t in titles:
                        print("   ", t)
                    break
    except SystemExit as exc:
        print("series matrix unavailable:", exc)


if __name__ == "__main__":
    main()
