"""G6 step 3: download GSE21158 (10 human cancer cell lines +/- IFN-alpha-2a).

Series design (from the GEO record): a panel of 10 human cancer cell lines
(ME-15, D10, LS174T, AsPC1, A549, MIAPACA, CaLu6, HCT116, PANC1, JUSO) treated
with 100 U/ml IFN-alpha-2a; control (C), 4 h (H) and 24 h (D), triplicates,
Illumina microarray. This gives baseline and induction in the same platform for
the same cell line - the cleanest available test of whether baseline expression
proxies interferon responsiveness.
"""

from __future__ import annotations

import gzip
import os
import re
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "g6_ifn_responsiveness_20260926")
RAW = os.path.join(OUT, "raw")
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21158"


def download(url: str, dest: str) -> str:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print("cached", dest)
        return dest
    req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
    with urlopen(req, timeout=300) as resp, open(dest, "wb") as fh:  # noqa: S310
        fh.write(resp.read())
    print("downloaded", dest, os.path.getsize(dest), "bytes")
    return dest


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    matrix = download(
        f"{BASE}/matrix/GSE21158_series_matrix.txt.gz",
        os.path.join(RAW, "GSE21158_series_matrix.txt.gz"),
    )
    with gzip.open(matrix, "rt", encoding="utf-8", errors="replace") as fh:
        head: list[str] = []
        for line in fh:
            head.append(line.rstrip("\n"))
            if line.startswith("!series_matrix_table_begin"):
                break
    titles = [
        re.sub(r'^"|"$', "", re.sub(r"^!Sample_title\t", "", l))
        for l in head
        if l.startswith("!Sample_title")
    ]
    platform = [l for l in head if l.startswith("!Sample_platform_id")]
    accession = [l for l in head if l.startswith("!Series_geo_accession")]
    print("series:", accession[0] if accession else "?")
    print("platform:", platform[0] if platform else "?")
    print("n sample titles:", len(titles))
    print("titles:", titles[:12])


if __name__ == "__main__":
    main()
