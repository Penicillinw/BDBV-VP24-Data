"""Confirm that the GSE342661 count matrix is really downloadable (read-only)."""

from __future__ import annotations

import os
import re
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
os.makedirs(RAW, exist_ok=True)

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE342nnn/GSE342661/suppl/"
FILES = ["GSE342661_counts_tpm.matrix.gz"]


def head(url):
    req = urllib.request.Request(url, method="HEAD",
                                 headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return fh.status, fh.headers.get("Content-Length")


for name in FILES:
    try:
        status, size = head(BASE + name)
        print(f"{name}: HTTP {status}, {size} bytes")
    except Exception as exc:                     # noqa: BLE001
        print(f"{name}: FAILED -> {exc}")

try:
    req = urllib.request.Request(BASE, headers={"User-Agent": "Mozilla/5.0"})
    html = urllib.request.urlopen(req, timeout=60).read().decode("utf-8", "replace")
    with open(os.path.join(RAW, "suppl_listing.html"), "w", encoding="utf-8") as fh:
        fh.write(html)
    listing = re.findall(r'href="([^"]+)"', html)
    print("listing:", [x for x in listing if x.lower().endswith(
        (".gz", ".txt", ".csv", ".tsv", ".xlsx", ".tar"))][:12])
except Exception as exc:                         # noqa: BLE001
    print("listing failed:", exc)
