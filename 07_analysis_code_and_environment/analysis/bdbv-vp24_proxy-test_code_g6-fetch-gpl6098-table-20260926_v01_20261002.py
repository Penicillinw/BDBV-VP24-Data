"""Download the complete GPL6098 platform table (probe ID -> Symbol) and check the
series-matrix probe ID format."""

from __future__ import annotations

import gzip
import os
import re
from urllib.request import Request, urlopen

ROOT = r"G:\本迪布焦研究"
RAW = os.path.join(ROOT, "analysis", "g6_ifn_responsiveness_20260926", "raw")
DEST = os.path.join(RAW, "GPL6098_platform_table.txt")

URL = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
       "?acc=GPL6098&targ=self&form=text&view=full")


def main() -> None:
    if not (os.path.exists(DEST) and os.path.getsize(DEST) > 2_000_000):
        req = Request(URL, headers={"User-Agent": "CodexResearch/1.0"})
        with urlopen(req, timeout=600) as resp, open(DEST, "wb") as fh:  # noqa: S310
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
    text = open(DEST, encoding="utf-8", errors="replace").read()
    print("platform table chars:", len(text), "complete:", "!platform_table_end" in text)
    header = None
    n = 0
    for line in text.splitlines():
        if line.startswith("!platform_table_begin"):
            continue
        if header is None and line.startswith("ID\t"):
            header = line.split("\t")
            print("columns:", header)
            continue
        if header and line and not line.startswith("!"):
            n += 1
    print("table rows:", n)

    with gzip.open(os.path.join(RAW, "GSE21158_series_matrix.txt.gz"), "rt",
                   encoding="utf-8", errors="replace") as fh:
        in_table = False
        probes = []
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                next(fh)
                continue
            if in_table:
                if line.startswith("!series_matrix_table_end"):
                    break
                probes.append(line.split("\t")[0].strip('"'))
                if len(probes) >= 12:
                    break
    print("series matrix first probes:", probes)


if __name__ == "__main__":
    main()
