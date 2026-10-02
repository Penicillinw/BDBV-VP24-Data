"""GSE342661: primary differentiated human nasal epithelial ALI cultures + IFN-beta1 /
IFN-lambda1-4. Fetch the design and locate the count matrix.
"""

from __future__ import annotations

import os
import re
from urllib.request import Request, urlopen

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "g6g_gse342661_20260926")
RAW = os.path.join(OUT, "raw")
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE342nnn/GSE342661"


def get(url: str, tries: int = 3, binary: bool = False):
    last = None
    for a in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=300) as r:  # noqa: S310
                data = r.read()
            return data if binary else data.decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            import time

            time.sleep(2 * (a + 1))
    return None if binary else f"ERROR: {last}"


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    design = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                 "?acc=GSE342661&targ=self&form=text&view=quick")
    with open(os.path.join(RAW, "design.txt"), "w", encoding="utf-8") as fh:
        fh.write(design)
    for key in ["title", "overall_design", "type", "platform_id", "summary",
                "contributor", "pubmed_id", "relation"]:
        vals = re.findall(rf"!Series_{key} = (.+)", design)
        for v in vals[:3]:
            print(f"{key}: {v[:400]}")
    gsm = re.findall(r"!Series_sample_id = (GSM\d+)", design) or re.findall(r"=GSM\d+", design)
    print("GSM count in quick view:", len(set(gsm)))

    print("\n-- supplementary listing --")
    for path in ["suppl/", "matrix/"]:
        listing = get(f"{BASE}/{path}")
        files = re.findall(r'href="([^"?/][^"]*)"', listing or "")
        print(path, files)

    print("\n-- series matrix header (if present) --")
    raw = get(f"{BASE}/matrix/GSE342661_series_matrix.txt.gz", binary=True)
    if raw and isinstance(raw, bytes) and raw[:2] == b"\x1f\x8b":
        with open(os.path.join(RAW, "GSE342661_series_matrix.txt.gz"), "wb") as fh:
            fh.write(raw)
        import gzip

        with gzip.open(os.path.join(RAW, "GSE342661_series_matrix.txt.gz"), "rt",
                       encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(("!Sample_title", "!Sample_characteristics",
                                    "!Sample_treatment_protocol", "!Sample_source_name")):
                    print(line.rstrip()[:600])
    else:
        print("series matrix not downloadable")


if __name__ == "__main__":
    main()
