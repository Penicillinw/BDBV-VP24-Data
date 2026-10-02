"""Read the GEO design records of the shortlisted primary-cell / organoid IFN datasets."""

from __future__ import annotations

import os
import re
import time
from urllib.request import Request, urlopen

OUT = r"G:\本迪布焦研究\analysis\g6b_primary_epithelial_20260926"
RAW = os.path.join(OUT, "raw")
CANDIDATES = ["GSE342661", "GSE46599", "GSE305982", "GSE148712", "GSE145496",
              "GSE208684", "GSE243230", "GSE206680"]


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
    return ""


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    for acc in CANDIDATES:
        text = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                   f"?acc={acc}&targ=self&form=text&view=quick")
        if not text:
            print(f"== {acc}: fetch failed")
            continue
        with open(os.path.join(RAW, f"{acc}_design.txt"), "w", encoding="utf-8") as fh:
            fh.write(text)
        def field(name: str) -> str:
            m = re.search(rf"!Series_{name} = (.+)", text)
            return m.group(1).strip() if m else ""

        print(f"\n===== {acc}: {field('title')}")
        print("  type:", field("type"), "| platform:", field("platform_id")[:40])
        print("  samples:", field("sample_id").count("GSM"))
        design = field("overall_design")
        print("  design:", design[:700])
        time.sleep(0.4)


if __name__ == "__main__":
    main()
