"""Probe GSE342661: primary differentiated human nasal epithelial cells (ALI) treated
with IFN-lambda1/2/3/4 or IFN-beta1 at 6/24/72 h.

Single cell type, but non-transformed differentiated epithelium and it contains the type
III arm. If it has several donors x conditions it can support a within-primary-cell-type
test of whether baseline tone predicts the response.
"""

from __future__ import annotations

import os
import re
from urllib.request import Request, urlopen

OUT = r"G:\本迪布焦研究\analysis\g6b_primary_epithelial_20260926\raw"
os.makedirs(OUT, exist_ok=True)


def get(url: str, tries: int = 4) -> str:
    last = None
    for a in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=120) as r:  # noqa: S310
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            import time

            time.sleep(2 * (a + 1))
    return f"ERROR {last}"


def main() -> None:
    text = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
               "?acc=GSE342661&targ=self&form=text&view=quick")
    with open(os.path.join(OUT, "GSE342661_quick.txt"), "w", encoding="utf-8") as fh:
        fh.write(text)
    for key in ["title", "overall_design", "type", "platform_id", "summary", "sample_id"]:
        m = re.search(rf"!Series_{key} = (.+)", text)
        if m:
            val = m.group(1)
            print(f"{key}: {val[:600]}")
    print("\nGSM count:", len(re.findall(r"=GSM\d+", text)))
    supp = get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE342nnn/GSE342661/suppl/")
    print("supplementary:", re.findall(r'href="([^"?/][^"]*)"', supp))


if __name__ == "__main__":
    main()
