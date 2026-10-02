"""Locate the actual data for GSE342661 (180 samples): per-GSM supplementary links,
SRA/BioProject route, and the GEO FTP mirrors."""

from __future__ import annotations

import os
import re
from urllib.request import Request, urlopen

RAW = r"G:\本迪布焦研究\analysis\g6g_gse342661_20260926\raw"


def get(url: str, tries: int = 3) -> str:
    last = None
    for a in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=180) as r:  # noqa: S310
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            import time

            time.sleep(2 * (a + 1))
    return f"ERROR: {last}"


def main() -> None:
    design = open(os.path.join(RAW, "design.txt"), encoding="utf-8").read()
    gsms = sorted(set(re.findall(r"(GSM\d+)", design)))
    print("GSMs found:", len(gsms), gsms[:5], "...", gsms[-3:] if gsms else "")

    for gsm in gsms[:3]:
        txt = get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
                  "&targ=self&form=text&view=brief")
        print(f"\n== {gsm}")
        for key in ["title", "supplementary_file", "source_name", "characteristics",
                    "treatment_protocol_ch1", "extraction_protocol"]:
            for v in re.findall(rf"!Sample_{key} = (.+)", txt)[:3]:
                print(f"   {key}: {v.strip()[:200]}")
        rel = re.findall(r"!Sample_relation = (.+)", txt)
        for v in rel[:3]:
            print("   relation:", v.strip()[:160])

    print("\n-- BioProject / SRA route --")
    for url in ["https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1508721",
                "https://www.ebi.ac.uk/ena/browser/api/xml/PRJNA1508721"]:
        t = get(url)
        print(url, "->", ("len %d" % len(t)) if not t.startswith("ERROR") else t[:80])
        if url.endswith("xml") and not t.startswith("ERROR"):
            with open(os.path.join(RAW, "PRJNA1508721.xml"), "w", encoding="utf-8") as fh:
                fh.write(t)
            print("   run accessions:", sorted(set(re.findall(r'<RUN accession="([^"]+)"', t)))[:5])


if __name__ == "__main__":
    main()
