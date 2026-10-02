"""G6c step 2: fetch per-sample design metadata (SOFT brief) for the candidate
primary-immune-cell interferon datasets shortlisted from
analysis/g6c_primary_immune_20260926/datasets_screened.tsv.

Writes analysis/g6c_primary_immune_20260926/raw/<GSE>_brief.txt and prints a
design digest (sample title + characteristics) so the inclusion decision can be
made from the depositor's own description before any value is inspected.
"""

from __future__ import annotations

import os
import re
import sys
import time
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926")
RAW = os.path.join(OUT, "raw")

CANDIDATES = [
    "GSE306664",   # 5 donors x 4 immune subsets x IFN-a/b/g/l, scRNA (fixed 10x)
    "GSE147310",   # IRF1: human monocytes and macrophages +/- IFN
    "GSE147306",
    "GSE147309",
    "GSE147311",
    "GSE118305",   # monocyte-derived macrophages, large stimulus panel
    "GSE158434",   # IFN-stimulated and HIV-1-infected macrophages
    "GSE190623",   # FcgammaRIIa attenuates IFN-I response in monocytes
    "GSE193152",   # lipid accumulation, downregulated type-I IFN responses
    "GSE192709",
    "GSE125352",   # IFN-beta-inducible genes in human MDMs
    "GSE30536",    # IFN-alpha2-treated macrophages +/- HIV (time course)
    "GSE125817",   # moDC responses to HIV and other innate stimuli
    "GSE157857",   # CD1c+ DC and monocytes, IFN-I, CITE-seq
    "GSE201250",   # moDC +/- HIV, MDA5 knockdown
    "GSE242722",   # interferon / antiviral genes in macrophages
    "GSE244130",   # JAK inhibitor sensitivity in IFN-g-primed macrophages
    "GSE236156",   # alveolar macrophages vs monocyte-derived, TB/IFN
]


def get(url: str, tries: int = 4) -> str:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=180) as resp:  # noqa: S310 (fixed https host)
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - transient network failures
            last = exc
            time.sleep(3 * (attempt + 1))
    print(f"FAILED {url}: {last}")
    return ""


def brief(acc: str) -> str:
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={acc}&targ=all&form=text&view=brief"
    )
    return get(url)


def digest(acc: str, text: str) -> None:
    titles = re.findall(r"^!Sample_title = (.*)$", text, flags=re.M)
    chars = re.findall(r"^!Sample_characteristics_ch1 = (.*)$", text, flags=re.M)
    orgs = re.findall(r"^!Sample_organism_ch1 = (.*)$", text, flags=re.M)
    series_type = re.findall(r"^!Series_type = (.*)$", text, flags=re.M)
    platforms = sorted(set(re.findall(r"^!Series_platform_id = (.*)$", text, flags=re.M)))
    print(f"\n=== {acc} ===")
    print(f"samples: {len(titles)}  platforms: {','.join(platforms)}  type: {','.join(series_type)}")
    if orgs:
        from collections import Counter

        print("organisms:", Counter(orgs).most_common(3))
    if chars:
        from collections import Counter

        for value, count in Counter(chars).most_common(14):
            print(f"  [{count:3d}] {value[:150]}")
    for title in titles[:24]:
        print(f"   - {title[:150]}")
    if len(titles) > 24:
        print(f"   ... +{len(titles) - 24} more")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(RAW, exist_ok=True)
    for acc in CANDIDATES:
        path = os.path.join(RAW, f"{acc}_brief.txt")
        if os.path.exists(path) and os.path.getsize(path) > 2000:
            text = open(path, encoding="utf-8", errors="replace").read()
        else:
            text = brief(acc)
            if text:
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(text)
        if text:
            digest(acc, text)


if __name__ == "__main__":
    main()
