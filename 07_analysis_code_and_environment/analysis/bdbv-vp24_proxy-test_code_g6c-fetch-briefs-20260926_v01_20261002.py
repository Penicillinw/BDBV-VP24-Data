"""G6c helper: fetch SOFT brief metadata for arbitrary GEO series and print a
design digest (sample titles + characteristics). Usage:

    python tools/g6c_fetch_briefs_20260926.py GSE314922 GSE314416 ...

Raw responses land in analysis/g6c_primary_immune_20260926/raw/<GSE>_brief.txt.
"""

from __future__ import annotations

import collections
import os
import re
import sys
import time
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926", "raw")


def get(url: str, tries: int = 4) -> str:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=180) as resp:  # noqa: S310
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(3 * (attempt + 1))
    print(f"FAILED {url}: {last}")
    return ""


def digest(acc: str, text: str) -> None:
    blocks = re.split(r"(?m)^\^SAMPLE = ", text)
    series = re.findall(r"(?m)^!Series_summary = (.*)$", text)
    design = re.findall(r"(?m)^!Series_overall_design = (.*)$", text)
    print(f"\n=== {acc} === samples: {len(blocks) - 1}")
    for line in (design or series)[:3]:
        print(f"   design: {line[:400]}")
    titles: list[str] = []
    char_counter: collections.Counter = collections.Counter()
    for block in blocks[1:]:
        ti = re.search(r"(?m)^!Sample_title = (.*)$", block)
        if ti:
            titles.append(ti.group(1))
        for value in re.findall(r"(?m)^!Sample_characteristics_ch1 = (.*)$", block):
            char_counter[value] += 1
    for value, count in char_counter.most_common(16):
        print(f"   [{count:4d}] {value[:140]}")
    for title in titles[:16]:
        print(f"      - {title[:140]}")
    if len(titles) > 16:
        print(f"      ... +{len(titles) - 16} more titles")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(RAW, exist_ok=True)
    for acc in sys.argv[1:]:
        path = os.path.join(RAW, f"{acc}_brief.txt")
        if os.path.exists(path) and os.path.getsize(path) > 2000:
            text = open(path, encoding="utf-8", errors="replace").read()
        else:
            text = get(
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                f"?acc={acc}&targ=all&form=text&view=brief"
            )
            if text:
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(text)
        if text:
            digest(acc, text)


if __name__ == "__main__":
    main()
