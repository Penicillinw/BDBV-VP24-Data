"""G6 step 2: read the sample design of the shortlisted series.

For each candidate series, pull every GSM title and summarise how many distinct
cell types/cell lines were profiled with and without interferon.
"""

from __future__ import annotations

import json
import os
import re
import time
from collections import Counter, defaultdict
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "g6_ifn_responsiveness_20260926")
RAW = os.path.join(OUT, "raw")
CANDIDATES = ["GSE21158", "GSE186610", "GSE330780", "GSE306664", "GSE327707", "GSE54648"]


def get(url: str, tries: int = 4) -> str:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=120) as resp:  # noqa: S310
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"failed: {url}\n{last}")


def series_samples(acc: str) -> list[str]:
    """Parse the GEO series text view for its GSM list and titles."""
    text = get(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={acc}&targ=self&form=text&view=quick"
    )
    with open(os.path.join(RAW, f"{acc}_series.txt"), "w", encoding="utf-8") as fh:
        fh.write(text)
    return re.findall(r"=GSM\d+", text)


def sample_titles(gsms: list[str]) -> dict[str, str]:
    """Titles come from the same series text view: !Sample_title lines are ordered."""
    return {}


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    report: dict[str, object] = {}
    for acc in CANDIDATES:
        text = get(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
            f"?acc={acc}&targ=self&form=text&view=quick"
        )
        with open(os.path.join(RAW, f"{acc}_series.txt"), "w", encoding="utf-8") as fh:
            fh.write(text)
        gsms = re.findall(r"=GSM\d+", text)
        title_lines = re.findall(r"!Sample_title = (.+)", text)
        titles = {g.lstrip("="): t.strip() for g, t in zip(gsms, title_lines)}
        if not titles:
            # fall back: pair the ordered title lines with ordered accessions
            accs = re.findall(r"!Sample_geo_accession = (GSM\d+)", text)
            titles = {a: t.strip() for a, t in zip(accs, title_lines)}
        with open(os.path.join(RAW, f"{acc}_sample_titles.json"), "w", encoding="utf-8") as fh:
            json.dump(titles, fh, ensure_ascii=False, indent=1)
        treated = [t for t in titles.values() if re.search(r"ifn|interferon|ifna|ifnb|ifnl", t, re.I)]
        control = [t for t in titles.values() if re.search(r"control|untreated|mock|vehicle|nt\b|med", t, re.I)]
        report[acc] = {
            "n_samples": len(titles),
            "n_treated_like_titles": len(treated),
            "n_control_like_titles": len(control),
            "example_titles": list(titles.values())[:18],
        }
        print(f"\n===== {acc}  ({len(titles)} samples)")
        for t in list(titles.values())[:18]:
            print("   ", t[:120])
    with open(os.path.join(OUT, "candidate_designs.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print("\nwrote", os.path.join(OUT, "candidate_designs.json"))


if __name__ == "__main__":
    main()
