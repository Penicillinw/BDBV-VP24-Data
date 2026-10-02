"""Is GSE342661 usable in this environment? Check per-run sizes via ENA, and look for a
processed matrix or code deposit linked to the study."""

from __future__ import annotations

import csv
import io
import json
import os
import re
from urllib.request import Request, urlopen

RAW = r"G:\本迪布焦研究\analysis\g6g_gse342661_20260926\raw"
OUT = r"G:\本迪布焦研究\analysis\g6g_gse342661_20260926"


def get(url: str, tries: int = 3):
    last = None
    for a in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=300) as r:  # noqa: S310
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            import time

            time.sleep(2 * (a + 1))
    return f"ERROR: {last}"


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1508721"
           "&result=read_run&fields=run_accession,sample_title,fastq_bytes,"
           "library_strategy,instrument_platform,read_count&format=tsv")
    t = get(url)
    if t.startswith("ERROR"):
        print("ENA filereport failed:", t[:200])
    else:
        with open(os.path.join(RAW, "ena_runs.tsv"), "w", encoding="utf-8") as fh:
            fh.write(t)
        rows = list(csv.DictReader(io.StringIO(t), delimiter="\t"))
        print("runs:", len(rows))
        print("columns:", list(rows[0].keys()) if rows else [])
        total = 0
        sizes = []
        for r in rows:
            try:
                b = sum(int(x) for x in str(r.get("fastq_bytes", "")).split(";") if x.strip().isdigit())
            except Exception:  # noqa: BLE001
                b = 0
            sizes.append(b)
            total += b
        print(f"total fastq bytes: {total/1e9:.1f} GB; median run {sorted(sizes)[len(sizes)//2]/1e9:.2f} GB")
        titles = [r.get("sample_title", "") for r in rows]
        conditions = sorted({re.sub(r"biol rep \d+", "", t_).strip() for t_ in titles})
        print("condition labels:", conditions[:24], f"({len(conditions)} total)")

        # subset feasibility: control vs one IFN at 6h and 24h, both routes, all 5 reps
        for label in ["6h, Control, a", "6h, IFNB1, a", "6h, IFNL3, a"]:
            sub = [r for r in rows if label.split(",")[0] in r.get("sample_title", "")
                   and label.split(",")[1].strip() in r.get("sample_title", "")]
            print(f"  subset [{label}]: {len(sub)} runs")

    print("\n-- search for a processed deposit / code --")
    for q in ["GSE342661", "PRJNA1508721", "nasal epithelial interferon lambda beta air-liquid interface Chu"]:
        r = get("https://api.github.com/search/repositories?q=" + q.replace(" ", "+"))
        if not r.startswith("ERROR"):
            try:
                d = json.loads(r)
                print(f"  github '{q}': {d.get('total_count')} repos",
                      [i['full_name'] for i in d.get('items', [])[:3]])
            except Exception:  # noqa: BLE001
                print(f"  github '{q}': unparsable")
    r = get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
            "query=%22air-liquid%20interface%22%20AND%20%22interferon%20lambda%22%20AND%20nasal"
            "&format=json&pageSize=5")
    if not r.startswith("ERROR"):
        try:
            d = json.loads(r)
            for it in d.get("resultList", {}).get("result", [])[:5]:
                print("  EPMC:", it.get("pmid"), it.get("pubYear"), it.get("title", "")[:90])
        except Exception:  # noqa: BLE001
            pass


if __name__ == "__main__":
    main()
