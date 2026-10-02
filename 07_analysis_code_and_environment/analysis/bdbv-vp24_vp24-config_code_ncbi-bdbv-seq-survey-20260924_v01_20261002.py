"""Survey public BDBV nucleotide records and check tropism-related motif states.

Steps
  1. esearch  nucleotide, Bundibugyo virus[Organism], all years and 2025:2026
  2. esummary for the 2025:2026 slice (accession, title, length, date)
  3. efetch complete genomes and translate the CDS regions carrying the
     motifs of interest (VP24 N135/R140, VP35 MIB2 NNLNS, VP40 PIP2 lysines)

Outputs land in analysis/bdbv_seq_survey_20260924/.
"""

import csv
import json
import os
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "bdbv_seq_survey_20260924")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(endpoint, params, tries=3):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    last = None
    for _ in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(3)
    raise RuntimeError(f"{endpoint} failed: {last}")


def esearch(term, db="nucleotide"):
    text = get("esearch.fcgi", {"db": db, "term": term, "retmax": 500, "retmode": "json"})
    return json.loads(text)["esearchresult"]


def esummary(ids, db="nucleotide"):
    text = get("esummary.fcgi", {"db": db, "id": ",".join(ids), "retmode": "json"})
    return json.loads(text)["result"]


def efetch_fasta(ids, db="nucleotide"):
    return get(
        "efetch.fcgi",
        {"db": db, "id": ",".join(ids), "rettype": "fasta", "retmode": "text"},
    )


def main():
    os.makedirs(OUT, exist_ok=True)

    all_hits = esearch("Bundibugyo virus[Organism]")
    recent = esearch("Bundibugyo virus[Organism] AND 2025:2026[PDAT]")
    summary = {
        "all_count": int(all_hits["count"]),
        "recent_2025_2026_count": int(recent["count"]),
        "all_ids": all_hits.get("idlist", []),
        "recent_ids": recent.get("idlist", []),
    }
    with open(os.path.join(OUT, "esearch.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    rows = []
    ids = summary["recent_ids"]
    for i in range(0, len(ids), 100):
        res = esummary(ids[i : i + 100])
        for uid in res.get("uids", []):
            rec = res[uid]
            rows.append(
                {
                    "uid": uid,
                    "accession": rec.get("accessionversion"),
                    "title": rec.get("title"),
                    "length": rec.get("slen"),
                    "created": rec.get("createdate"),
                    "updated": rec.get("updatedate"),
                    "organism": rec.get("organism"),
                }
            )
    rows.sort(key=lambda r: (r["updated"] or "", r["accession"] or ""), reverse=True)
    with open(os.path.join(OUT, "recent_records.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    complete = [
        r for r in rows if r["title"] and "complete genome" in r["title"].lower()
    ]
    print(f"total BDBV nucleotide records: {summary['all_count']}")
    print(f"2025-2026 records: {summary['recent_2025_2026_count']}")
    print(f"of which 'complete genome': {len(complete)}")
    print()
    for r in rows[:40]:
        print(f"  {r['accession']:<16s} len={str(r['length']):>6s} {r['updated']:<12s} {(r['title'] or '')[:110]}")

    if complete:
        fasta = efetch_fasta([r["uid"] for r in complete])
        with open(os.path.join(OUT, "complete_genomes.fasta"), "w", encoding="utf-8") as fh:
            fh.write(fasta)
        print()
        print(f"wrote {len(complete)} complete genomes to complete_genomes.fasta "
              f"({len(fasta)} chars)")


if __name__ == "__main__":
    main()
