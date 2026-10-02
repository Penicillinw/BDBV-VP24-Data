"""Independent cross-check of the BDBV VP24 census claim (coordinator copy).

Enumerates every NCBI protein record matching `Bundibugyo virus[Organism] AND
VP24[All Fields]`, records its length, collection date, isolate and DBSOURCE
nucleotide accession, and writes a TSV so the count can be re-derived.

Read-only w.r.t. project data; writes only into
analysis/vp24_bdbv_census_20260925/coordinator_crosscheck/.
"""

from __future__ import annotations

import csv
import json
import os
import re
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "vp24_bdbv_census_20260925", "coordinator_crosscheck")
QUERY = "Bundibugyo virus[Organism] AND VP24[All Fields]"
EBOV_REF_ACC = "NC_002549.1"
POSITIONS = (83, 135, 140, 141, 184, 217)


def get(url: str) -> str:
    req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
    with urlopen(req, timeout=180) as resp:  # noqa: S310 (fixed https host)
        return resp.read().decode("utf-8", "replace")


def ebov_vp24() -> str:
    """VP24 protein of the EBOV type reference, via the CDS-protein FASTA route."""
    text = get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={EBOV_REF_ACC}&rettype=fasta_cds_aa&retmode=text"
    )
    for record in text.split(">")[1:]:
        lines = record.splitlines()
        if "[gene=VP24]" in lines[0]:
            return "".join(lines[1:]).strip().replace("*", "")
    raise SystemExit("EBOV VP24 CDS not found")


def fasta_proteins(accessions: list[str]) -> dict[str, str]:
    text = get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=protein&id={','.join(accessions)}&rettype=fasta&retmode=text"
    )
    out: dict[str, str] = {}
    for record in text.split(">")[1:]:
        lines = record.splitlines()
        acc = lines[0].split()[0]
        seq = "".join(lines[1:]).strip().replace("*", "")
        out[acc] = seq
        out[acc.split(".")[0]] = seq  # accept unversioned lookup keys
    return out


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    search = json.loads(
        get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=protein&retmax=500&retmode=json&term={quote(QUERY)}"
        )
    )
    ids = search["esearchresult"]["idlist"]
    total = search["esearchresult"]["count"]

    text = get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=protein&id={','.join(ids)}&rettype=gp&retmode=text"
    )

    rows = []
    for rec in text.split("\n//"):
        locus = re.search(r"^LOCUS\s+(\S+)\s+(\d+) aa", rec, re.M)
        if not locus:
            continue
        date = re.findall(r'/collection_date="([^"]+)"', rec)
        isolate = re.findall(r'/isolate="([^"]+)"', rec)
        country = re.findall(r'/country="([^"]+)"', rec)
        dbsource = re.search(r"DBSOURCE\s+accession\s+(\S+)", rec)
        definition = re.search(r"^DEFINITION\s+(.+?)\.?$", rec, re.M)
        rows.append(
            {
                "protein_accession": locus.group(1),
                "length_aa": int(locus.group(2)),
                "class": "full_length" if int(locus.group(2)) == 251 else "partial",
                "collection_date": date[0] if date else "",
                "country": country[0] if country else "",
                "isolate": isolate[0] if isolate else "",
                "nucleotide_accession": dbsource.group(1) if dbsource else "",
                "definition": (definition.group(1).strip() if definition else ""),
            }
        )

    rows.sort(key=lambda r: (r["length_aa"], r["collection_date"], r["protein_accession"]))

    ref = ebov_vp24()
    seqs = fasta_proteins([r["protein_accession"] for r in rows if r["class"] == "full_length"])
    for r in rows:
        seq = seqs.get(r["protein_accession"], "")
        r["same_length_as_ebov"] = len(seq) == len(ref)
        n_sub = (
            sum(1 for a, b in zip(seq, ref) if a != b) if r["class"] == "full_length" else ""
        )
        r["n_sub_vs_EBOV_unaligned"] = n_sub
        for pos in POSITIONS:
            r[f"pos{pos}"] = seq[pos - 1] if seq and len(seq) >= pos else ""
        r["four_position_motif"] = (
            "".join(r[f"pos{p}"] for p in (83, 135, 140, 141))
            if r["class"] == "full_length"
            else "NA"
        )

    tsv = os.path.join(OUT, "bdbv_vp24_protein_records.tsv")
    with open(tsv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    full = [r for r in rows if r["class"] == "full_length"]
    partial = [r for r in rows if r["class"] == "partial"]
    years = sorted({r["collection_date"][-4:] for r in full if r["collection_date"]})
    motifs = {}
    for r in full:
        motifs[r["four_position_motif"]] = motifs.get(r["four_position_motif"], 0) + 1
    summary = {
        "query": QUERY,
        "esearch_count": total,
        "records_parsed": len(rows),
        "full_length_records": len(full),
        "partial_records": len(partial),
        "partial_accessions": [r["protein_accession"] for r in partial],
        "full_length_collection_years": years,
        "full_length_collection_dates": sorted({r["collection_date"] for r in full}),
        "distinct_isolates": sorted({r["isolate"] for r in full}),
        "ebov_reference": EBOV_REF_ACC,
        "ebov_vp24_length": len(ref),
        "ebov_positions_reported": list(POSITIONS),
        "four_position_motif_counts": motifs,
        "motif_definition": "positions 83/135/140/141 in EBOV VP24 numbering",
        "all_full_length_are_SQHA": set(motifs) == {"SQHA"},
        "n_sub_vs_EBOV_unaligned_range": sorted(
            {r["n_sub_vs_EBOV_unaligned"] for r in full}
        ),
    }
    with open(os.path.join(OUT, "census_crosscheck_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print()
    print(
        f"{'accession':14s} {'aa':>4s} {'date':12s} {'sub':>4s} {'83':>3s} {'135':>4s} "
        f"{'140':>4s} {'141':>4s} {'184':>4s} {'217':>4s} {'nt':12s} isolate"
    )
    for r in rows:
        print(
            f"{r['protein_accession']:14s} {r['length_aa']:>4d} {r['collection_date']:12s} "
            f"{str(r['n_sub_vs_EBOV_unaligned']):>4s} {r['pos83']:>3s} {r['pos135']:>4s} "
            f"{r['pos140']:>4s} {r['pos141']:>4s} {r['pos184']:>4s} {r['pos217']:>4s} "
            f"{r['nucleotide_accession']:12s} {r['isolate'][:32]}"
        )


if __name__ == "__main__":
    main()
