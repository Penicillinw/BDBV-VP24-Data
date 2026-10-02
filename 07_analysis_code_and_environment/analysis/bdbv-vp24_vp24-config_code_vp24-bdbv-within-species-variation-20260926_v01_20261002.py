"""How much does BDBV VP24 actually vary within the species?

Question: the manuscript calls S83/Q135/H140/A141 "invariant". Invariance across a
small, near-identical set of sequences carries almost no statistical information
unless we know how many *other* positions vary. This script quantifies that: it
downloads every full-length BDBV VP24 protein record, counts the variable sites,
and reports what an invariance-at-four-sites observation is worth.

Read-only w.r.t. project data; writes only into
analysis/vp24_bdbv_census_20260925/coordinator_crosscheck/.
"""

from __future__ import annotations

import csv
import json
import os
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "vp24_bdbv_census_20260925", "coordinator_crosscheck")
RECORDS = os.path.join(OUT, "bdbv_vp24_protein_records.tsv")
FOUR = (83, 135, 140, 141)
LENGTH = 251


def fetch(accessions: list[str]) -> dict[str, str]:
    body = ("db=protein&rettype=fasta&retmode=text&id=" + ",".join(accessions)).encode()
    req = Request(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        data=body,
        headers={"User-Agent": "CodexResearch/1.0"},
    )
    with urlopen(req, timeout=300) as resp:  # noqa: S310 (fixed https host)
        text = resp.read().decode("utf-8", "replace")
    seqs: dict[str, str] = {}
    for record in text.split(">")[1:]:
        lines = record.splitlines()
        seqs[lines[0].split()[0]] = "".join(lines[1:]).strip().replace("*", "")
    return seqs


def main() -> None:
    with open(RECORDS, encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["class"] == "full_length"]
    seqs = fetch([r["protein_accession"] for r in rows])
    usable = {k: v for k, v in seqs.items() if len(v) == LENGTH}

    variable = [
        i + 1
        for i in range(LENGTH)
        if len({s[i] for s in usable.values()}) > 1
    ]
    residues = {
        str(p): sorted({s[p - 1] for s in usable.values()}) for p in variable
    }
    distinct = sorted(set(usable.values()))
    rate = len(variable) / LENGTH if usable else 0.0

    result = {
        "n_full_length_records": len(usable),
        "n_distinct_sequences": len(distinct),
        "vp24_length": LENGTH,
        "n_variable_positions": len(variable),
        "variable_positions": variable,
        "variable_position_residues": residues,
        "four_positions_are_variable": [p for p in FOUR if p in variable],
        "site_polymorphism_rate": round(rate, 5),
        "p_four_named_sites_invariant_if_drawn_at_random": round((1 - rate) ** 4, 4),
    }
    with open(os.path.join(OUT, "bdbv_within_species_variation.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)

    print(json.dumps(result, ensure_ascii=False, indent=2))
    print()
    print("Reading: with this many variable sites, invariance of any four named")
    print("positions is the expected outcome rather than evidence of constraint.")


if __name__ == "__main__":
    main()
