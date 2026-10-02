"""B5 step 2 - ontology-driven map from HPA cell-type names to Cell Ontology.

Route (frozen before any agreement statistic was computed):
  1. parse cl.obo (downloaded from purl.obolibrary.org, snapshot kept in raw/)
  2. for every HPA single-cell-type name, propose CL candidates by
       a. exact match on the CL term label
       b. normalised match (lower-case, plural -> singular, punctuation)
       c. match on a CL synonym
       d. token-overlap fallback (reported, never accepted automatically)
  3. restrict to CL terms that are ancestors-or-selves of a label that is
     actually present in a local atlas (there is no point mapping to a term
     that no independent dataset annotates)
  4. write the candidate table for explicit curation; the curated file
     (label_mapping_frozen.tsv) is hand-checked and is the frozen artefact.

Outputs: raw/hpa_cl_candidates.tsv, raw/cl_terms_used.tsv
"""

from __future__ import annotations

import csv
import json
import os
import re
from collections import defaultdict

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
RAW = os.path.join(OUT, "raw")
OBO = os.path.join(RAW, "cl.obo")
WIDE = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925", "hpa_ifn_landscape_wide.tsv")

# Curated shorthand used by HPA (frozen; these are HPA's own column values).
HPA_ALIASES = {
    "b-cells": "B cell",
    "t-cells": "T cell",
    "nk-cells": "natural killer cell",
    "pdcs": "plasmacytoid dendritic cell, human",
    "cdc": "conventional dendritic cell",
    "microglia": "microglial cell",
    "erythrocytes": "erythrocyte",
    "platelets": "platelet",
    "hepatocytes": "hepatocyte",
    "colonocytes": "colonocyte",
    "podocytes": "podocyte",
    "melanocytes": "melanocyte",
    "cardiomyocytes": "cardiac muscle cell",
    "somatotrophs": "somatotroph",
    "lactotrophs": "lactotroph",
    "corticotrophs": "corticotroph",
    "thyrotrophs": "thyrotroph",
    "gonadotrophs": "gonadotroph",
}


def parse_obo(path: str) -> dict[str, dict[str, object]]:
    terms: dict[str, dict[str, object]] = {}
    cur: dict[str, object] | None = None
    in_term = False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                if cur and cur.get("id"):
                    terms[str(cur["id"])] = cur
                cur = {"name": "", "is_a": [], "synonyms": [], "obsolete": False, "parents_alt": []}
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                if cur and cur.get("id"):
                    terms[str(cur["id"])] = cur
                cur = None
                in_term = False
                continue
            if not in_term or cur is None or not line:
                continue
            if line.startswith("id: "):
                cur["id"] = line[4:].strip()
            elif line.startswith("name: "):
                cur["name"] = line[6:].strip()
            elif line.startswith("is_a: "):
                cur["is_a"].append(line[6:].split("!")[0].strip())
            elif line.startswith("relationship: part_of "):
                cur["parents_alt"].append(line.split()[2].strip())
            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "(.*?)" (\w+)', line)
                if m:
                    cur["synonyms"].append(m.group(1))
            elif line.startswith("is_obsolete: true"):
                cur["obsolete"] = True
    if cur and cur.get("id"):
        terms[str(cur["id"])] = cur
    return terms


def norm(text: str) -> str:
    text = text.lower().strip()
    text = text.replace("alpha-beta", "alphabeta").replace("-", " ")
    text = re.sub(r"[^a-z0-9 ]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    words = []
    for w in text.split(" "):
        if w.endswith("s") and not w.endswith("ss") and len(w) > 3:
            w = w[:-1]
        words.append(w)
    return " ".join(sorted(words))


def main() -> None:
    terms = parse_obo(OBO)
    print(f"CL terms parsed: {len(terms)}")

    # CL terms actually used by the atlases
    used = list(csv.DictReader(open(os.path.join(RAW, "atlas_cl_union.tsv"), encoding="utf-8"), delimiter="\t"))
    used_ids = {r["cl_id"] for r in used}

    # ancestors (reflexive transitive closure) over is_a + part_of
    parent: dict[str, set[str]] = defaultdict(set)
    for tid, rec in terms.items():
        for p in list(rec["is_a"]) + list(rec["parents_alt"]):  # type: ignore[arg-type]
            parent[tid].add(p)

    def ancestors(tid: str) -> set[str]:
        seen = {tid}
        stack = [tid]
        while stack:
            cur = stack.pop()
            for p in parent.get(cur, ()):  # type: ignore[arg-type]
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        return seen

    cache: dict[str, set[str]] = {}

    def anc_cached(tid: str) -> set[str]:
        if tid not in cache:
            cache[tid] = ancestors(tid)
        return cache[tid]

    # index CL names/synonyms -> ids (excluding obsolete)
    by_norm: dict[str, set[str]] = defaultdict(set)
    by_norm_syn: dict[str, set[str]] = defaultdict(set)
    for tid, rec in terms.items():
        if rec["obsolete"]:
            continue
        if rec["name"]:
            by_norm[norm(str(rec["name"]))].add(tid)
        for s in rec["synonyms"]:  # type: ignore[union-attr]
            by_norm_syn[norm(str(s))].add(tid)

    hpa = [r["cell_type"] for r in csv.DictReader(open(WIDE, encoding="utf-8"), delimiter="\t")]
    rows = []
    for name in hpa:
        query = HPA_ALIASES.get(name, name)
        nq = norm(query)
        cands: dict[str, str] = {}
        for cid in by_norm.get(nq, set()):
            cands[cid] = "exact_label"
        for cid in by_norm_syn.get(nq, set()):
            cands.setdefault(cid, "synonym")
        # token-overlap fallback: CL names sharing >=2 tokens with the query
        qtok = set(nq.split())
        for tid, rec in terms.items():
            if rec["obsolete"] or not rec["name"]:
                continue
            ntok = set(norm(str(rec["name"])).split())
            if len(qtok & ntok) >= 2 and len(qtok) <= 4:
                cands.setdefault(tid, "token_overlap")
        cand_list = []
        for cid, how in cands.items():
            desc = anc_cached(cid) & used_ids
            cand_list.append(
                {
                    "hpa_cell_type": name,
                    "cl_id": cid,
                    "cl_label": terms.get(cid, {}).get("name", ""),
                    "match_type": how,
                    "atlas_descendant_terms": ";".join(sorted(desc)),
                    "n_atlas_terms_covered": len(desc),
                }
            )
        cand_list.sort(key=lambda r: (-int(r["n_atlas_terms_covered"]), r["match_type"]))
        rows.extend(cand_list)

    with open(os.path.join(RAW, "hpa_cl_candidates.tsv"), "w", encoding="utf-8", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    used_rows = [
        {
            "cl_id": r["cl_id"],
            "cl_label": r["cl_labels"],
            "total_cells": r["total_cells"],
            "n_atlases": r["n_atlases"],
            "atlas_ids_in_union": r["atlases"],
        }
        for r in used
    ]
    with open(os.path.join(RAW, "cl_terms_used.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(used_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(used_rows)

    n_with = len({r["hpa_cell_type"] for r in rows if int(r["n_atlas_terms_covered"]) > 0})
    print(f"HPA names: {len(hpa)}; with >=1 atlas-covering candidate: {n_with}")
    print("wrote candidates for curation")


if __name__ == "__main__":
    main()
