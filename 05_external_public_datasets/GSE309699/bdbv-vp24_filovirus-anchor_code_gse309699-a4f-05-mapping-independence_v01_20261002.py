"""Does the registered OBO `is_a` parsing defect reach `label_mapping.tsv`'s cl_id?

Criteria frozen before running (C1-C5):

  C1  Two parsers of `raw/cl.obo` are built here, sharing one loop; the ONLY
      difference is how an `is_a:` value is cleaned:
        naive   : value.split("!")[0].strip()            (as in b5_common.parse_obo)
        cleaned : value.split("!")[0], then strip a trailing `{...}` qualifier
      The term *metadata* tables must then be compared, not assumed equal:
          term id -> (name, tuple(sorted(synonyms)), is_obsolete)
      Verdict C1 = IDENTICAL if and only if the two tables are equal as mappings.

  C2  The set of term ids must be equal (no term may be dropped or added).

  C3  An `is_a` value that is not a term id in the file is a silent edge loss.
      Count them under the naive parser; under the cleaned parser the count of
      such values must be reported as well (it need not be zero - curated
      ontology files also carry non-CL parents).

  C4  Recompute the ancestor-closure damage independently: number of terms whose
      reflexive transitive closure changes, and the total number of
      (term, ancestor) pairs gained.  These are cross-checks of the numbers
      reported by the sibling audit (1320 edges / 3756 terms / 25241 pairs).

  C5  Decision rule stated in advance: if C1 is IDENTICAL and C2 holds, then
      `label_mapping.tsv`'s `cl_id`, `mapping_route` and `evidence` columns
      cannot be affected by the defect - they read only CURATED, the exact-label
      index and the synonym index, none of which calls onto.ancestors().  Only
      `testable` and the `covering_*` columns can move, and those were already
      re-derived by the sibling audit (79 -> 81 testable types, pairing count
      unchanged at 162).  Verdict is written to mapping_independence.json.

Read-only with respect to every other directory.  Writes only
`mapping_independence.json` next to this script.
"""

from __future__ import annotations

import json
import os
import re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OBO = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926", "raw", "cl.obo")

QUALIFIER = re.compile(r"\s*\{[^}]*\}\s*$")


def parse(path: str, clean_qualifier: bool):
    """Single loop; the flag is the only behavioural difference."""
    meta: dict[str, tuple] = {}
    isa: dict[str, list[str]] = defaultdict(list)
    alt: dict[str, list[str]] = defaultdict(list)
    cur: dict | None = None
    in_term = False

    def flush(rec):
        if not rec or not rec.get("id"):
            return
        tid = rec["id"]
        meta[tid] = (rec["name"], tuple(sorted(rec["syn"])), rec["obs"])
        isa[tid] = list(rec["isa"])
        alt[tid] = list(rec["alt"])

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                flush(cur)
                cur = {"id": None, "name": "", "isa": [], "syn": [], "obs": False, "alt": []}
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                flush(cur)
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
                value = line[6:].split("!")[0].strip()
                if clean_qualifier:
                    value = QUALIFIER.sub("", value).strip()
                cur["isa"].append(value)
            elif line.startswith("relationship: part_of "):
                cur["alt"].append(line.split()[2].strip())
            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "(.*?)" (\w+)', line)
                if m:
                    cur["syn"].append(m.group(1))
            elif line.startswith("is_obsolete: true"):
                cur["obs"] = True
    flush(cur)
    return meta, isa, alt


def closure_counts(ids, isa, alt):
    parent: dict[str, set] = defaultdict(set)
    for tid in ids:
        for p in isa.get(tid, []) + alt.get(tid, []):
            parent[tid].add(p)
    out = {}
    for tid in ids:
        seen = {tid}
        stack = [tid]
        while stack:
            cur = stack.pop()
            for p in parent.get(cur, ()):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        out[tid] = seen
    return out


def main() -> None:
    naive_meta, naive_isa, naive_alt = parse(OBO, clean_qualifier=False)
    clean_meta, clean_isa, clean_alt = parse(OBO, clean_qualifier=True)

    c1 = naive_meta == clean_meta
    c2 = set(naive_meta) == set(clean_meta)
    ids = sorted(clean_meta)

    bad_naive = [
        (tid, v)
        for tid in ids
        for v in naive_isa.get(tid, [])
        if v not in clean_meta
    ]
    bad_clean = [
        (tid, v)
        for tid in ids
        for v in clean_isa.get(tid, [])
        if v not in clean_meta
    ]
    # duplicated ids that only become valid once the qualifier is stripped
    repaired = [
        (tid, v, QUALIFIER.sub("", v).strip())
        for tid, v in bad_naive
        if QUALIFIER.sub("", v).strip() in clean_meta
    ]

    anc_naive = closure_counts(ids, naive_isa, naive_alt)
    anc_clean = closure_counts(ids, clean_isa, clean_alt)
    changed_terms = [t for t in ids if anc_naive[t] != anc_clean[t]]
    gained = sum(len(anc_clean[t] - anc_naive[t]) for t in ids)
    lost = sum(len(anc_naive[t] - anc_clean[t]) for t in ids)

    # C5 support: the exact-label and synonym indices, built the way b5_03 does
    def norm(s: str) -> str:
        s = s.lower().strip()
        s = re.sub(r"[^a-z0-9 ]+", " ", s)
        s = re.sub(r"\s+", " ", s).strip()
        if s.endswith("s") and not s.endswith("ss"):
            s = s[:-1]
        return s

    def indexes(meta):
        by_label, by_syn = defaultdict(set), defaultdict(set)
        for tid, (name, syns, obs) in meta.items():
            if obs:
                continue
            by_label[norm(name)].add(tid)
            for s in syns:
                by_syn[norm(s)].add(tid)
        return {k: frozenset(v) for k, v in by_label.items()}, {
            k: frozenset(v) for k, v in by_syn.items()
        }

    n_lab, n_syn = indexes(naive_meta)
    c_lab, c_syn = indexes(clean_meta)

    verdict = {
        "C1_term_metadata_identical": bool(c1),
        "C2_term_id_set_identical": bool(c2),
        "C1_label_index_identical": n_lab == c_lab,
        "C1_synonym_index_identical": n_syn == c_syn,
        "C3_isa_values_not_a_term_id_naive": len(bad_naive),
        "C3_isa_values_not_a_term_id_cleaned": len(bad_clean),
        "C3_values_repaired_by_stripping_qualifier": len(repaired),
        "C4_terms_with_changed_ancestor_closure": len(changed_terms),
        "C4_ancestor_pairs_gained": gained,
        "C4_ancestor_pairs_lost": lost,
        "C5_mapping_columns_affected": (
            "NO - cl_id / mapping_route / evidence read only CURATED, the exact-label "
            "index and the synonym index, none of which uses onto.ancestors(); both "
            "indices are identical between the two parsers."
            if (c1 and c2 and n_lab == c_lab and n_syn == c_syn)
            else "YES - investigate"
        ),
    }

    out = os.path.join(HERE, "mapping_independence.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(verdict, fh, indent=1, ensure_ascii=False)

    for k, v in verdict.items():
        print(f"{k}: {v}")
    print("wrote", out)


if __name__ == "__main__":
    main()
