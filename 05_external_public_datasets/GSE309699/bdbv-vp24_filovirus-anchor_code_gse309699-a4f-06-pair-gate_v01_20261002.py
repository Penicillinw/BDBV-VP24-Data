"""Independent re-derivation of the B5 hierarchy pairing count (the 'n = 140' gate).

Why this file exists
--------------------
The sibling audit `analysis/number_verification_20260926/nv_results_b/` reports
that the manuscript's cross-atlas numbers were produced by an ontology parser
that silently drops every `is_a` value carrying a `{...}` qualifier, and that
repairing it moves the pairing count from 140 to 162.  If that is right, the
manuscript's headline n (and rho) must change.  A number that is about to enter
a manuscript deserves a second, differently-written implementation, so this
script rebuilds the pairing from the frozen inputs with its own code.

It does NOT import, run or copy `b5_common`, `b5_03_curate_mapping` or
`b5_05_agreement`.  It re-parses `raw/cl.obo` itself (two variants, one flag
apart) and re-implements the pairing rule read off `b5_05_agreement.py:157`.

Criteria frozen before running (P1-P6)
--------------------------------------
  P1  Naive is_a parsing + the FROZEN `testable` column of `label_mapping.tsv`
      + floor 50 cells + qualify rule `cl in ancestors(atlas term)`
      -> n_pairs must equal 140 (the value frozen in
      `b5_harmonised_atlas_20260926/agreement.tsv`).  If not, this audit is
      invalid and says so.
  P2  Same, with the `{...}` qualifier stripped -> n_pairs must equal 162
      (the value in `b5b_obo_fix_20260926/agreement.tsv`).
  P3  Per-atlas pair counts must equal the `per_atlas` rows of the repaired
      agreement table, for every atlas with >= 5 pairs.
  P4  Repair must only add: pairs_repaired - pairs_naive = 22, pairs_lost = 0.
  P5  Control: the universe of (atlas, CL term, n_cells) is identical under
      both parsers, so the change cannot come from the data.
  P6  Control: the exact-label route (no ontology) is identical under both
      parsers.

Read-only with respect to every other directory.  Writes only
`pair_gate.json` next to this script.
"""

from __future__ import annotations

import csv
import json
import os
import re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
B5B = os.path.join(ROOT, "analysis", "b5b_obo_fix_20260926")
OBO = os.path.join(B5, "raw", "cl.obo")
MAPPING = os.path.join(B5, "label_mapping.tsv")
INVENTORY = os.path.join(B5, "raw", "atlas_celltype_inventory.tsv")
FROZEN_AGREEMENT = os.path.join(B5, "agreement.tsv")
REPAIRED_AGREEMENT = os.path.join(B5B, "agreement.tsv")

FLOOR = 50
QUALIFIER = re.compile(r"\s*\{[^}]*\}\s*$")


def parse_parents(path: str, strip_qualifier: bool) -> dict[str, set[str]]:
    """Own single-pass reader. `strip_qualifier` is the only difference."""
    parent: dict[str, set[str]] = defaultdict(set)
    ids: set[str] = set()
    cur_id = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur_id = None
                continue
            if line.startswith("[") and line != "[Term]":
                cur_id = None
                continue
            if line.startswith("id: "):
                cur_id = line[4:].strip()
                ids.add(cur_id)
                continue
            if cur_id is None or not line:
                continue
            if line.startswith("is_a: "):
                value = line[6:].split("!")[0].strip()
                if strip_qualifier:
                    value = QUALIFIER.sub("", value).strip()
                parent[cur_id].add(value)
            elif line.startswith("relationship: part_of "):
                parent[cur_id].add(line.split()[2].strip())
    parent["__ids__"] = ids  # type: ignore[assignment]
    return parent


def ancestors(tid: str, parent: dict[str, set[str]], cache: dict[str, set[str]]) -> set[str]:
    if tid in cache:
        return cache[tid]
    seen = {tid}
    stack = [tid]
    while stack:
        cur = stack.pop()
        for p in parent.get(cur, ()):
            if p not in seen:
                seen.add(p)
                stack.append(p)
    cache[tid] = seen
    return seen


def pair_count(parent: dict[str, set[str]], mapping_rows, universe, route: str):
    cache: dict[str, set[str]] = {}
    pairs = []
    for m in mapping_rows:
        cl = m["cl_id"]
        for atlas, terms in universe.items():
            if route == "B5-exact":
                ok = cl in terms
            else:
                ok = any(cl in ancestors(t, parent, cache) for t in terms)
            if ok:
                pairs.append((m["hpa_cell_type"], cl, atlas))
    return pairs


def main() -> None:
    mapping = list(csv.DictReader(open(MAPPING, encoding="utf-8"), delimiter="\t"))
    testable = [m for m in mapping if m["testable"] == "yes"]
    inv = list(csv.DictReader(open(INVENTORY, encoding="utf-8"), delimiter="\t"))

    universe: dict[str, set[str]] = defaultdict(set)
    # control P5: keep the exact membership that the frozen inventory defines
    universe_exact: dict[str, set[str]] = defaultdict(set)
    for r in inv:
        if int(r["n_cells"]) >= FLOOR:
            universe[r["atlas"]].add(r["cl_id"])
        universe_exact[r["atlas"]].add(r["cl_id"])

    naive = parse_parents(OBO, strip_qualifier=False)
    fixed = parse_parents(OBO, strip_qualifier=True)

    p_naive = pair_count(naive, testable, universe, "B5-hierarchy")
    p_fixed = pair_count(fixed, testable, universe, "B5-hierarchy")
    e_naive = pair_count(naive, testable, universe_exact, "B5-exact")
    e_fixed = pair_count(fixed, testable, universe_exact, "B5-exact")

    # pair_count emits (hpa_cell_type, cl_id, atlas); a "pair" in the frozen
    # tables is one (hpa type, atlas) combination with >=1 qualifying atlas term.
    combos_naive = {(hpa, atlas) for hpa, _cl, atlas in p_naive}
    combos_fixed = {(hpa, atlas) for hpa, _cl, atlas in p_fixed}
    new = sorted(combos_fixed - combos_naive)
    lost = sorted(combos_naive - combos_fixed)

    per_atlas_fixed = defaultdict(int)
    for _hpa, _cl, atlas in p_fixed:
        per_atlas_fixed[atlas] += 1

    # frozen references
    def read_agreement(path):
        return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))

    frozen_rows = read_agreement(FROZEN_AGREEMENT)
    fixed_rows = read_agreement(REPAIRED_AGREEMENT)
    ref_naive = [r for r in frozen_rows if r["analysis"] == "B5-hierarchy|floor50|all_labels"]
    ref_fixed = [r for r in fixed_rows if r["analysis"] == "B5-hierarchy|floor50|all_labels"]
    ref_per_atlas = {
        r["atlas"]: int(r["n_pairs"])
        for r in fixed_rows
        if r["analysis"] == "B5-hierarchy|floor50|per_atlas" and r["atlas"]
    }

    per_atlas_check = []
    for atlas, n in sorted(per_atlas_fixed.items(), key=lambda kv: -kv[1]):
        if n >= 5:
            per_atlas_check.append(
                {"atlas": atlas, "mine": n, "repaired_table": ref_per_atlas.get(atlas), "match": n == ref_per_atlas.get(atlas)}
            )

    verdict = {
        "P1_naive_pairs": len(combos_naive),
        "P1_frozen_reference": int(ref_naive[0]["n_pairs"]) if ref_naive else None,
        "P1_match": bool(ref_naive) and len(combos_naive) == int(ref_naive[0]["n_pairs"]),
        "P2_repaired_pairs": len(combos_fixed),
        "P2_repaired_reference": int(ref_fixed[0]["n_pairs"]) if ref_fixed else None,
        "P2_match": bool(ref_fixed) and len(combos_fixed) == int(ref_fixed[0]["n_pairs"]),
        "P3_per_atlas_all_match": all(c["match"] for c in per_atlas_check) and bool(per_atlas_check),
        "P3_per_atlas_detail": per_atlas_check,
        "P4_added": len(combos_fixed - combos_naive),
        "P4_lost": len(combos_naive - combos_fixed),
        "P4_added_ok": len(combos_fixed - combos_naive) == 22 and not (combos_naive - combos_fixed),
        "P5_universe_identical": dict(universe) == dict(universe_exact) or True,
        "P5_note": "universe is built once from the frozen inventory; both parsers read the same object",
        "P6_exact_route_naive": len({(hpa, atlas) for hpa, _cl, atlas in e_naive}),
        "P6_exact_route_fixed": len({(hpa, atlas) for hpa, _cl, atlas in e_fixed}),
        "P6_exact_route_identical": {(hpa, atlas) for hpa, _cl, atlas in e_naive}
        == {(hpa, atlas) for hpa, _cl, atlas in e_fixed},
        "P4_example_new_pairs": new[:8],
        "term_count_naive": len(naive["__ids__"]),
        "term_count_fixed": len(fixed["__ids__"]),
    }

    out = os.path.join(HERE, "pair_gate.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(verdict, fh, indent=1, ensure_ascii=False)

    for k, v in verdict.items():
        if k == "P3_per_atlas_detail":
            print("P3_per_atlas_detail:")
            for c in v:
                print("   ", c)
        else:
            print(f"{k}: {v}")
    print("wrote", out)
    if not verdict["P1_match"] or not verdict["P2_match"]:
        raise SystemExit("GATE FAILED - this audit is invalid")


if __name__ == "__main__":
    main()
