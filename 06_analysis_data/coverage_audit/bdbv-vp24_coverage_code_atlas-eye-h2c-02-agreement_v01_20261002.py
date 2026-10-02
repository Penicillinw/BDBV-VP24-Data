"""H2c step 2 - re-run the low-end cross-atlas test with the retina connection
repaired, as explicitly labelled sensitivity arms.

Everything is held at the frozen rule set of ``h2_lowend_agreement_hierarchy.py``
and ``b5_05_agreement.py``:

  * atlas pseudobulk from raw integer counts -> CPM -> the three Delta
    components as the mean of log10(CPM+1).  The pseudobulk file is produced by
    h2's own ``h2_pseudobulk_new_atlases.py``; the count layer actually used is
    carried in its ``count_source`` column and is reported here per atlas
    instead of being assumed;
  * connection through Cell Ontology ids only; an atlas term matches a reference
    type when it is that term or a MORE SPECIFIC term below it (reflexive
    is_a + part_of closure);
  * cell-count floor 50 cells per (atlas, CL term);
  * z-scoring within atlas;
  * 4,000 label permutations, seed 20260926;
  * reference Delta from analysis/t4_ifn_landscape_20260925/restriction_gradient.tsv
    (frozen; never recomputed here).

Pre-registered arms (declared before any statistic in this file was computed):

  R0  frozen B5 mapping, unchanged.  GATE: must reproduce n = 38, rho = 0.568,
      permutation p = 0.0005 of h2b_agreement.json.  If the gate fails, the
      harness is wrong and no arm below is interpreted until it is fixed.
  R1  R0 + photoreceptor cells only (rod -> CL:0000604, cone -> CL:0000573).
      This is the arm that speaks to the unclosed "photoreceptor" item.
  R2  R1 + retinal interneurons (amacrine, bipolar, horizontal).
  R3  R2 + Mueller glia  = the whole retina family.

The added mappings are curated additions written by this task (see PROPOSED in
``h2c_01_label_diagnosis.py``).  They are NOT part of the frozen B5 artefact and
are therefore never mixed into R0.

HARNESS FIDELITY.  The gate is reproduced with the frozen script's own ancestor
loader, copied verbatim, because a first attempt that reused
``b5_common.Ontology`` did NOT reproduce it.  The discrepancy is itself a
finding and is reported in FINDINGS.md:

  (a) ``b5_common.parse_obo`` keeps the ``{is_inferred="true"}`` qualifier on
      qualified ``is_a`` lines, so the parent id becomes the corrupt string
      ``CL:0000630 {is_inferred="true"}``.  Every edge carrying such a qualifier
      is therefore dropped from the CL graph, which cost B5 and my first pass the
      mesothelial-cell < Sertoli-cell edge.  The frozen h2 loaders use
      ``line.split()[1]`` and are not affected.
  (b) the frozen h2 script drops the pairs of any atlas with fewer than three
      matched reference types from the pooled statistic (bare ``continue``),
      which is why the frozen pooled n is 38 and not 40.  Both conventions are
      reported here so the reader can see how much of any change is due to the
      repair and how much to that filter.

Also reported, per the standing evidence rules: the whole arms x floor grid,
achieved power / minimum detectable |rho|, leave-one-atlas-out for the repaired
arm, and a collateral-match audit proving each added term connects exactly the
intended reference type.

Writes only inside analysis/h2c_photoreceptor_20260926/.
Outputs: h2c_agreement.json, h2c_R0_pairs.tsv, h2c_pairs_R3.tsv,
         h2c_newly_connected.tsv
"""

from __future__ import annotations

import csv
import json
import math
import os
import statistics as st
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
H2 = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926")
REF = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                   "restriction_gradient.tsv")
PB = os.path.join(H2, "raw", "new_atlas_pseudobulk.tsv")
LAB = os.path.join(H2, "raw", "new_atlas_labels.tsv")
MAP = os.path.join(B5, "label_mapping.tsv")

sys.path.insert(0, B5)
from b5_common import min_detectable_rho  # noqa: E402

DRAWS, SEED = 4000, 20260926
FLOORS = [50, 200, 500]          # 50 is the frozen main floor

# Curated repair (this task); identical content to PROPOSED in h2c_01.
REPAIR = {
    "rod photoreceptor cells": "CL:0000604",
    "cone photoreceptor cells": "CL:0000573",
    "retinal amacrine cells": "CL:0000561",
    "retinal bipolar cells": "CL:0000748",
    "retinal horizontal cells": "CL:0000745",
    "m\u00fcller glia": "CL:0000636",
}
ARM_ADDITIONS = {
    "R0_frozen": [],
    "R1_photoreceptor": ["rod photoreceptor cells", "cone photoreceptor cells"],
    "R2_plus_interneurons": ["rod photoreceptor cells", "cone photoreceptor cells",
                             "retinal amacrine cells", "retinal bipolar cells",
                             "retinal horizontal cells"],
    "R3_whole_retina_family": list(REPAIR),
}


def load_tsv(path: str) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- ancestor loader copied verbatim from h2_lowend_agreement_hierarchy.py ----
def load_ancestors(obo_path: str) -> dict[str, set[str]]:
    parents: dict[str, set[str]] = {}
    cur = None
    with open(obo_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = None
            elif line.startswith("id: CL:"):
                cur = line.split("id: ")[1].strip()
                parents.setdefault(cur, set())
            elif cur and line.startswith("is_a: "):
                parents[cur].add(line.split()[1])
            elif cur and line.startswith("relationship: part_of "):
                parents[cur].add(line.split()[2])
    return parents


def closure(term: str, parents: dict[str, set[str]],
            memo: dict[str, set[str]]) -> set[str]:
    if term in memo:
        return memo[term]
    out = {term}
    for p in parents.get(term, ()):
        out |= closure(p, parents, memo)
    memo[term] = out
    return out


# --- frozen statistics, copied from the same script --------------------------
def spearman(x: list[float], y: list[float]) -> float:
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1
        return r
    rx, ry = rank(x), rank(y)
    mx, my = st.fmean(rx), st.fmean(ry)
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(len(x)))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float("nan")


def perm_p(x: list[float], y: list[float], draws: int = DRAWS,
           seed: int = SEED) -> float:
    obs, rng, yy = abs(spearman(x, y)), seed, list(y)
    n, ge = len(yy), 0
    for _ in range(draws):
        rng = (1103515245 * rng + 12345) % (2 ** 31)
        j = rng % n
        rng = (1103515245 * rng + 12345) % (2 ** 31)
        k = rng % n
        yy[j], yy[k] = yy[k], yy[j]
        if abs(spearman(x, yy)) >= obs:
            ge += 1
    return (ge + 1) / (draws + 1)


def compute_matches(ref, t2h, by_atlas, parents, memo, added, floor,
                    drop_atlas=None):
    """atlas -> {reference type: [(z, n_cells, atlas label), ...]} under the rule."""
    t2h = {t: list(v) for t, v in t2h.items()}
    for hpa in added:
        t2h.setdefault(REPAIR[hpa], []).append(hpa)
    result: dict[str, dict[str, list[tuple[float, int, str]]]] = {}
    for atlas, recs in sorted(by_atlas.items()):
        if atlas == drop_atlas:
            continue
        vals = [r["delta"] for r in recs]
        mu, sd = st.fmean(vals), (st.pstdev(vals) or 1.0)
        matches: dict[str, list[tuple[float, int, str]]] = defaultdict(list)
        for r in recs:
            if r["n_cells"] < floor:
                continue
            anc = closure(r["cl_id"], parents, memo) if r["cl_id"].startswith("CL:") \
                else {r["cl_id"]}
            for t in anc:
                for hpa in t2h.get(t, []):
                    if hpa in ref:
                        matches[hpa].append(((r["delta"] - mu) / sd,
                                             r["n_cells"], r["cl_label"]))
        result[atlas] = matches
    return result


def build_pairs(ref, t2h, by_atlas, parents, memo, added, floor,
                min_matches_filter: bool, drop_atlas=None):
    """Return (pairs, per_atlas, excluded_pairs) for the frozen connection rule.

    ``min_matches_filter`` True reproduces the frozen h2 script exactly: an atlas
    whose match count is below three contributes its per-atlas row but none of
    its pairs to the pooled statistic.
    """
    allm = compute_matches(ref, t2h, by_atlas, parents, memo, added, floor,
                           drop_atlas)
    pairs, per_atlas, excluded = [], [], []
    for atlas, matches in allm.items():
        per_atlas.append({"atlas": atlas, "n_matched": len(matches),
                          "rho": None})
        rows = []
        for h in sorted(matches):
            v = matches[h]
            rows.append({"atlas": atlas, "reference_cell_type": h,
                         "reference_delta": round(ref[h], 3),
                         "atlas_delta_z": round(st.fmean(x[0] for x in v), 3),
                         "atlas_labels": " | ".join(sorted({x[2] for x in v})),
                         "n_cells": sum(x[1] for x in v)})
        if len(matches) >= 3:
            rho = spearman([ref[h] for h in sorted(matches)],
                           [st.fmean(v[0] for v in matches[h])
                            for h in sorted(matches)])
            per_atlas[-1]["rho"] = None if math.isnan(rho) else round(rho, 3)
            pairs.extend(rows)
        elif min_matches_filter:
            excluded.extend(rows)
        else:
            pairs.extend(rows)
    return pairs, per_atlas, excluded


def pooled(pairs) -> dict:
    if len(pairs) < 3:
        return {"n": len(pairs), "rho": None, "p_permutation": None}
    xs = [p["reference_delta"] for p in pairs]
    ys = [p["atlas_delta_z"] for p in pairs]
    return {"n": len(pairs), "rho": round(spearman(xs, ys), 3),
            "p_permutation": round(perm_p(xs, ys), 4),
            "min_detectable_abs_rho_power80": round(min_detectable_rho(len(pairs)), 3)}


def write_tsv(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    parents = load_ancestors(os.path.join(B5, "raw", "cl.obo"))
    memo: dict[str, set[str]] = {}
    ref = {r["cell_type"]: float(r["delta_restriction"]) for r in load_tsv(REF)}
    t2h: dict[str, list[str]] = defaultdict(list)
    for r in load_tsv(MAP):
        t2h[r["cl_id"]].append(r["hpa_cell_type"])
    by_atlas: dict[str, list[dict]] = defaultdict(list)
    lab = {(r["atlas"], r["cl_id"]): r["cl_label"] for r in load_tsv(LAB)}
    for r in load_tsv(PB):
        by_atlas[r["atlas"]].append({"cl_id": r["cl_id"],
                                     "cl_label": lab.get((r["atlas"], r["cl_id"]),
                                                         r["cl_id"]),
                                     "n_cells": int(r["n_cells"]),
                                     "delta": float(r["delta"]),
                                     "count_source": r["count_source"]})

    out: dict[str, object] = {
        "rule": ("frozen h2/b5 rule set: raw counts -> CPM -> equal-weight "
                 "log10(CPM+1) components; ontology-id connection with "
                 "is_a + part_of ancestry; floor 50 cells; within-atlas z; "
                 "4000 label permutations, seed 20260926"),
        "arms_frozen_convention": {}, "arms_no_min_matches_filter": {},
        "floor_sensitivity": {}, "leave_one_atlas_out_R3": {},
        "collateral_audit": {},
        "count_layer_by_atlas": {a: sorted({r["count_source"] for r in recs})
                                 for a, recs in by_atlas.items()},
    }

    for conv, flag in (("arms_frozen_convention", True),
                       ("arms_no_min_matches_filter", False)):
        for arm, added in ARM_ADDITIONS.items():
            pairs, per_atlas, excl = build_pairs(ref, t2h, by_atlas, parents,
                                                 memo, added, 50, flag)
            out[conv][arm] = {"additions": {h: REPAIR[h] for h in added},
                              "per_atlas": per_atlas, "pooled": pooled(pairs),
                              "n_pairs_excluded_by_min_matches_3": len(excl),
                              "excluded_pairs": excl}

    # headline artefacts, on the frozen convention
    fz_arms = out["arms_frozen_convention"]  # type: ignore[index]
    write_tsv(os.path.join(HERE, "h2c_R0_pairs.tsv"),
              build_pairs(ref, t2h, by_atlas, parents, memo, [], 50, True)[0])
    write_tsv(os.path.join(HERE, "h2c_pairs_R3.tsv"),
              build_pairs(ref, t2h, by_atlas, parents, memo, list(REPAIR), 50,
                          True)[0])

    # gate diff against the frozen h2b output
    frozen_path = os.path.join(H2, "h2b_pairs.tsv")
    if os.path.exists(frozen_path):
        fz = {(r["atlas"], r["reference_cell_type"]) for r in load_tsv(frozen_path)}
        mine = {(p["atlas"], p["reference_cell_type"]) for p in
                build_pairs(ref, t2h, by_atlas, parents, memo, [], 50, True)[0]}
        out["gate_diff_vs_h2b"] = {
            "frozen_n": len(fz), "recomputed_n": len(mine),
            "only_in_recomputed": sorted(f"{a}::{h}" for a, h in mine - fz),
            "only_in_frozen": sorted(f"{a}::{h}" for a, h in fz - mine),
        }

    # floor sensitivity on the frozen convention, R0 vs R3
    for arm in ("R0_frozen", "R3_whole_retina_family"):
        out["floor_sensitivity"][arm] = {}
        for f in FLOORS:
            p, _, _ = build_pairs(ref, t2h, by_atlas, parents, memo,
                                  ARM_ADDITIONS[arm], f, True)
            out["floor_sensitivity"][arm][f"floor_{f}"] = pooled(p)

    # leave-one-atlas-out on the repaired arm (frozen convention)
    for atlas in sorted(by_atlas):
        p, _, _ = build_pairs(ref, t2h, by_atlas, parents, memo, list(REPAIR),
                              50, True, drop_atlas=atlas)
        out["leave_one_atlas_out_R3"][atlas.replace(".h5ad", "")] = pooled(p)

    # newly connected pairs (in R3, absent from R0), frozen convention
    base = {(p["atlas"], p["reference_cell_type"]) for p in
            build_pairs(ref, t2h, by_atlas, parents, memo, [], 50, True)[0]}
    fam = build_pairs(ref, t2h, by_atlas, parents, memo, list(REPAIR), 50,
                      True)[0]
    new = []
    for p in fam:
        if (p["atlas"], p["reference_cell_type"]) in base:
            continue
        q = dict(p)
        q["same_sign"] = (p["reference_delta"] > 0) == (p["atlas_delta_z"] > 0)
        new.append(q)
    write_tsv(os.path.join(HERE, "h2c_newly_connected.tsv"), new)
    out["newly_connected_pairs"] = new

    # pairs connected by the repair but dropped by the frozen min-matches rule
    p_all, _, _ = build_pairs(ref, t2h, by_atlas, parents, memo, list(REPAIR),
                              50, False)
    out["connected_but_excluded_by_min_matches_rule"] = [
        {"atlas": p["atlas"], "reference_cell_type": p["reference_cell_type"],
         "atlas_delta_z": p["atlas_delta_z"], "n_cells": p["n_cells"]}
        for p in p_all if (p["atlas"], p["reference_cell_type"]) not in
        {(q["atlas"], q["reference_cell_type"]) for q in fam}]

    # collateral audit: atlas terms each repaired term reaches, and whether any
    # other reference type shares the same CL id
    for hpa, cl in REPAIR.items():
        reached = []
        for atlas, recs in by_atlas.items():
            for r in recs:
                if cl in (closure(r["cl_id"], parents, memo)
                          if r["cl_id"].startswith("CL:") else {r["cl_id"]}) \
                        and r["n_cells"] >= 50:
                    reached.append(f"{atlas.replace('.h5ad', '')}:{r['cl_id']}"
                                   f"({r['n_cells']})")
        out["collateral_audit"][hpa] = {
            "cl_id": cl, "atlas_terms_reached": reached,
            "other_reference_types_sharing_this_cl":
                [h for h, c in REPAIR.items() if h != hpa and c == cl]}

    out["arms_frozen_convention"]["R0_frozen"]["gate"] = (
        "target from h2b_agreement.json: n=38, rho=0.568, p=0.0005")

    # --- degenerate atlases: within-atlas z is not on an absolute scale ------
    # An atlas contributing one CL term has z == 0 for every term; an atlas
    # contributing exactly two has z == +1 / -1.  For those, only the ORDER of
    # the matched reference types is interpretable, never the signed distance.
    deg = {}
    for arm, added in (("R0_frozen", []), ("R3_whole_retina_family", list(REPAIR))):
        for atlas, matches in compute_matches(ref, t2h, by_atlas, parents,
                                              memo, added, 50).items():
            if len(matches) > 2:
                continue
            names = sorted(matches)
            zs = {h: round(st.fmean(x[0] for x in matches[h]), 3) for h in names}
            all_equal = len(set(zs.values())) == 1
            ordinal = None
            if len(names) == 2 and not all_equal:
                ordinal = (sorted(names, key=lambda h: zs[h]) ==
                           sorted(names, key=lambda h: ref[h]))
            deg[f"{arm}::{atlas.replace('.h5ad', '')}"] = {
                "n_matched": len(names),
                "reference_types": names,
                "atlas_z": zs,
                "reference_delta": {h: round(ref[h], 3) for h in names},
                "pair_order_concordant_with_reference": ordinal,
                "reading": ("single-term atlas: z is identically 0 - no "
                            "variance-based information"
                            if len(names) == 1 else
                            "NO variance-based information: every matched type "
                            "sits on the same z (the atlas carries only one CL "
                            "term), so ties make the order meaningless too"
                            if all_equal else
                            "two-term atlas: z is identically +/-1 - only the "
                            "ORDER of the two types is interpretable, not the "
                            "signed distance"),
            }
    out["degenerate_atlas_readings"] = deg

    # Written under two names on purpose.  `h2c_agreement.json` is the
    # conventional name, but an unowned, incomplete script
    # (`h2c_photoreceptor_match.py`, present in this directory and not written by
    # this task) also declares that filename as its output, so the authoritative
    # copy is additionally kept under `h2c_agreement_retina.json` and never
    # silently overwritten by a foreign script.
    payload = json.dumps(out, ensure_ascii=False, indent=1)
    for name in ("h2c_agreement_retina.json", "h2c_agreement.json"):
        with open(os.path.join(HERE, name), "w", encoding="utf-8") as fh:
            fh.write(payload)

    print("gate diff vs h2b_pairs.tsv:",
          json.dumps(out.get("gate_diff_vs_h2b", {}), ensure_ascii=False))
    for conv in ("arms_frozen_convention", "arms_no_min_matches_filter"):
        print(f"\n[{conv}]")
        for arm in ARM_ADDITIONS:
            print(f"  {arm:26s} {json.dumps(out[conv][arm]['pooled'])}")
    print("\nnewly connected (frozen convention)")
    for p in new:
        print(f"  {p['atlas'][:26]:28s} {p['reference_cell_type']:26s} "
              f"ref {p['reference_delta']:+7.3f} atlas {p['atlas_delta_z']:+7.3f} "
              f"n={p['n_cells']:>6d} same_sign={p['same_sign']}")
    print("\nconnected but excluded by the frozen min-matches rule")
    for p in out["connected_but_excluded_by_min_matches_rule"]:
        print(f"  {p['atlas'][:26]:28s} {p['reference_cell_type']:26s} "
              f"atlas {p['atlas_delta_z']:+7.3f} n={p['n_cells']:>6d}")
    print("\ndegenerate atlases (only order is interpretable)")
    for k, v in deg.items():
        print(f"  {k}  n={v['n_matched']}  {v['atlas_z']}  "
              f"order_ok={v['pair_order_concordant_with_reference']}")


if __name__ == "__main__":
    main()
