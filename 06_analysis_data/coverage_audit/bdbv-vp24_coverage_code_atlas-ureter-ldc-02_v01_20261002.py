"""Why the two frozen anchors did not reproduce exactly.

Read-only diagnostics for the two failures printed by ldc_01:

  F1  19 atlases, convention A: the y values match b5b field for field
      (max |diff| = 0.0 over all 162 pairs) but the pooled rho is 0.4901
      against b5b's 0.5191.  A pooled Spearman is invariant under any strictly
      monotone transform of x, so the only way for y to be identical and rho to
      differ is that the two x vectors are not monotone-equivalent: b5b used
      Delta rebuilt from the HPA wide table, this run used the tabulated
      `delta_restriction`.  Measured here.

  F2   7 atlases, convention B: 41 pairs here against h2b's 38, and the
      per-pair z values differ by up to 0.598.  h2b z-scored each new atlas
      over every row of its pseudobulk table; this run restricts the atlas's
      universe to the labels the atlas inventory lists, at or above the floor.
      Measured here.

Nothing is written by this script.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
B5B = os.path.join(ROOT, "analysis", "b5b_obo_fix_20260926")
H2 = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926")
sys.path.insert(0, B5B)

from b5_common import load_hpa, hpa_delta, spearman, zscore  # noqa: E402

REF = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                   "restriction_gradient.tsv")
P07 = os.path.join(H2, "raw", "new_atlas_pseudobulk.tsv")
LBL07 = os.path.join(H2, "raw", "new_atlas_labels.tsv")
B5B_PAIRS = os.path.join(B5B, "b5_matched_pairs.tsv")
H2B_PAIRS = os.path.join(H2, "h2b_pairs.tsv")


def edge(name, text):
    print(f"\n=== {name} ===\n{text}")


def main():
    ref = {r["cell_type"]: float(r["delta_restriction"])
           for r in csv.DictReader(open(REF, encoding="utf-8"), delimiter="\t")}
    wide = load_hpa()
    rebuilt = {ct: hpa_delta(row) for ct, row in wide.items()}
    shared = sorted(set(ref) & set(rebuilt))
    rho = spearman([ref[c] for c in shared], [rebuilt[c] for c in shared])
    edge("F1a  reference Delta: tabulated vs rebuilt from the wide table",
         f"shared types           : {len(shared)}\n"
         f"spearman               : {rho:.6f}\n"
         f"ties in tabulated      : "
         f"{len(shared) - len({ref[c] for c in shared})}\n"
         f"ties in rebuilt        : "
         f"{len(shared) - len({round(rebuilt[c], 6) for c in shared})}")
    pairs = list(csv.DictReader(open(B5B_PAIRS, encoding="utf-8"),
                                delimiter="\t"))
    xs_tab = [ref[p["hpa_cell_type"]] for p in pairs]
    xs_reb = [zscore(rebuilt)[p["hpa_cell_type"]] for p in pairs]
    ys = [float(p["y_all_labels"]) for p in pairs]
    edge("F1b  the same 162 pairs under the two x conventions",
         f"x = tabulated delta_restriction : rho = {spearman(xs_tab, ys):.4f}\n"
         f"x = rebuilt wide (b5b)          : rho = {spearman(xs_reb, ys):.4f}\n"
         f"b5b agreement.tsv               : rho = 0.5191")
    diff = [(a, b) for a, b in zip(sorted(ref, key=ref.get),
                                   sorted(rebuilt, key=rebuilt.get)) if a != b]
    edge("F1c  first order disagreements between the two reference vectors",
         "\n".join(f"  {a:38s} vs {b}" for a, b in diff[:12])
         + f"\n  ... {len(diff)} / {len(shared)} positions differ" if diff
         else "  none")

    rows = defaultdict(list)
    for r in csv.DictReader(open(P07, encoding="utf-8"), delimiter="\t"):
        rows[r["atlas"]].append(r)
    uni = defaultdict(list)
    for r in csv.DictReader(open(LBL07, encoding="utf-8"), delimiter="\t"):
        if (r.get("cl_id") or "").strip():
            uni[r["atlas"]].append(r)
    lines = [f"  {'atlas':44s} {'pb':>4} {'pb<50':>6} {'lbl':>4} "
             f"{'lbl>=50':>8} {'same set':>9}"]
    for atlas in sorted(rows):
        pb = rows[atlas]
        lab = uni.get(atlas, [])
        pb50 = {r["cl_id"] for r in pb if int(r["n_cells"]) >= 50}
        lb50 = {r["cl_id"] for r in lab if int(r["n_cells"]) >= 50}
        lines.append(
            f"  {atlas[:44]:44s} {len(pb):>4} "
            f"{sum(1 for r in pb if int(r['n_cells']) < 50):>6} {len(lab):>4} "
            f"{len(lb50):>8} {str(pb50 == lb50):>9}")
    edge("F2a  z universe of the seven new atlases: pseudobulk vs inventory",
         "\n".join(lines))

    import importlib.util

    # Loaded by path rather than by name: the working directory holds non-ASCII
    # characters and the by-name finder has been unreliable for it.
    spec = importlib.util.spec_from_file_location(
        "ldc_01", os.path.join(HERE, "ldc_01_combined_agreement.py"))
    L = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(L)

    onto = L.Ontology(L.OBO)
    _ref, _rank = L.load_reference()
    mapping, _testable = L.load_mapping()
    pb07 = L.load_pseudobulk(P07)
    uni_inv = L.load_universe(LBL07)
    uni_pb = {(a, c): v["n_cells"] for (a, c), v in pb07.items()}
    atlases = {a for a, _ in pb07}
    # h2b z-scored over every pseudobulk row; this run over the inventory.
    stats_pb = L.atlas_stats(pb07, uni_pb, 0)
    stats_inv = L.atlas_stats(pb07, uni_inv, L.FLOOR)
    pairs_pb = L.build_pairs(onto, pb07, mapping, uni_pb, stats_pb, 0, atlases)
    pairs_inv = L.build_pairs(onto, pb07, mapping, uni_inv, stats_inv,
                              L.FLOOR, atlases)
    k_pb = {(p["atlas"], p["hpa_cell_type"]) for p in pairs_pb}
    k_inv = {(p["atlas"], p["hpa_cell_type"]) for p in pairs_inv}
    z_pb = {(p["atlas"], p["hpa_cell_type"]): p["atlas_delta_z"]
            for p in pairs_pb}
    worst = []
    for p in pairs_inv:
        key = (p["atlas"], p["hpa_cell_type"])
        if key in z_pb:
            worst.append((abs(z_pb[key] - p["atlas_delta_z"]), key))
    worst.sort(reverse=True)
    # h2b recorded the *unweighted mean of the member labels' z values*; this
    # run follows b5b and z-scores the library-weighted pooled Delta.  Both are
    # computed here on the same seven atlases.
    h2b_style = {}
    by_atlas = defaultdict(list)
    for (a, c), v in pb07.items():
        by_atlas[a].append((c, v))
    import statistics as st
    zs = {}
    for a, recs in by_atlas.items():
        vals = [v["delta"] for _, v in recs]
        mu, sd = st.fmean(vals), (st.pstdev(vals) or 1.0)
        for c, v in recs:
            zs[(a, c)] = (v["delta"] - mu) / sd
    for p in pairs_inv:
        members = p["atlas_terms"].split(";")
        h2b_style[(p["atlas"], p["hpa_cell_type"])] = \
            st.fmean(zs[(p["atlas"], t)] for t in members)
    h2b_diff = [(abs(h2b_style[k] - p["atlas_delta_z"]), k) for p in pairs_inv
                if (k := (p["atlas"], p["hpa_cell_type"])) in h2b_style]
    h2b_diff.sort(reverse=True)
    frozen7 = list(csv.DictReader(open(H2B_PAIRS, encoding="utf-8"),
                                  delimiter="\t"))
    shared7 = [(p, r) for p in pairs_inv for r in frozen7
               if (p["atlas"], p["hpa_cell_type"])
               == (r["atlas"], r["reference_cell_type"])]
    rho_pooled = L.spearman([_ref[p["hpa_cell_type"]] for p, _ in shared7],
                            [p["atlas_delta_z"] for p, _ in shared7])
    rho_h2bstyle = L.spearman([_ref[p["hpa_cell_type"]] for p, _ in shared7],
                              [h2b_style[(p["atlas"], p["hpa_cell_type"])]
                               for p, _ in shared7])
    rho_frozen = L.spearman([float(r["reference_delta"]) for _, r in shared7],
                            [float(r["atlas_delta_z"]) for _, r in shared7])
    edge("F2b  the seven new atlases: pseudobulk universe vs inventory "
         "universe",
         f"pairs, pseudobulk universe : {len(k_pb)}\n"
         f"pairs, inventory universe  : {len(k_inv)}\n"
         f"h2b published              : {len(open(H2B_PAIRS, encoding='utf-8').readlines()) - 1}\n"
         f"rows only under the inventory universe:\n"
         + "\n".join(f"    {a}  {b}" for a, b in sorted(k_inv - k_pb))
         + "\n  largest z differences (both universes):\n"
         + "\n".join(f"    {d:.4f}  {a}  {b}" for d, (a, b) in worst[:5]))
    edge("F2c  pooling rule, not the z universe, drives the h2b gap",
         f"unweighted mean of member z (h2b rule) : "
         f"rho = {rho_h2bstyle:.4f}\n"
         f"z of library-weighted pooled Delta     : rho = {rho_pooled:.4f}\n"
         f"h2b frozen table, as published         : rho = {rho_frozen:.4f}\n"
         f"largest |difference| between the two pooling rules: "
         f"{h2b_diff[0][0]:.4f} on {h2b_diff[0][1][1]} / "
         f"{h2b_diff[0][1][0][:28]}")


if __name__ == "__main__":
    main()
