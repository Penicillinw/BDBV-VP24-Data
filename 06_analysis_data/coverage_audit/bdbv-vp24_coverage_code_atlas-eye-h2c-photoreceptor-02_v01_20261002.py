"""H2c: close the photoreceptor half of the low-end coverage gap.

Diagnosis first (read-only): the HPA cell types "rod photoreceptor cells" and
"cone photoreceptor cells" were already present in the B5 candidate table
(`analysis/b5_harmonised_atlas_20260926/raw/hpa_cl_candidates.tsv`) mapped to
CL:0000604 (retinal rod cell) and CL:0000573 (retinal cone cell) by
token_overlap -- but with `n_atlas_terms_covered = 0`, i.e. none of the 19
original atlases contained a cell carrying those ontology terms.  They were
therefore untestable and never entered the curated `label_mapping.tsv`, which is
why the hierarchy-aware pass (`h2_lowend_agreement_hierarchy.py`) could not
connect them: the failure is on the REFERENCE side of the join, not on the atlas
side.

This script supplies the two mapping rows (see `extension_mapping.tsv`) and
repeats the frozen hierarchy-aware cross-atlas test with everything else held at
the frozen values: raw counts, CPM, the three Delta components as the mean of
log10(CPM+1), z-scoring within atlas, a cell-count floor of 50, and a
label-permutation p value (4000 draws, seed 20260926).

Both the frozen run (n = 38 pairs, rho = 0.568) and the extended run are
reported, plus a leave-one-atlas-out sensitivity for the two eye atlases that
now contribute photoreceptors.

Outputs: h2c_pairs.tsv, h2c_agreement.json
"""

from __future__ import annotations

import csv
import json
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OBO = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926", "raw", "cl.obo")
REF = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                   "restriction_gradient.tsv")
MAP = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926",
                   "label_mapping.tsv")
EXT = os.path.join(HERE, "extension_mapping.tsv")
H2 = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926")
PB = os.path.join(H2, "raw", "new_atlas_pseudobulk.tsv")
LAB = os.path.join(H2, "raw", "new_atlas_labels.tsv")
FLOOR = 50
DRAWS = 4000
SEED = 20260926
GENES = {"IFN_I": ["IFNAR1", "IFNAR2"], "IFN_III": ["IFNLR1", "IL10RB"],
         "ISG": ["ISG15", "MX1"]}


def load_ancestors():
    parents, cur = {}, None
    with open(OBO, encoding="utf-8") as fh:
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


def closure(term, parents, memo):
    if term in memo:
        return memo[term]
    out = {term}
    for p in parents.get(term, ()):
        out |= closure(p, parents, memo)
    memo[term] = out
    return out


def spearman(x, y):
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


def perm_p(x, y, draws=DRAWS, seed=SEED):
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


def component_means(row):
    """Recompute the three Delta components from the CPM columns, log10(CPM+1)."""
    vals = {}
    for name, genes in GENES.items():
        vals[name] = st.fmean(math.log10(float(row["cpm_" + g]) + 1.0)
                              for g in genes)
    return vals


def main():
    parents = load_ancestors()
    memo = {}
    ref = {r["cell_type"]: float(r["delta_restriction"])
           for r in csv.DictReader(open(REF, encoding="utf-8"), delimiter="\t")}

    mapping_rows = list(csv.DictReader(open(MAP, encoding="utf-8"), delimiter="\t"))
    ext_rows = list(csv.DictReader(open(EXT, encoding="utf-8"), delimiter="\t"))
    frozen_map = {}
    for r in mapping_rows:
        frozen_map.setdefault(r["cl_id"], []).append(r["hpa_cell_type"])
    ext_map = {}
    for r in ext_rows:
        ext_map.setdefault(r["cl_id"], []).append(r["hpa_cell_type"])
    unmapped_hpa = sorted(set(ref) - {r["hpa_cell_type"] for r in mapping_rows})

    lab = {(r["atlas"], r["cl_id"]): r["cl_label"]
           for r in csv.DictReader(open(LAB, encoding="utf-8"), delimiter="\t")}
    rows = list(csv.DictReader(open(PB, encoding="utf-8"), delimiter="\t"))

    # integrity check: stored delta vs a recompute from the CPM columns
    worst = 0.0
    for r in rows:
        recomputed = st.fmean(component_means(r).values())
        worst = max(worst, abs(recomputed - float(r["delta"])))

    by_atlas = {}
    for r in rows:
        by_atlas.setdefault(r["atlas"], []).append(r)

    def run(term_to_hpa, min_matches):
        """min_matches reproduces the frozen rule of the hierarchy-aware pass:
        an atlas enters the pooled test only if it contributes at least three
        matched labels (a within-atlas correlation is otherwise not defined)."""
        pairs = []
        skipped = {}
        for atlas, recs in sorted(by_atlas.items()):
            vals = [float(r["delta"]) for r in recs]
            mu, sd = st.fmean(vals), (st.pstdev(vals) or 1.0)
            z = {r["cl_id"]: (float(r["delta"]) - mu) / sd for r in recs}
            matches = {}
            for r in recs:
                if int(r["n_cells"]) < FLOOR:
                    continue
                cl = r["cl_id"]
                anc = closure(cl, parents, memo) if cl.startswith("CL:") else {cl}
                for t in anc:
                    for hpa in term_to_hpa.get(t, []):
                        if hpa in ref:
                            matches.setdefault(hpa, []).append(
                                (z[cl], int(r["n_cells"]), lab.get((atlas, cl), cl)))
            if len(matches) < min_matches:
                skipped[atlas] = {h: round(st.fmean(x[0] for x in v), 3)
                                  for h, v in sorted(matches.items())}
                continue
            for h in sorted(matches):
                v = matches[h]
                pairs.append({"atlas": atlas, "reference_cell_type": h,
                              "reference_delta": round(ref[h], 3),
                              "atlas_delta_z": round(st.fmean(x[0] for x in v), 3),
                              "atlas_labels": " | ".join(sorted({x[2] for x in v})),
                              "n_cells": sum(x[1] for x in v)})
        xs = [p["reference_delta"] for p in pairs]
        ys = [p["atlas_delta_z"] for p in pairs]
        summary = {"n": len(pairs), "n_atlases": len({p["atlas"] for p in pairs}),
                   "rho": round(spearman(xs, ys), 3),
                   "p_permutation": round(perm_p(xs, ys), 4),
                   "min_matches_per_atlas": min_matches,
                   "skipped_atlases": skipped}
        return pairs, summary

    frozen_pairs, frozen = run(frozen_map, 3)
    _, frozen_all = run(frozen_map, 0)
    merged = {k: list(v) for k, v in frozen_map.items()}
    for k, v in ext_map.items():
        merged.setdefault(k, []).extend(v)
    ext_pairs, extended = run(merged, 3)
    _, extended_all = run(merged, 0)

    # self-check: the frozen rule set must reproduce the recorded hierarchy pass
    recorded = {"n": 38, "rho": 0.568, "p_permutation": 0.0005}
    self_check = {
        "recorded_h2b": recorded,
        "reproduced": {k: frozen[k] for k in recorded},
        "matches": all(frozen[k] == v for k, v in recorded.items()),
    }

    added = [p for p in ext_pairs
             if p not in frozen_pairs]

    # leave-one-atlas-out on the atlases that now contribute photoreceptors
    loo = {}
    target_atlases = sorted({p["atlas"] for p in added})
    for atlas in target_atlases:
        keep = [p for p in ext_pairs if p["atlas"] != atlas]
        xs = [p["reference_delta"] for p in keep]
        ys = [p["atlas_delta_z"] for p in keep]
        if len(keep) > 2:
            loo[atlas] = {"n": len(keep), "rho": round(spearman(xs, ys), 3),
                          "p_permutation": round(perm_p(xs, ys), 4)}
    no_photo_atlas = [p for p in ext_pairs
                      if p["atlas"] != "7b75b2c4_Photoreceptor_cells_of_the_human_fovea_and_perip.h5ad"]
    xs = [p["reference_delta"] for p in no_photo_atlas]
    ys = [p["atlas_delta_z"] for p in no_photo_atlas]
    drop_fovea = {"n": len(no_photo_atlas), "rho": round(spearman(xs, ys), 3),
                  "p_permutation": round(perm_p(xs, ys), 4)}

    out = {
        "diagnosis": {
            "unmapped_reference_cell_types": len(unmapped_hpa),
            "note": ("photoreceptor rows were candidates in the B5 candidate table "
                     "but had zero atlas coverage at curation time"),
        },
        "integrity": {"max_abs_deviation_stored_delta_vs_recomputed": round(worst, 12)},
        "self_check": self_check,
        "frozen_run": frozen,
        "frozen_run_all_atlases": frozen_all,
        "extended_run": extended,
        "extended_run_all_atlases": extended_all,
        "added_pairs": sorted(added, key=lambda p: p["reference_delta"]),
        "leave_one_atlas_out": loo,
        "drop_fovea_atlas": drop_fovea,
        "pairs": sorted(ext_pairs, key=lambda p: p["reference_delta"]),
    }
    with open(os.path.join(HERE, "h2c_agreement.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(HERE, "h2c_pairs.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["atlas", "reference_cell_type",
                                       "reference_delta", "atlas_delta_z",
                                       "atlas_labels", "n_cells"])
        w.writeheader()
        w.writerows(out["pairs"])

    print(json.dumps({k: out[k] for k in
                      ("self_check", "frozen_run", "extended_run",
                       "frozen_run_all_atlases", "extended_run_all_atlases",
                       "leave_one_atlas_out", "drop_fovea_atlas", "integrity")},
                     ensure_ascii=False, indent=1))
    print("\nadded pairs:")
    for p in out["added_pairs"]:
        print(f"  {p['reference_cell_type']:34s} ref {p['reference_delta']:+7.3f} "
              f"atlas {p['atlas_delta_z']:+7.3f} n={p['n_cells']:>7d} "
              f"[{p['atlas'][:22]}] {p['atlas_labels'][:34]}")


if __name__ == "__main__":
    main()
