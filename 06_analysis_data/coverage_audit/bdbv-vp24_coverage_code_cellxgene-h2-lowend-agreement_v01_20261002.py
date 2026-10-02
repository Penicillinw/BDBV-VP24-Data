"""H2b: hierarchy-aware version of the low-end cross-atlas test.

The first pass of this test joined atlas labels to the reference through exact
Cell Ontology identifiers, so the urothelial and photoreceptor datasets failed to
match even though they contain the compartments in question: their labels are
subtypes (basal / intermediate / umbrella cell of urothelium; foveal and
peripheral photoreceptor) whose Cell Ontology terms sit BELOW the reference
terms. This script repeats the test with the same ancestry rule the original
cross-atlas test used (is_a and part_of both treated as ancestry), so a
reference term matches an atlas label when the atlas label is the term or one of
its descendants.

Everything else is held at the values of the frozen rule set: raw counts, CPM,
the three Delta components as the mean of log10(CPM+1), z-scoring within atlas,
a cell-count floor of 50, and a label-permutation p value.

Outputs: h2b_agreement.json, h2b_pairs.tsv
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
PB = os.path.join(HERE, "raw", "new_atlas_pseudobulk.tsv")
LAB = os.path.join(HERE, "raw", "new_atlas_labels.tsv")
LABELS = os.path.join(HERE, "raw", "new_atlas_pseudobulk.tsv")
FLOOR = 50


def load_ancestors():
    parents = {}
    cur = None
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


def perm_p(x, y, draws=4000, seed=20260926):
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


def main():
    parents = load_ancestors()
    memo = {}
    ref = {}
    for r in csv.DictReader(open(REF, encoding="utf-8"), delimiter="\t"):
        ref[r["cell_type"]] = float(r["delta_restriction"])
    term_to_hpa = {}
    for r in csv.DictReader(open(MAP, encoding="utf-8"), delimiter="\t"):
        term_to_hpa.setdefault(r["cl_id"], []).append(r["hpa_cell_type"])
    lab = {(r["atlas"], r["cl_id"]): r["cl_label"]
           for r in csv.DictReader(open(LAB, encoding="utf-8"), delimiter="\t")}
    rows = list(csv.DictReader(open(PB, encoding="utf-8"), delimiter="\t"))
    by_atlas = {}
    for r in rows:
        by_atlas.setdefault(r["atlas"], []).append(r)

    per_atlas, pairs = [], []
    for atlas, recs in sorted(by_atlas.items()):
        vals = [float(r["delta"]) for r in recs]
        mu, sd = st.fmean(vals), (st.pstdev(vals) or 1.0)
        z = {r["cl_id"]: (float(r["delta"]) - mu) / sd for r in recs}
        matches = {}
        for r in recs:
            if int(r["n_cells"]) < FLOOR:
                continue
            cl = r["cl_id"]
            anc = closure(cl, parents, memo) if cl in parents or cl.startswith("CL:") else {cl}
            for t in anc:
                for hpa in term_to_hpa.get(t, []):
                    if hpa in ref:
                        matches.setdefault(hpa, []).append((z[cl], int(r["n_cells"]),
                                                            lab.get((atlas, cl), cl)))
        if len(matches) < 3:
            per_atlas.append({"atlas": atlas, "n_matched": len(matches), "rho": None})
            continue
        names = sorted(matches)
        xs = [ref[h] for h in names]
        ys = [st.fmean(v[0] for v in matches[h]) for h in names]
        rho = spearman(xs, ys)
        per_atlas.append({"atlas": atlas, "n_matched": len(names),
                          "rho": None if math.isnan(rho) else round(rho, 3)})
        for h in names:
            v = matches[h]
            pairs.append({"atlas": atlas, "reference_cell_type": h,
                          "reference_delta": round(ref[h], 3),
                          "atlas_delta_z": round(st.fmean(x[0] for x in v), 3),
                          "atlas_labels": " | ".join(sorted({x[2] for x in v})),
                          "n_cells": sum(x[1] for x in v)})
    xs = [p["reference_delta"] for p in pairs]
    ys = [p["atlas_delta_z"] for p in pairs]
    pooled = {"n": len(pairs), "rho": round(spearman(xs, ys), 3),
              "p_permutation": round(perm_p(xs, ys), 4)} if len(pairs) > 2 else None
    out = {"rule": "hierarchy-aware (is_a + part_of ancestry), floor 50 cells",
           "per_atlas": per_atlas, "pooled": pooled,
           "pairs": sorted(pairs, key=lambda p: p["reference_delta"])}
    with open(os.path.join(HERE, "h2b_agreement.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(HERE, "h2b_pairs.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["atlas", "reference_cell_type",
                                       "reference_delta", "atlas_delta_z",
                                       "atlas_labels", "n_cells"])
        w.writeheader()
        w.writerows(out["pairs"])
    print(json.dumps({"per_atlas": per_atlas, "pooled": pooled},
                     ensure_ascii=False, indent=1))
    print("\npairs:")
    for p in out["pairs"]:
        print(f"  {p['reference_cell_type']:40s} ref {p['reference_delta']:+7.3f} "
              f"atlas {p['atlas_delta_z']:+7.3f} n={p['n_cells']:>7d} "
              f"[{p['atlas'][:24]}] {p['atlas_labels'][:44]}")


if __name__ == "__main__":
    main()
