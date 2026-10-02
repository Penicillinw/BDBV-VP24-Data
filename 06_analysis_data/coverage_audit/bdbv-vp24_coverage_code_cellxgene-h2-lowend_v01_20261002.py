"""H2: does the Delta ordering replicate in independent atlases that DO contain the
low-end compartments the earlier cross-atlas test could not cover?

Background
----------
The manuscript's partition places retinal pigment epithelium, photoreceptors,
Sertoli and Leydig cells and choroid plexus epithelium at the low end of the
interferon-tone axis, and urothelium at the high end. The nineteen atlases used
for the original cross-atlas test contained none of those cell types, so those
predictions were untested rather than contradicted. Seven further CELLxGENE
datasets that do contain them were downloaded (see raw/download_manifest.json);
their pseudobulk scores were computed by h2_pseudobulk_new_atlases.py using the
same conventions as the original test.

Frozen rules
------------
R1  the atlas score is Delta from raw counts, z-scored across the labels of that
    atlas (the primary convention of the original cross-atlas test).
R2  labels are joined through Cell Ontology identifiers, not names: a CL id is
    matched to the reference through the mapping table of the original test.
R3  the primary cell-count floor is 50 cells per label; the floor of 1 is
    reported as a sensitivity analysis, not as the primary result.
R4  a compartment is called replicated only if the matched pair points the same
    way as the reference (high reference Delta with high atlas Delta, or low with
    low); the pooled correlation is reported with its n and a label permutation
    p value.
"""

from __future__ import annotations

import csv
import json
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
REF = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                   "restriction_gradient.tsv")
MAP = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926",
                   "label_mapping.tsv")
PB = os.path.join(HERE, "raw", "new_atlas_pseudobulk.tsv")
FLOOR = 50


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1
        return r
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = st.fmean(rx), st.fmean(ry)
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float("nan")


def perm_p(x, y, draws=4000, seed=20260926):
    obs = abs(spearman(x, y))
    rng = seed
    yy = list(y)
    n = len(yy)
    ge = 0
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
    ref = {}
    for r in csv.DictReader(open(REF, encoding="utf-8"), delimiter="\t"):
        ref[r["cell_type"]] = float(r["delta_restriction"])
    cl_to_hpa = {}
    for r in csv.DictReader(open(MAP, encoding="utf-8"), delimiter="\t"):
        cl_to_hpa.setdefault(r["cl_id"], []).append(r["hpa_cell_type"])

    rows = list(csv.DictReader(open(PB, encoding="utf-8"), delimiter="\t"))
    by_atlas = {}
    for r in rows:
        by_atlas.setdefault(r["atlas"], []).append(r)

    per_atlas, pairs, rejected = [], [], {}
    for atlas, recs in sorted(by_atlas.items()):
        vals = [float(r["delta"]) for r in recs]
        mu, sd = st.fmean(vals), (st.pstdev(vals) or 1.0)
        z = {r["cl_id"]: (float(r["delta"]) - mu) / sd for r in recs}
        matched = []
        for r in recs:
            if int(r["n_cells"]) < FLOOR:
                rejected.setdefault(atlas, []).append(
                    {"cl_id": r["cl_id"], "n_cells": int(r["n_cells"])})
                continue
            for hpa in cl_to_hpa.get(r["cl_id"], []):
                if hpa in ref:
                    matched.append((hpa, ref[hpa], z[r["cl_id"]],
                                    int(r["n_cells"])))
        if not matched:
            per_atlas.append({"atlas": atlas, "n_matched": 0, "rho": None})
            continue
        rho = spearman([m[1] for m in matched], [m[2] for m in matched]) \
            if len(matched) > 2 else None
        per_atlas.append({"atlas": atlas, "n_matched": len(matched),
                          "rho": None if rho is None or math.isnan(rho) else round(rho, 3)})
        for m in matched:
            pairs.append({"atlas": atlas, "cell_type": m[0],
                          "reference_delta": round(m[1], 3),
                          "atlas_delta_z": round(m[2], 3), "n_cells": m[3]})

    pools = []
    for a in per_atlas:
        if a["n_matched"]:
            pools.extend([p for p in pairs if p["atlas"] == a["atlas"]])
    pooled = None
    if len(pools) > 2:
        pooled = {"n": len(pools),
                  "rho": round(spearman([p["reference_delta"] for p in pools],
                                        [p["atlas_delta_z"] for p in pools]), 3)}
        pooled["p_permutation"] = round(
            perm_p([p["reference_delta"] for p in pools],
                   [p["atlas_delta_z"] for p in pools]), 4)
    out = {"primary_floor_cells": FLOOR,
           "per_atlas": per_atlas,
           "pooled": pooled,
           "pairs": sorted(pairs, key=lambda p: p["reference_delta"]),
           "labels_below_floor": rejected}
    with open(os.path.join(HERE, "h2_lowend_agreement.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: out[k] for k in ("primary_floor_cells", "per_atlas", "pooled")},
                     ensure_ascii=False, indent=1))
    print("\nmatched pairs (reference Delta, atlas Delta z):")
    for p in out["pairs"]:
        print(f"  {p['cell_type']:42s} ref {p['reference_delta']:+7.3f} "
              f"atlas {p['atlas_delta_z']:+7.3f}  n={p['n_cells']:>7d}  "
              f"[{p['atlas'][:26]}]")


if __name__ == "__main__":
    main()
