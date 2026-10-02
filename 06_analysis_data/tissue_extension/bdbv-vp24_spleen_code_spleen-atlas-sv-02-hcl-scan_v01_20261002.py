"""SRB-VERIFY step 2: independent single-pass re-derivation of the HCL
pseudobulk values (criteria V2/V3 inputs).

Independence rule: no import from analysis/spleen_reference_rebuild_20260926/,
no use of its pseudobulk products as inputs.  Everything is recomputed from the
read-only HCL atlas (analysis/h2_lowend_atlas_hunt_20260926/raw/atlas_h5ad/
2adb1f8a_...h5ad) with an independently written CSR accumulator.

Two groupings are built in the same pass:
  A  whole-atlas pseudobulk by `cell_type_ontology_term_id`   (anchor side)
  B  spleen-restricted pseudobulk by `cell_type_ontology_term_id`

Also regenerates the 200 expression-matched six-gene modules from scratch
(same frozen rule: each gene drawn from HPA genes whose mean log10(nCPM+1)
across the 154 reference types is within +/-0.10 of its cognate score gene,
seed 20260926) and compares them with SRB's raw/hpa_gene_background_matched.tsv.

Writes:
    raw/sv_hcl_whole_atlas_pseudobulk.tsv
    raw/sv_hcl_spleen_pseudobulk.tsv
    raw/sv_hcl_scan.json
    raw/sv_modules.json
"""

from __future__ import annotations

import bisect
import csv
import json
import math
import os
import random

import h5py
import numpy as np
from scipy import sparse

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUT = os.path.join(ROOT, "analysis", "srb_verify_20260926")
RAW = os.path.join(OUT, "raw")
SRB_RAW = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926", "raw")
H2RAW = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926", "raw")
HCL = os.path.join(H2RAW, "atlas_h5ad",
                   "2adb1f8a_Construction_of_a_human_cell_landscape_at_single.h5ad")

SCORE_GENES = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]
N_MODULES = 200
WINDOW = 0.10
SEED = 20260926
BLOCK = 20000
FLOOR = 50


def decode_group(node):
    """Decode an AnnData categorical group or a plain dataset to object array."""
    if isinstance(node, h5py.Group):
        codes = node["codes"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c)
                         for c in node["categories"][:]], dtype=object)
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = ""
        return out
    arr = node[:]
    return np.array([a.decode() if isinstance(a, bytes) else str(a) for a in arr], dtype=object)


def build_modules(bg):
    """Frozen rule regeneration of the 200 expression-matched six-gene modules."""
    rng = random.Random(SEED)
    targets = [bg[g] for g in SCORE_GENES]
    pool = sorted(((g, v) for g, v in bg.items() if g not in SCORE_GENES), key=lambda kv: kv[1])
    levels = [v for _, v in pool]
    modules, used, tries = [], set(), 0
    while len(modules) < N_MODULES and tries < N_MODULES * 200:
        tries += 1
        m = []
        for t in targets:
            lo = bisect.bisect_left(levels, t - WINDOW)
            hi = bisect.bisect_right(levels, t + WINDOW)
            if hi <= lo:
                continue
            for _ in range(50):
                g = pool[rng.randrange(lo, hi)][0]
                if g not in m and g not in used:
                    m.append(g)
                    break
        if len(m) == 6:
            modules.append(m)
            used.update(m)
    return modules, targets, tries


def main() -> None:
    # ---- regenerate modules from my own background table (written by sv_01) ----
    bg = {}
    with open(os.path.join(RAW, "sv_hpa_gene_background.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bg[r["gene"]] = float(r["mean_log10_nCPM"])
    print(f"background genes (own re-derivation): {len(bg)}")
    modules, targets, tries = build_modules(bg)
    print(f"modules built: {len(modules)} (tries {tries})")

    with open(os.path.join(SRB_RAW, "hpa_gene_background_matched.tsv"), encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    per_idx = {}
    for r in rows:
        per_idx.setdefault(int(r["module_index"]), []).append(r["matched_gene"])
    srb_mod_list = [per_idx[i] for i in sorted(per_idx)]
    same_set = (len(srb_mod_list) == len(modules)
                and all(sorted(a) == sorted(b) for a, b in zip(modules, srb_mod_list)))
    n_overlap = sum(1 for a, b in zip(modules, srb_mod_list) if sorted(a) == sorted(b))
    print(f"module list agreement with SRB: exact-set {same_set} ; identical module objects {n_overlap}/{len(modules)}")

    panel = SCORE_GENES + [g for m in modules for g in m]
    with open(os.path.join(RAW, "sv_modules.json"), "w", encoding="utf-8") as fh:
        json.dump({"n_modules": len(modules), "seed": SEED, "window": WINDOW,
                   "gene_means": dict(zip(SCORE_GENES, targets)),
                   "identical_module_objects": n_overlap,
                   "modules": modules,
                   "srb_modules": srb_mod_list}, fh, indent=1)

    # ---- single pass over the atlas ----
    with h5py.File(HCL, "r") as f:
        genes = [str(x) for x in decode_group(f["raw/var"]["feature_name"])]
        tissue = decode_group(f["obs"]["tissue"])
        cl = decode_group(f["obs"]["cell_type_ontology_term_id"])
        node = f["raw/X"]
        n_cells, n_genes = (int(v) for v in node.attrs["shape"])
        assert len(genes) == n_genes and len(tissue) == n_cells

        first = {}
        for i, g in enumerate(genes):
            first.setdefault(g, i)
        pos = np.array([first.get(g, -1) for g in panel], dtype=np.int64)
        keep = np.flatnonzero(pos >= 0)
        cols = pos[keep]
        missing = [panel[k] for k in range(len(panel)) if pos[k] < 0]
        col_of = {panel[k]: j for j, k in enumerate(keep)}
        print(f"panel {len(panel)} genes -> {len(keep)} found in HCL, {len(missing)} missing")

        lut = np.full(n_genes, -1, dtype=np.int64)
        lut[cols] = np.arange(len(keep), dtype=np.int64)

        labels_a, invA = np.unique(cl, return_inverse=True)
        spleen_mask = (tissue == "spleen")
        print(f"spleen cells by tissue label: {int(spleen_mask.sum())}")
        labels_b, invB = np.unique(cl[spleen_mask], return_inverse=True)
        nA, nB = len(labels_a), len(labels_b)
        invB_full = np.full(n_cells, -1, dtype=np.int64)
        invB_full[spleen_mask] = invB

        accA = np.zeros((nA, len(keep)), dtype=np.float64)
        accB = np.zeros((nB, len(keep)), dtype=np.float64)
        libA = np.zeros(nA, dtype=np.float64)
        libB = np.zeros(nB, dtype=np.float64)

        onehotA = sparse.csr_matrix((np.ones(n_cells), (np.arange(n_cells), invA)),
                                    shape=(n_cells, nA))
        data, indices, indptr = node["data"], node["indices"], node["indptr"]
        for i0 in range(0, n_cells, BLOCK):
            i1 = min(n_cells, i0 + BLOCK)
            ip = indptr[i0:i1 + 1]
            d0, d1 = int(ip[0]), int(ip[-1])
            if d1 <= d0:
                continue
            d = data[d0:d1].astype(np.float64)
            ix = indices[d0:d1]
            rel = ip - d0
            rows = np.repeat(np.arange(i1 - i0, dtype=np.int64), np.diff(rel))
            # library sizes (all genes) per label
            libA += np.bincount(invA[i0:i1][rows], weights=d, minlength=nA)
            gb = invB_full[i0:i1][rows]
            okb = gb >= 0
            if okb.any():
                libB += np.bincount(gb[okb], weights=d[okb], minlength=nB)
            # selected genes only
            cg = lut[ix]
            m = cg >= 0
            if m.any():
                sub = sparse.csr_matrix((d[m], (rows[m], cg[m])),
                                        shape=(i1 - i0, len(keep)))
                accA += np.asarray((onehotA[i0:i1].T @ sub).todense())
                sb = invB_full[i0:i1]
                okr = sb >= 0
                if okr.any():
                    onehotB = sparse.csr_matrix(
                        (np.ones(int(okr.sum())), (np.flatnonzero(okr), sb[okr])),
                        shape=(i1 - i0, nB))
                    accB += np.asarray((onehotB.T @ sub).todense())

    countsA = np.bincount(invA, minlength=nA)
    countsB = np.bincount(invB, minlength=nB)

    def cpm_table(acc, lib, labels, gene_keys, counts):
        out = []
        for i, lab in enumerate(labels):
            row = {"label": lab, "n_cells": int(counts[i])}
            for g in gene_keys:
                row[g] = 1e6 * acc[i, col_of[g]] / lib[i] if lib[i] > 0 else float("nan")
            out.append(row)
        return out

    rowsA = cpm_table(accA, libA, labels_a, SCORE_GENES, countsA)
    with open(os.path.join(RAW, "sv_hcl_whole_atlas_pseudobulk.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        fh.write("label\tn_cells\t" + "\t".join(SCORE_GENES) + "\n")
        for r in rowsA:
            fh.write(r["label"] + "\t" + str(r["n_cells"]) + "\t"
                     + "\t".join(f"{r[g]:.4f}" for g in SCORE_GENES) + "\n")

    rowsB = cpm_table(accB, libB, labels_b, SCORE_GENES, countsB)
    with open(os.path.join(RAW, "sv_hcl_spleen_pseudobulk.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        fh.write("label\tn_cells\t" + "\t".join(SCORE_GENES) + "\n")
        for r in rowsB:
            fh.write(r["label"] + "\t" + str(r["n_cells"]) + "\t"
                     + "\t".join(f"{r[g]:.4f}" for g in SCORE_GENES) + "\n")

    np.savez_compressed(os.path.join(RAW, "sv_hcl_accumulators.npz"),
                        accA=accA, accB=accB, libA=libA, libB=libB,
                        labels_a=np.array(labels_a, dtype=object),
                        labels_b=np.array(labels_b, dtype=object),
                        panel=np.array([panel[k] for k in keep], dtype=object),
                        keep=np.array(keep))

    # raw integer check + spleen row-block equivalence
    with h5py.File(HCL, "r") as f:
        node = f["raw/X"]
        d0 = node["data"][:200000]
    intfrac = float(np.mean(d0 == np.round(d0)))

    spleen_idx = np.flatnonzero(spleen_mask)
    info = {
        "n_cells": n_cells, "n_genes": n_genes,
        "spleen_cells_by_label": int(spleen_mask.sum()),
        "spleen_row_block_first_last": [int(spleen_idx[0]), int(spleen_idx[-1])],
        "spleen_is_one_contiguous_block": bool(spleen_idx[-1] - spleen_idx[0] + 1 == len(spleen_idx)),
        "spleen_block_matches_SRB_range_191895_207700": bool(
            spleen_idx[0] == 191895 and spleen_idx[-1] == 207700),
        "n_whole_atlas_labels": nA, "n_spleen_labels": nB,
        "raw_integral_fraction_first_200k": intfrac,
        "panel_genes_requested": len(panel), "panel_genes_found": len(keep),
        "panel_genes_missing": missing,
        "spleen_labels_ge_floor": int(sum(1 for r in rowsB if r["n_cells"] >= FLOOR)),
        "spleen_labels_all": [{"label": r["label"], "n_cells": r["n_cells"]} for r in rowsB],
        "module_objects_identical_to_SRB": n_overlap,
    }
    with open(os.path.join(RAW, "sv_hcl_scan.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in info.items() if k != "spleen_labels_all"},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
