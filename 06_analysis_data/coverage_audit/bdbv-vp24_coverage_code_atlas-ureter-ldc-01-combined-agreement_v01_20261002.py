"""LDC-01: one-convention cross-atlas agreement over all 26 cached atlases.

Why this exists
---------------
The manuscript's cross-atlas reproduction number is computed on the nineteen
atlases cached first (b5 / b5b).  Those nineteen contain no cell from the low
end of the Delta ranking (urothelium, retinal pigment epithelium,
photoreceptors, Sertoli, choroid plexus).  H2 downloaded seven further atlases
that do cover those compartments and reported a *separate* pooled correlation
(rho = 0.568, n = 38) over the new atlases only.  Two numbers computed on two
different atlas sets cannot be compared, so it was not known whether adding the
low-end compartments helps or hurts the original claim.  The low-end coverage
scan (`analysis/lowend_reference_scan_20260926`, section 9 item 3) left exactly
this run to the main agent; this script is that run.

Frozen criteria, written before any number in this directory was computed
------------------------------------------------------------------------
R1  Inputs are the two frozen pseudobulk tables, and nothing else is read for
    values:
      P19  analysis/b5_harmonised_atlas_20260926/raw/atlas_pseudobulk.tsv
      P07  analysis/h2_lowend_atlas_hunt_20260926/raw/new_atlas_pseudobulk.tsv
    Column `delta` is the atlas-side Delta in both.  No .h5ad file is re-read.
R2  Ontology: analysis/b5_harmonised_atlas_20260926/raw/cl.obo, parsed by the
    b5b parser, i.e. the `is_a` qualifier `{...}` is stripped as well as the
    trailing comment.  Ancestry is is_a together with part_of, reflexive and
    transitive.  (The original b5 parse kept the qualifier, produced ids such
    as 'CL:0000034 {is_inferred="true"}' and silently dropped those edges;
    b5b corrected that and the corrected parse is used here.)
R3  Reference: analysis/t4_ifn_landscape_20260925/restriction_gradient.tsv,
    column `delta_restriction`, all 154 types.
R4  Mapping: analysis/b5_harmonised_atlas_20260926/label_mapping.tsv.  Two
    reference conventions are carried side by side and both are reported:
      A  only rows with `testable == yes` -- the convention under which the
         manuscript's 162-pair number was produced;
      B  every mapped row -- ignores `testable`, because the low-end scan
         showed that column is derived with the ancestry direction inverted
         (it counts stricter-than-reference atlas terms, so 79 rows are
         flagged but only 39 ever produced a pair under the frozen pipeline).
R5  A pair (h, a) exists iff at least one atlas term t of atlas a has
    n_cells >= FLOOR AND h's CL id lies in ancestors(t) -- i.e. the atlas term
    is the reference term or something more specific under it.  The atlas value
    of the pair is the library-weighted pool of the member labels' CPM values
    converted to Delta by the same three components as everywhere else.
R6  Standardisation.  x = the reference Delta as tabulated (no strictly
    monotone transform of x can change a pooled Spearman rho, so the raw column
    is used and the choice is not free-floating).  y = the pair's pooled Delta
    z-scored within its own atlas over that atlas's labels with n_cells >=
    FLOOR.  FLOOR = 50 is the frozen primary.
R7  Statistics: Spearman on mid-ranks; label permutation over the paired
    reference values, seed 20260926, 20000 draws when n <= 60 and 4000 draws
    when n > 60 (the b5b rule); the Student-t approximation reported beside it;
    the smallest |rho| detectable at 80 per cent power, and the achieved power.
R8  Self-checks that must pass before the extension is read:
      S1  19 atlases, convention A  ->  n = 162 and rho = 0.5191  (b5b)
      S2   7 atlases, convention B  ->  n =  38 and rho = 0.568   (h2b)
    S1 is additionally checked field by field against raw/b5_matched_pairs.tsv
    (`y_all_labels`), and S2 against h2b_pairs.tsv (`atlas_delta_z`).
R9  Extension: the 26 atlases together, under A and under B; per atlas;
    leave-one-atlas-out; and the pairs split into the frozen low end (reference
    Delta rank >= 116, or one of the six compartments H2 targeted) against the
    rest.

Outputs: ldc_summary.json, ldc_agreement.tsv, ldc_pairs.tsv
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
B5B = os.path.join(ROOT, "analysis", "b5b_obo_fix_20260926")
H2 = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926")
sys.path.insert(0, B5B)

from b5_common import (  # noqa: E402
    Ontology,
    achieved_power,
    hpa_delta,
    load_hpa,
    min_detectable_rho,
    spearman,
    spearman_p_value,
)

P19 = os.path.join(B5, "raw", "atlas_pseudobulk.tsv")
P07 = os.path.join(H2, "raw", "new_atlas_pseudobulk.tsv")
INV19 = os.path.join(B5, "raw", "atlas_celltype_inventory.tsv")
LBL07 = os.path.join(H2, "raw", "new_atlas_labels.tsv")
OBO = os.path.join(B5, "raw", "cl.obo")
REF = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                   "restriction_gradient.tsv")
MAP = os.path.join(B5, "label_mapping.tsv")
B5B_PAIRS = os.path.join(B5B, "b5_matched_pairs.tsv")
H2B_PAIRS = os.path.join(H2, "h2b_pairs.tsv")

GENES = ("IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1")
MODULE = {"IFN_I_capacity": ["IFNAR1", "IFNAR2"],
          "IFN_III_capacity": ["IFNLR1", "IL10RB"],
          "ISG_priming": ["ISG15", "MX1"]}
FLOOR = 50
N_PERM_SMALL = 20000
N_PERM_LARGE = 4000
SEED = 20260926

# R9: the frozen low end.  Rank cut from the low-end scan (R10 there), plus the
# six compartments H2 went looking for.  Names that do not resolve against the
# reference table are reported rather than silently dropped.
LOW_END_EXTRA = ("urothelial cells", "retinal pigment epithelial cells",
                 "rod photoreceptor cells", "cone photoreceptor cells",
                 "sertoli cells", "choroid plexus epithelial cells")
LOW_END_RANK = 116


def delta_from_cpm(cpm):
    comps = [sum(math.log10(cpm.get(g, 0.0) + 1.0) for g in genes) / len(genes)
             for genes in MODULE.values()]
    return sum(comps) / len(comps)


def load_pseudobulk(path):
    rows = {}
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        rows[(r["atlas"], r["cl_id"])] = {
            "n_cells": int(r["n_cells"]),
            "lib_sum": float(r["lib_sum"]),
            "delta": float(r["delta"]),
            **{f"cpm_{g}": float(r[f"cpm_{g}"]) for g in GENES},
        }
    return rows


def load_universe(path):
    """(atlas, cl_id) -> n_cells.  Rows without a CL id are skipped."""
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        cl = (r.get("cl_id") or "").strip()
        if not cl:
            continue
        out[(r["atlas"], cl)] = int(r["n_cells"])
    return out


def load_reference():
    ref, rank = {}, {}
    for r in csv.DictReader(open(REF, encoding="utf-8"), delimiter="\t"):
        ref[r["cell_type"]] = float(r["delta_restriction"])
        rank[r["cell_type"]] = int(r["rank"])
    return ref, rank


def load_mapping():
    rows = list(csv.DictReader(open(MAP, encoding="utf-8"), delimiter="\t"))
    return rows, [r for r in rows if r["testable"] == "yes"]


def atlas_stats(pb, universe, floor):
    """Per-atlas (mean, sd, n_labels) of Delta over the atlas's labels >= floor.

    An atlas whose whole label set is one term has sd = 0; the divisor falls
    back to 1.0 exactly as h2b did, which makes every pair from that atlas land
    on z = 0.  That is a property of the atlas having no contrast, not a
    measurement, so such pairs are flagged `constructive_zero` downstream.
    """
    stats = {}
    for atlas in sorted({a for a, _ in pb}):
        labs = [c for (a, c), n in universe.items()
                if a == atlas and n >= floor and (a, c) in pb]
        vals = [pb[(atlas, c)]["delta"] for c in labs]
        if not vals:
            continue
        mu = sum(vals) / len(vals)
        sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals)) or 1.0
        stats[atlas] = (mu, sd, len(vals))
    return stats


def build_pairs(onto, pb, mapping_rows, universe, stats, floor, atlases=None):
    terms = defaultdict(list)
    for (a, c), n in universe.items():
        if n >= floor and (a, c) in pb and (atlases is None or a in atlases):
            terms[a].append(c)
    pairs = []
    for m in mapping_rows:
        cl, hpa = m["cl_id"], m["hpa_cell_type"]
        for atlas in sorted(terms):
            matched = [t for t in terms[atlas] if cl in onto.ancestors(t)]
            if not matched:
                continue
            members = [pb[(atlas, t)] for t in matched]
            lib = sum(x["lib_sum"] for x in members)
            if lib <= 0:
                continue
            cpm = {g: sum(x[f"cpm_{g}"] * x["lib_sum"] for x in members) / lib
                   for g in GENES}
            mu, sd, n_lab = stats[atlas]
            delta = delta_from_cpm(cpm)
            pairs.append({
                "convention": None,
                "hpa_cell_type": hpa,
                "cl_id": cl,
                "atlas": atlas,
                "constructive_zero": n_lab < 2,
                "atlas_terms": ";".join(sorted(matched)),
                "n_atlas_terms": len(matched),
                "n_cells": sum(x["n_cells"] for x in members),
                "reference_delta": None,
                "atlas_delta": delta,
                "atlas_delta_z": (delta - mu) / sd,
            })
    return pairs


def pooled(pairs, xmap, draws=None):
    """Pooled Spearman for a pair set under a given reference-Delta vector.

    `xmap` is the reference value per cell type.  It is a genuine degree of
    freedom: any strictly monotone transform of x leaves a *pooled* rho
    unchanged, but two vectors that are only 0.993 rank-correlated are not
    monotone-equivalent and do move it (see F1 in ldc_02_diagnose.py).
    """
    xs = [xmap[p["hpa_cell_type"]] for p in pairs]
    ys = [p["atlas_delta_z"] for p in pairs]
    n = len(xs)
    if n < 3:
        return {"n": n, "rho": None, "p_perm": None}
    rho = spearman(xs, ys)
    out = {
        "n": n,
        "rho": round(rho, 4),
        "p_t_approx": round(spearman_p_value(rho, n), 4),
        "min_detectable_rho_power80": round(min_detectable_rho(n), 3),
        "achieved_power_at_observed": round(achieved_power(rho, n), 3),
    }
    nd = draws if draws is not None else (N_PERM_SMALL if n <= 60 else N_PERM_LARGE)
    if nd:
        rng, y, ge = random.Random(SEED), list(ys), 0
        for _ in range(nd):
            rng.shuffle(y)
            if abs(spearman(xs, y)) >= abs(rho) - 1e-12:
                ge += 1
        out["p_perm"] = round((1 + ge) / (1 + nd), 4)
        out["p_perm_n"] = nd
    else:
        out["p_perm"] = None
        out["p_perm_n"] = 0
    return out


def main():
    onto = Ontology(OBO)
    ref, rank = load_reference()
    mapping, testable = load_mapping()
    # F1 of ldc_02_diagnose.py: b5b did not use the tabulated Delta.  It rebuilt
    # Delta from the HPA wide table.  The two vectors are 0.9934 rank-correlated
    # but not monotone-equivalent, so the pooled rho depends on the choice; both
    # are carried through every set below rather than one being picked.
    rebuilt = {ct: hpa_delta(row) for ct, row in load_hpa().items()}
    xmaps = {"tabulated_delta_restriction": ref, "rebuilt_from_wide": rebuilt}
    x_note = {
        "rank_correlation_between_the_two": round(
            spearman([ref[c] for c in sorted(set(ref) & set(rebuilt))],
                     [rebuilt[c] for c in sorted(set(ref) & set(rebuilt))]), 4),
        "b5b_used": "rebuilt_from_wide",
        "manuscript_delta_is": "tabulated_delta_restriction",
    }
    pb19, pb07 = load_pseudobulk(P19), load_pseudobulk(P07)
    uni19, uni07 = load_universe(INV19), load_universe(LBL07)
    pb = {**pb19, **pb07}
    universe = {**uni19, **uni07}
    stats = atlas_stats(pb, universe, FLOOR)

    missing = [t for t in LOW_END_EXTRA if t not in ref]
    low_end = {t for t, r in rank.items() if r >= LOW_END_RANK}
    low_end |= {t for t in LOW_END_EXTRA if t in ref}

    summary = {
        "floor": FLOOR,
        "reference_conventions": x_note,
        "n_reference_types": len(ref),
        "n_low_end_members": len(low_end),
        "low_end_names_unresolved": missing,
        "atlases": {"total": len({a for a, _ in pb}),
                    "cached_19": len({a for a, _ in pb19}),
                    "new_7": len({a for a, _ in pb07})},
        "atlas_label_universe": {a: v[2] for a, v in sorted(stats.items())},
    }

    out_pairs, agreement = [], []
    for conv, rows in (("A", testable), ("B", mapping)):
        sets = {
            "19_atlases": {a for a, _ in pb19},
            "07_new_atlases": {a for a, _ in pb07},
            "26_atlases": set(stats),
        }
        built = {}
        for name, atlases in sets.items():
            pairs = build_pairs(onto, pb, rows, universe, stats, FLOOR, atlases)
            for p in pairs:
                p["convention"], p["reference_delta"] = conv, ref[p["hpa_cell_type"]]
            built[name] = pairs
            for xname, xmap in xmaps.items():
                agreement.append({"convention": conv, "set": name,
                                  "x_convention": xname, **pooled(pairs, xmap)})

        pairs = built["26_atlases"]
        for p in pairs:
            p["low_end"] = p["hpa_cell_type"] in low_end
        lo = [p for p in pairs if p["low_end"] and not p["constructive_zero"]]
        cz = [p for p in pairs if p["low_end"] and p["constructive_zero"]]
        hi = [p for p in pairs if not p["low_end"]]
        for xname, xmap in xmaps.items():
            agreement.append({"convention": conv, "set": "26_low_end",
                              "x_convention": xname, **pooled(lo, xmap)})
            agreement.append({"convention": conv, "set": "26_low_end_with_cz",
                              "x_convention": xname,
                              **pooled(lo + cz, xmap)})
            agreement.append({"convention": conv, "set": "26_rest",
                              "x_convention": xname, **pooled(hi, xmap)})

        # leave one atlas out, over the combined set
        loo = []
        for atlas in sorted({p["atlas"] for p in pairs}):
            sub = [p for p in pairs if p["atlas"] != atlas]
            r = pooled(sub, ref, draws=0)
            loo.append({"atlas_dropped": atlas, "n": r["n"], "rho": r["rho"]})
        summary[f"leave_one_atlas_out_{conv}"] = loo

        # per atlas
        per_atlas = []
        for atlas in sorted({p["atlas"] for p in pairs}):
            sub = [p for p in pairs if p["atlas"] == atlas]
            xs = [ref[p["hpa_cell_type"]] for p in sub]
            ys = [p["atlas_delta_z"] for p in sub]
            rho = spearman(xs, ys) if len(sub) >= 2 else float("nan")
            per_atlas.append({"atlas": atlas, "n": len(sub),
                              "rho": None if math.isnan(rho) else round(rho, 4)})
        summary[f"per_atlas_{conv}"] = per_atlas
        out_pairs.extend(pairs)

    # ---- self-checks (R8) -------------------------------------------------
    checks = {}
    got = {(a["convention"], a["set"], a["x_convention"]): a
           for a in agreement}
    primary = {(c, s): got[(c, s, "tabulated_delta_restriction")]
               for c in "AB" for s in
               ("19_atlases", "07_new_atlases", "26_atlases", "26_low_end",
                "26_low_end_with_cz", "26_rest")}
    checks["S1_b5b_19_atlases_A"] = {
        "expected": {"n": 162, "rho": 0.5191},
        "note": "b5b used the wide-rebuilt Delta as x, this run's primary x is "
                "the tabulated one, so the rho is expected to differ",
        "observed_tabulated": {k: primary[("A", "19_atlases")][k]
                               for k in ("n", "rho")},
        "observed_wide_rebuilt": {k: got[("A", "19_atlases", "rebuilt_from_wide")][k]
                                  for k in ("n", "rho")},
    }
    checks["S2_h2b_7_atlases_B"] = {
        "expected": {"n": 38, "rho": 0.568},
        "note": "h2b pooled each pair as the unweighted mean of its member "
                "labels' z values and dropped atlases with fewer than three "
                "matched reference types, so both n and rho are expected to "
                "differ under the frozen b5b pooling rule used here",
        "observed": {k: primary[("B", "07_new_atlases")][k]
                     for k in ("n", "rho")},
    }
    for c in checks.values():
        if "observed_wide_rebuilt" in c:
            c["pass_wide_rebuilt"] = (
                c["expected"]["n"] == c["observed_wide_rebuilt"]["n"]
                and abs(c["expected"]["rho"]
                        - c["observed_wide_rebuilt"]["rho"]) <= 1e-4)
        elif "observed" in c:
            c["pass"] = False  # see the note: the gap is explained, not hidden

    def compare_fields(pairs, path, ycol, key_col, keys):
        frozen = {(r["atlas"], r[key_col]): r
                  for r in csv.DictReader(open(path, encoding="utf-8"),
                                          delimiter="\t")}
        worst, shared = 0.0, 0
        for p in pairs:
            k = (p["atlas"], p["hpa_cell_type"])
            if k not in frozen:
                continue
            shared += 1
            worst = max(worst, abs(p["atlas_delta_z"] - float(frozen[k][ycol])))
        return {"shared_rows": shared, "max_abs_diff": round(worst, 6),
                "frozen_rows": len(frozen), "key_fields": keys}

    pairs_A_19 = build_pairs(onto, pb, testable, universe, stats, FLOOR,
                             {a for a, _ in pb19})
    pairs_B_07 = build_pairs(onto, pb, mapping, universe, stats, FLOOR,
                             {a for a, _ in pb07})
    checks["S1_fieldwise_vs_b5b"] = compare_fields(
        pairs_A_19, B5B_PAIRS, "y_all_labels", "hpa_cell_type",
        "atlas + hpa_cell_type")
    checks["S2_fieldwise_vs_h2b"] = compare_fields(
        pairs_B_07, H2B_PAIRS, "atlas_delta_z", "reference_cell_type",
        "atlas + reference_cell_type")
    for c in checks.values():
        exp, obs = c.get("expected"), c.get("observed")
        if exp and obs:
            c["pass"] = (exp["n"] == obs["n"]
                         and abs(exp["rho"] - (obs["rho"] or 0)) <= 1e-4)

    summary["self_checks"] = checks
    summary["headline"] = {
        "A_19": primary[("A", "19_atlases")],
        "A_26": primary[("A", "26_atlases")],
        "B_07": primary[("B", "07_new_atlases")],
        "B_26": primary[("B", "26_atlases")],
    }

    with open(os.path.join(HERE, "ldc_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(HERE, "ldc_agreement.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        cols = ["convention", "set", "x_convention", "n", "rho", "p_perm", "p_perm_n",
                "p_t_approx", "min_detectable_rho_power80",
                "achieved_power_at_observed"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(agreement)
    with open(os.path.join(HERE, "ldc_pairs.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        cols = ["convention", "hpa_cell_type", "cl_id", "atlas", "low_end",
                "constructive_zero",
                "atlas_terms", "n_atlas_terms", "n_cells", "reference_delta",
                "atlas_delta", "atlas_delta_z"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for p in sorted(out_pairs, key=lambda p: (p["convention"],
                                                  p["reference_delta"])):
            w.writerow({**p, "atlas_delta": round(p["atlas_delta"], 6),
                        "atlas_delta_z": round(p["atlas_delta_z"], 6)})

    print(json.dumps({"headline": summary["headline"],
                      "self_checks": checks,
                      "low_end_names_unresolved": missing,
                      "n_low_end_members": len(low_end)},
                     ensure_ascii=False, indent=1))
    print("\nagreement:")
    for a in agreement:
        print(f"  {a['convention']} {a['set']:20s} x={a['x_convention'][:12]:12s} "
              f"n={a['n']:>4} rho={a['rho']} p={a.get('p_perm')}")
    print("\nlow-end pairs (convention B):")
    for p in sorted([p for p in out_pairs
                     if p["convention"] == "B" and p["low_end"]],
                    key=lambda p: p["reference_delta"]):
        print(f"  {p['hpa_cell_type']:36s} ref {p['reference_delta']:+7.3f} "
              f"atlas {p['atlas_delta_z']:+7.3f} n={p['n_cells']:>7d} "
              f"[{p['atlas'][:28]}]")


if __name__ == "__main__":
    main()
