"""X5 - from-scratch recomputation of the Fig. 1a denominator chain.

Input (read-only):
  analysis/vp24_pan_species_20260925/coordinator_crosscheck/
      vp24_six_species_positions.tsv

Definitions (frozen in _coord/_round_L_20261002/out/W6/w6_denominators.tsv):
  retrieved            every row of the positions table
  excluded_fragment    a row that does not span the four-position window
                       (its called residues at 83/135/140/141 are all absent)
  assessed             retrieved - excluded_fragment
  uncalled_at_some_position
                       an assessed row in which at least one of the four
                       positions is X (unresolved)
  called_all_four      assessed - uncalled_at_some_position
  non-Bundibugyo       called_all_four summed over the five non-BDBV species
"""

from __future__ import annotations

import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"G:\本迪布焦研究"
SRC = os.path.join(ROOT, "analysis", "vp24_pan_species_20260925",
                   "coordinator_crosscheck", "vp24_six_species_positions.tsv")
FOUR = ["pos83", "pos135", "pos140", "pos141"]
ORDER = ["BDBV", "EBOV", "SUDV", "TAFV", "RESTV", "BOMV"]

d = pd.read_csv(SRC, sep="\t")
for c in FOUR:
    d[c] = d[c].fillna("-").astype(str).str.strip().str.upper()
d["length_aa"] = pd.to_numeric(d["length_aa"], errors="coerce")

rows = []
for sp in ORDER:
    s = d[d["species"] == sp]
    retrieved = len(s)
    # a fragment is a row whose four window positions are all uncalled AND
    # which is shorter than the window (does not span it); a row that spans the
    # window but has an unresolved call carries X, not an absent residue.
    excl = s[(s[FOUR] == "-").all(axis=1) & (s["length_aa"] < 200)]
    assessed = s.drop(excl.index)
    uncalled = assessed[(assessed[FOUR] == "X").any(axis=1)]
    called = assessed[~assessed.index.isin(uncalled.index)]
    rows.append({
        "species": sp, "retrieved": retrieved,
        "excluded_fragment": len(excl), "assessed": len(assessed),
        "uncalled_at_some_position": len(uncalled),
        "called_all_four": len(called),
        "note": ""})

out = pd.DataFrame(rows, columns=["species", "retrieved", "excluded_fragment",
                                  "assessed", "uncalled_at_some_position",
                                  "called_all_four", "note"])
tot = {c: int(out[c].sum()) for c in
       ["retrieved", "excluded_fragment", "assessed",
        "uncalled_at_some_position", "called_all_four"]}
tot["species"] = "TOTAL"
tot["note"] = (f"non-Bundibugyo called_all_four = "
               f"{tot['called_all_four'] - int(out.loc[out['species']=='BDBV','called_all_four'].iloc[0])}")
for c in ["retrieved", "excluded_fragment", "assessed",
          "uncalled_at_some_position", "called_all_four"]:
    tot[c] = tot[c]
out = pd.concat([out, pd.DataFrame([tot])[out.columns]], ignore_index=True)

out.to_csv(os.path.join(HERE, "X5_fig1a_denominator_recompute.tsv"), sep="\t",
           index=False)
print(out.to_string(index=False))
print()
print("caption chain: 3522 / 32 / 3490 / 19 / 3471 / 3440")
nonbdbv = tot["called_all_four"] - int(
    out.loc[out["species"] == "BDBV", "called_all_four"].iloc[0])
print("recomputed   :", tot["retrieved"], "/", tot["excluded_fragment"], "/",
      tot["assessed"], "/", tot["uncalled_at_some_position"], "/",
      tot["called_all_four"], "/", nonbdbv)
