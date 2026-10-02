"""X2 - authoritative recount of the Fig. 1b inter-chain polar contacts.

Recomputed from scratch from the raw coordinates; nothing is read from any
existing contact table.

  analysis/t4_vp24_struct_seqmap_20260925/4U2X.pdb
      chains A, B, C = EBOV VP24, chains D, E, F = human KPNA5 (cognate pairs
      A-D, B-E, C-F).  H atoms absent (X-ray, 3.15 A).
  figures/pymol_publication/complexes/vp24_kpna_complex.pdb
      chain A = the BDBV VP24 homology model (EBOV full-length numbering),
      chain B = the second copy of that model, chain Z = KPNA5 (335-506).

For every cross-chain heavy-atom pair within 4.0 A the script records the
distance and whether both atoms are polar (element N or O).  The caption's
criterion is "heavy atom, inter-chain, polar, < 3.6 A".
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"G:\本迪布焦研究"
PDB_4U2X = os.path.join(ROOT, "analysis", "t4_vp24_struct_seqmap_20260925",
                        "4U2X.pdb")
PDB_MODEL = os.path.join(ROOT, "figures", "pymol_publication", "complexes",
                         "vp24_kpna_complex.pdb")
COGNATE = [("A", "D"), ("B", "E"), ("C", "F")]
FOUR = (135, 140, 141, 184)


def load(path, keep_h=False):
    out = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom = line[12:16].strip()
        el = (line[76:78].strip().upper() or atom[0].upper())
        if not keep_h and (el == "H" or atom.startswith("H")):
            continue
        out.setdefault(line[21], []).append({
            "resseq": int(line[22:26]), "resname": line[17:20].strip(),
            "atom": atom, "elem": el,
            "xyz": np.array([float(line[30:38]), float(line[38:46]),
                             float(line[46:54])])})
    return out


def cross_pairs(chains, vp24_chain, partner_chain, cutoff):
    a, b = chains[vp24_chain], chains[partner_chain]
    rows = []
    for x in a:
        for y in b:
            d = float(np.linalg.norm(x["xyz"] - y["xyz"]))
            if d < cutoff:
                rows.append({
                    "chain_pair": f"{vp24_chain}-{partner_chain}",
                    "viral_residue": x["resseq"], "viral_resname": x["resname"],
                    "viral_atom": x["atom"], "partner_residue": y["resseq"],
                    "partner_resname": y["resname"], "partner_atom": y["atom"],
                    "distance_A": round(d, 3),
                    "is_polar": (x["elem"] in ("N", "O")
                                 and y["elem"] in ("N", "O")),
                    "cutoff_3.6": d < 3.6, "cutoff_4.0": True})
    return rows


rows = []
u = load(PDB_4U2X)
for v, k in COGNATE:
    rows += cross_pairs(u, v, k, 4.0)
m = load(PDB_MODEL)
rows += cross_pairs(m, "A", "Z", 4.0)

df = pd.DataFrame(rows)[["chain_pair", "viral_residue", "viral_resname",
                         "viral_atom", "partner_residue", "partner_resname",
                         "partner_atom", "distance_A", "is_polar", "cutoff_3.6",
                         "cutoff_4.0"]]
df.to_csv(os.path.join(HERE, "X2_contacts_authoritative.tsv"), sep="\t",
          index=False)

p = df[df["is_polar"]]
print("=== all cross-chain heavy-atom pairs < 4.0 A ===")
print(df.groupby("chain_pair").size().to_string())
print()
print("=== polar (N/O ... N/O) counts ===")
for cp in df["chain_pair"].unique():
    sub = p[p["chain_pair"] == cp]
    print(f"  {cp}: <4.0 A = {len(sub):3d}   <3.6 A = "
          f"{int(sub['cutoff_3.6'].sum()):3d}")
print()
print("=== 4U2X A-D polar pairs < 4.0 A, sorted ===")
ad = p[(p["chain_pair"] == "A-D")].sort_values("distance_A")
print(ad[["viral_residue", "viral_resname", "viral_atom", "partner_residue",
          "partner_resname", "partner_atom", "distance_A", "cutoff_3.6"]]
      .to_string(index=False))
print()
print("=== 4U2X A-D polar pairs < 3.6 A (the caption criterion) ===")
ad36 = ad[ad["cutoff_3.6"]]
print(f"  count = {len(ad36)}")
for _, r in ad36.iterrows():
    print(f"  {r.viral_resname}{r.viral_residue} {r.viral_atom} -> "
          f"{r.partner_resname}{r.partner_residue} {r.partner_atom}  "
          f"{r.distance_A:.3f}")
print()
print("=== the contested R140 NH2 ... D475 O pair (4U2X) ===")
tgt = df[(df["viral_residue"] == 140) & (df["viral_atom"].isin(["NH1", "NH2"]))
         & (df["partner_residue"] == 475)]
print(tgt.to_string(index=False) if len(tgt) else "  NOT FOUND within 4.0 A")

# exhaustive check restricted to the four variant positions
print()
print("=== 4U2X, polar pairs < 4.0 A limited to VP24 positions 135/140/141/184 ===")
fourv = p[(p["viral_residue"].isin(FOUR)) & (p["chain_pair"] == "A-D")]
print(fourv[["viral_residue", "viral_atom", "partner_residue",
             "partner_atom", "distance_A", "cutoff_3.6"]].sort_values(
    "distance_A").to_string(index=False))
print()
print("=== BDBV model chain A-Z polar pairs < 4.0 A limited to 135/140/141/184 ===")
mf = p[(p["chain_pair"] == "A-Z") & (p["viral_residue"].isin(FOUR))]
print(mf[["viral_residue", "viral_atom", "partner_residue", "partner_atom",
          "distance_A", "cutoff_3.6"]].sort_values("distance_A").to_string(
    index=False))
print()
print("=== BDBV model chain A-Z: ALL polar pairs < 3.6 A ===")
allm = p[(p["chain_pair"] == "A-Z") & (p["cutoff_3.6"])].sort_values(
    "distance_A")
print(f"  count = {len(allm)}")
print(allm[["viral_residue", "viral_atom", "partner_residue", "partner_atom",
            "distance_A"]].to_string(index=False))
