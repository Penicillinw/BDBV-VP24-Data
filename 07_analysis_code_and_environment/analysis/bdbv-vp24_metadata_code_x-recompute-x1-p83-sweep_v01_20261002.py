"""X1 - exhaustive sweep of every natural definition of "P83 to the second
VP24 copy" in the Fig. 1b caption.

Read-only over the workspace. Writes X1_p83_distance_sweep.tsv next to it.

Structures
----------
4U2X.pdb            experimental EBOV VP24 / KPNA5 complex; chains A, B, C are
                    EBOV VP24 (residue 83 = PRO) and D, E, F are human KPNA5.
vp24_kpna_complex.pdb   the shipped Fig. 1b homology model: chain A and chain B
                    are the BDBV VP24 model (residue numbering already on EBOV
                    full-length numbering, 11-231; position 83 = SER) and chain Z
                    is the KPNA5 chain (335-506).
bdbv_vp24_kpna5_from4u2x.pdb       template-transfer BDBV VP24 on KPNA5
bdbv_vp24_kpna1_from4u2x_paralog.pdb  template-transfer BDBV VP24 on KPNA1
                    (both renumbered 1-216 on the construct; their position 83 in
                    *EBOV* numbering is construct residue 68 - the model offset is
                    checked in-script, not assumed).

Metrics reported for every pair that contains the P83-bearing VP24 chain:
  min_heavy_atom          minimum distance between any P83 heavy atom and any
                          partner-chain heavy atom
  ca_ca                   P83 CA to the partner's *same-numbered* residue CA
  ca_to_partner_centroid  P83 CA to the partner-chain heavy-atom centroid
  residue_centroid_pair   P83 heavy-atom centroid to the partner's same-numbered
                          residue heavy-atom centroid
  chain_centroid_pair     whole-chain heavy-atom centroid to whole-chain
                          heavy-atom centroid
  bb_N / bb_C / bb_O / bb_CB   corresponding backbone/side-chain atom pairs
  p83_to_nearest_B_atom   P83 to the nearest atom of the literal chain labelled B
"""

from __future__ import annotations

import itertools
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"G:\本迪布焦研究"

PDBS = {
    "4U2X": os.path.join(ROOT, "analysis", "t4_vp24_struct_seqmap_20260925",
                         "4U2X.pdb"),
    "model_vp24_kpna_complex": os.path.join(
        ROOT, "figures", "pymol_publication", "complexes",
        "vp24_kpna_complex.pdb"),
    "tmpl_kpna5": os.path.join(ROOT, "data", "template_transfer_20260919",
                               "bdbv_vp24_kpna5_from4u2x.pdb"),
    "tmpl_kpna1": os.path.join(ROOT, "data", "template_transfer_20260919",
                               "bdbv_vp24_kpna1_from4u2x_paralog.pdb"),
}


def parse_pdb(path):
    """Return {chain: {'xyz': (n,3) float, 'elem': [..], 'atom': [..],
    'resseq': np.array int, 'resname': [..]}} keeping heavy atoms only."""
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom = line[12:16].strip()
            element = line[76:78].strip().upper()
            if not element:
                element = atom[0].upper()
            if element == "H" or atom.startswith("H"):
                continue
            rows.append((line[21], int(line[22:26]), line[26].strip(),
                         line[17:20].strip(), atom,
                         np.array([float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])])))
    chains = {}
    for ch, resseq, icode, resname, atom, xyz in rows:
        d = chains.setdefault(ch, {"xyz": [], "elem": [], "atom": [],
                                   "resseq": [], "resname": [], "icode": []})
        d["xyz"].append(xyz)
        d["elem"].append(element)
        d["atom"].append(atom)
        d["resseq"].append(resseq)
        d["resname"].append(resname)
        d["icode"].append(icode)
    for ch, d in chains.items():
        d["xyz"] = np.asarray(d["xyz"])
        d["resseq"] = np.asarray(d["resseq"])
    return chains


CHAIN = {}
for name, path in PDBS.items():
    CHAIN[name] = parse_pdb(path)


def res_mask(chains, ch, resseq):
    c = chains[ch]
    return (c["resseq"] == resseq) & np.array(
        [i == "" or i == " " for i in c["icode"]])


def resname_of(chains, ch, resseq):
    m = res_mask(chains, ch, resseq)
    return chains[ch]["resname"][int(np.argmax(m))] if m.any() else None


def centroid(xyz):
    return xyz.mean(axis=0)


def mindist(a, b):
    # (n,3) x (m,3)
    d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
    return float(d.min())


def atom_xyz(chains, ch, resseq, atom):
    c = chains[ch]
    m = res_mask(chains, ch, resseq) & np.array(
        [a == atom for a in c["atom"]])
    if not m.any():
        return None
    return c["xyz"][int(np.argmax(m))]


ROWS = []


def add(structure, pair, metric, value, note):
    ROWS.append({"structure": structure, "pair": pair, "metric": metric,
                 "value_A": "" if value is None else round(float(value), 3),
                 "note": note})


def pair_metrics(structure, chains, vp24_chain, partner_chain, p83_resseq,
                 partner_resseq83):
    """All metrics for one P83(chain vp24) -> partner chain pair."""
    c = chains[vp24_chain]
    p = chains[partner_chain]
    vp24_xyz = c["xyz"]
    p83_xyz = vp24_xyz[res_mask(chains, vp24_chain, p83_resseq)]
    p83_ca = atom_xyz(chains, vp24_chain, p83_resseq, "CA")
    pcent = centroid(p83_xyz)
    ccent = centroid(vp24_xyz)
    part_cent = centroid(p["xyz"])
    tag = f"{vp24_chain}-{partner_chain}"

    add(structure, tag, "min_heavy_atom", mindist(p83_xyz, p["xyz"]),
        "P83 any heavy atom to partner any heavy atom")

    if p83_ca is not None:
        add(structure, tag, "ca_to_partner_centroid",
            np.linalg.norm(p83_ca - part_cent),
            "P83 CA to partner heavy-atom centroid")
        add(structure, tag, "ca_to_nearest_partner_atom",
            float(np.sqrt(((p["xyz"] - p83_ca) ** 2).sum(-1)).min()),
            "P83 CA to nearest partner heavy atom")

    if partner_resseq83 is not None:
        p83p = p["xyz"][res_mask(chains, partner_chain, partner_resseq83)]
        pcent_p = centroid(p83p)
        add(structure, tag, "residue_centroid_pair",
            np.linalg.norm(pcent - pcent_p),
            f"P83 centroid to {partner_chain}{partner_resseq83} centroid")
        ca_p = atom_xyz(chains, partner_chain, partner_resseq83, "CA")
        if p83_ca is not None and ca_p is not None:
            add(structure, tag, "ca_ca", np.linalg.norm(p83_ca - ca_p),
                f"P83 CA to {partner_chain}{partner_resseq83} CA")
        for atom in ("N", "C", "O", "CB"):
            a1 = atom_xyz(chains, vp24_chain, p83_resseq, atom)
            a2 = atom_xyz(chains, partner_chain, partner_resseq83, atom)
            if a1 is not None and a2 is not None:
                add(structure, tag, f"bb_{atom}", np.linalg.norm(a1 - a2),
                    f"P83 {atom} to {partner_chain}{partner_resseq83} {atom}")

    add(structure, tag, "chain_centroid_pair",
        np.linalg.norm(ccent - part_cent),
        f"centroid(chain {vp24_chain}) to centroid(chain {partner_chain})")
    if p83_ca is not None:
        add(structure, tag, "p83ca_to_own_chain_centroid",
            np.linalg.norm(p83_ca - ccent), "P83 CA to own-chain centroid")


# ---------------------------------------------------------------- 4U2X
u = CHAIN["4U2X"]
for b1, b2 in itertools.combinations(("A", "B", "C"), 2):
    pair_metrics("4U2X(EBOV VP24)", u, b1, b2, 83, 83)
# P83 -> nearest KPNA5 atom (cognate and non-cognate)
for vp24c, kpnac in (("A", "D"), ("A", "E"), ("A", "F")):
    add("4U2X(EBOV VP24)", f"{vp24c}-{kpnac}", "min_heavy_atom",
        mindist(u[vp24c]["xyz"][res_mask(u, vp24c, 83)], u[kpnac]["xyz"]),
        "P83 to nearest KPNA5 heavy atom")
# every VP24 copy's own P83 -> nearest KPNA5 (the diagnostic that explains 55.0)
for vp24c, kpnac in (("B", "E"), ("C", "F"), ("B", "D"), ("C", "D"), ("C", "E")):
    add("4U2X(EBOV VP24)", f"{vp24c}-{kpnac}-per-copy-p83",
        "copy_p83_to_nearest_partner_atom",
        mindist(u[vp24c]["xyz"][res_mask(u, vp24c, 83)], u[kpnac]["xyz"]),
        f"P83 of copy {vp24c} to nearest atom of KPNA5 chain {kpnac}")

# literal "nearest chain B atom" from each VP24 chain
for src in ("A", "B", "C"):
    add("4U2X(EBOV VP24)", f"{src}-B", "p83_to_nearest_B_atom",
        mindist(u[src]["xyz"][res_mask(u, src, 83)], u["B"]["xyz"]),
        "P83 to nearest atom of literal chain B")

# ------------------------------------------------- shipped homology model
m = CHAIN["model_vp24_kpna_complex"]
pair_metrics("model(vp24_kpna_complex)", m, "A", "B", 83, 83)
add("model(vp24_kpna_complex)", "A-Z", "min_heavy_atom",
    mindist(m["A"]["xyz"][res_mask(m, "A", 83)], m["Z"]["xyz"]),
    "P83 to nearest KPNA5(Z) heavy atom")
add("model(vp24_kpna_complex)", "B-Z-per-copy-p83",
    "copy_p83_to_nearest_partner_atom",
    mindist(m["B"]["xyz"][res_mask(m, "B", 83)], m["Z"]["xyz"]),
    "P83 of the SECOND VP24 copy (chain B) to nearest KPNA5(Z) atom "
    "- the traceable source of the caption's '55.0'")
add("model(vp24_kpna_complex)", "A-B", "p83_to_nearest_B_atom",
    mindist(m["A"]["xyz"][res_mask(m, "A", 83)], m["B"]["xyz"]),
    "P83 to nearest atom of literal chain B (second VP24 copy)")

# ------------------------------------------- template-transfer BDBV models
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def seqdict(chains, ch):
    return {int(r): AA3.get(n, "X")
            for r, n in zip(chains[ch]["resseq"], chains[ch]["resname"])}


def best_offset(template_seq, ref_seq):
    """Offset o such that construct c corresponds to EBOV residue c + o."""
    best, bestn = None, -1
    for o in range(-40, 60):
        n = sum(1 for c, a in template_seq.items()
                if ref_seq.get(c + o) == a)
        if n > bestn:
            best, bestn = o, n
    return best, bestn, len(template_seq)


REF = seqdict(CHAIN["4U2X"], "A")     # EBOV-numbered reference sequence
for name, tag in (("tmpl_kpna5", "bdbv_vp24_kpna5_from4u2x"),
                  ("tmpl_kpna1", "bdbv_vp24_kpna1_from4u2x_paralog")):
    t = CHAIN[name]
    off, nmatch, ntot = best_offset(seqdict(t, "V"), REF)
    c83 = 83 - off                       # construct resseq == EBOV 83
    add(name, "V", "numbering_offset_check", None,
        f"construct = EBOV numbering minus offset {off} "
        f"({nmatch}/{ntot} residues match 4U2X chain A); EBOV 83 -> "
        f"construct {c83} = {resname_of(t, 'V', c83)}")
    p83_xyz = t["V"]["xyz"][res_mask(t, "V", c83)]
    add(name, f"V-p83(c{c83})", "min_heavy_atom_to_partner",
        mindist(p83_xyz, t["H"]["xyz"]),
        f"BDBV position 83 (construct residue {c83}, "
        f"{resname_of(t, 'V', c83)}) to nearest partner heavy atom")
    add(name, f"V-p83(c{c83})", "p83_to_partner_centroid",
        np.linalg.norm(centroid(p83_xyz) - centroid(t["H"]["xyz"])),
        "BDBV 83 centroid to partner heavy-atom centroid")
    # the very same file also carries a KPNA1/KPNA5 chain H; no second VP24 copy
    add(name, "V", "second_vp24_copy", None,
        "file contains one VP24 chain (V) only - no second copy to measure")

df = pd.DataFrame(ROWS, columns=["structure", "pair", "metric", "value_A",
                                 "note"])
out = os.path.join(HERE, "X1_p83_distance_sweep.tsv")
df.to_csv(out, sep="\t", index=False)

print("wrote", out, len(df), "rows")
print()
cand = df[pd.to_numeric(df["value_A"], errors="coerce").notna()].copy()
cand["value_A"] = cand["value_A"].astype(float)
hit = cand[(cand["value_A"] >= 54.9) & (cand["value_A"] <= 55.1)]
print("rows with value in [54.9, 55.1]:")
print(hit.to_string(index=False) if len(hit) else "  NONE")
print()
print("closest 12 rows to 55.0 A:")
near = cand.assign(d=(cand["value_A"] - 55.0).abs()).sort_values("d").head(12)
print(near[["structure", "pair", "metric", "value_A", "note"]].to_string(
    index=False))
