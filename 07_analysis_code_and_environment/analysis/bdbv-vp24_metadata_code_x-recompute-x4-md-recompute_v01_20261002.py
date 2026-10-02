"""X4 - inventory and from-scratch recomputation of the pilot MD workspace.

Inputs (read-only):
  analysis/t4_vp24_openmm_20260925/<system>/{A,B,D,mutant,minimised}.pdb
  analysis/t4_vp24_openmm_20260925/openmm_results.json
  analysis/t4_p83_stability_20260925/{vp24_WT,vp24_AP83}.pdb
  analysis/t4_p83_stability_20260925/{pilot_run,pilot_bench}.json

Numbering: the native files (mutant.pdb, WT.pdb) carry 4U2X/EBOV full-length
numbering (VP24 16-231, KPNA5 335-506).  The prepared files (A.pdb, B.pdb,
minimised.pdb) are renumbered on the construct, where construct c corresponds
to EBOV residue c + 15 (checked in-script, not assumed).
"""

from __future__ import annotations

import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"G:\本迪布焦研究"
OPENMM = os.path.join(ROOT, "analysis", "t4_vp24_openmm_20260925")
P83 = os.path.join(ROOT, "analysis", "t4_p83_stability_20260925")
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def parse(path, heavy_only=False):
    """{chain: {'xyz','elem','restype','resseq'}}"""
    ch = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        at = line[12:16].strip()
        el = (line[76:78].strip().upper() or at[0].upper())
        if heavy_only and (el == "H" or at.startswith("H")):
            continue
        d = ch.setdefault(line[21], {"xyz": [], "elem": [], "restype": [],
                                     "resseq": [], "atom": []})
        d["xyz"].append([float(line[30:38]), float(line[38:46]),
                         float(line[46:54])])
        d["elem"].append(el)
        d["restype"].append(AA3.get(line[17:20].strip(), "X"))
        d["resseq"].append(int(line[22:26]))
        d["atom"].append(at)
    for d in ch.values():
        d["xyz"] = np.asarray(d["xyz"], dtype=float)
    return ch


def res_ca(chains, chain, resseq):
    d = chains[chain]
    idx = [i for i, (r, a) in enumerate(zip(d["resseq"], d["atom"]))
           if r == resseq and a == "CA"]
    return d["xyz"][idx[0]] if idx else None


def kabsch_rmsd(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    Pr = (R @ Pc.T).T
    return float(np.sqrt(((Pr - Qc) ** 2).sum(1).mean()))


def contact_pairs(chains_a, xyz_a, xyz_b, cutoff):
    """Return list of (i,j) atom pairs within cutoff, chunked."""
    out = []
    for i0 in range(0, len(xyz_a), 400):
        blk = xyz_a[i0:i0 + 400]
        dd = np.sqrt(((blk[:, None, :] - xyz_b[None, :, :]) ** 2).sum(-1))
        ii, jj = np.where(dd < cutoff)
        out += list(zip(ii + i0, jj))
    return out


# ------------------------------------------------------------- inventory
systems = sorted(d for d in os.listdir(OPENMM)
                 if os.path.isdir(os.path.join(OPENMM, d)))
states = ["A.pdb", "B.pdb", "D.pdb", "mutant.pdb", "minimised.pdb"]
inv = []
for sysname in systems:
    for st in states:
        path = os.path.join(OPENMM, sysname, st)
        if not os.path.exists(path):
            continue
        ch = parse(path)
        nall = sum(len(d["xyz"]) for d in ch.values())
        chh = parse(path, heavy_only=True)
        nheavy = sum(len(d["xyz"]) for d in chh.values())
        comp = ";".join(
            f"{c}:{len(set(d['resseq']))}res[{min(d['resseq'])}-"
            f"{max(d['resseq'])}]" for c, d in sorted(ch.items()))
        inv.append({"system": sysname, "state": st, "n_atoms": nall,
                    "n_heavy_atoms": nheavy,
                    "file_bytes": os.path.getsize(path),
                    "chain_composition": comp})
inv_df = pd.DataFrame(inv)

# ------------------------------------------- mutations vs WT (4U2X numbering)
wt_native = parse(os.path.join(OPENMM, "WT", "mutant.pdb"), heavy_only=True)
wt_seq = {}
for r, t in zip(wt_native["A"]["resseq"], wt_native["A"]["restype"]):
    wt_seq.setdefault(r, t)
muts = {}
for sysname in systems:
    p = os.path.join(OPENMM, sysname, "mutant.pdb")
    if not os.path.exists(p):
        muts[sysname] = "no mutant.pdb"
        continue
    ch = parse(p, heavy_only=True)
    seq = {}
    for r, t in zip(ch["A"]["resseq"], ch["A"]["restype"]):
        seq.setdefault(r, t)
    diff = [f"{wt_seq[r]}{r}{seq[r]}" for r in sorted(set(seq) & set(wt_seq))
            if seq[r] != wt_seq[r]]
    muts[sysname] = ",".join(diff) if diff else "none"
inv_df["mutations_vs_WT_4U2X_numbering"] = inv_df.apply(
    lambda r: muts[r["system"]] if r["state"] == "mutant.pdb" else "", axis=1)
inv_df.to_csv(os.path.join(HERE, "X4_md_system_inventory.tsv"), sep="\t",
              index=False)
print("=== X4 inventory ===")
print(inv_df[["system", "state", "n_atoms", "n_heavy_atoms",
              "chain_composition", "mutations_vs_WT_4U2X_numbering"]]
      .to_string(index=False))
print()

# ------------------------------------------ recomputed metrics (openmm dir)
WT_MIN = parse(os.path.join(OPENMM, "WT", "minimised.pdb"), heavy_only=True)
wt_heavy_a = WT_MIN["A"]
wt_ca = np.array([res_ca(WT_MIN, "A", c) for c in range(62, 76)])  # EBOV 77-90

# native contacts (heavy atoms within 0.8 nm), defined on the WT complex
wt_all = np.vstack([WT_MIN["A"]["xyz"], WT_MIN["B"]["xyz"]])


def res_ids(chains):
    """Ordered {chain: [resseq,...]} in file order (first appearance)."""
    out = {}
    for c, d in chains.items():
        seen = []
        for r in d["resseq"]:
            if not seen or seen[-1] != r:
                seen.append(r)
        out[c] = seen
    return out


def keyed_atoms(chains):
    """{(chain, res_index, atom_name): xyz} using the residue *order*."""
    order = {c: {r: i for i, r in enumerate(ids)}
             for c, ids in res_ids(chains).items()}
    out = {}
    for c, d in chains.items():
        for r, a, xyz in zip(d["resseq"], d["atom"], d["xyz"]):
            out[(c, order[c][r], a)] = xyz
    return out


def residue_contacts(chains, cutoff, pairs_only=("A", "B")):
    """Set of (chain_i,res_i,chain_j,res_j) residue pairs with any heavy-atom
    pair within cutoff (used both for native sets and interface counts)."""
    out = set()
    for ci, cj in [tuple(sorted(pairs_only))]:
        if ci not in chains or cj not in chains:
            continue
        di, dj = chains[ci], chains[cj]
        for i0 in range(0, len(di["xyz"]), 300):
            blk = di["xyz"][i0:i0 + 300]
            dd = np.sqrt(((blk[:, None, :] - dj["xyz"][None, :, :]) ** 2).sum(-1))
            ii, jj = np.where(dd < cutoff)
            for a, b in zip(ii, jj):
                out.add((ci, di["resseq"][i0 + a], cj, dj["resseq"][b]))
    return out


# native contacts are defined on the WT minimised complex
native_res_pairs = residue_contacts(WT_MIN, 8.0, ("A", "B"))
# map the WT native residue pairs onto residue *order* so they can be looked up
# in every system irrespective of the mutation-driven atom/renumbering changes
wt_order = {c: {r: i for i, r in enumerate(res_ids(WT_MIN)[c])}
            for c in ("A", "B")}


def order_pairs(pairs):
    return {(c1, wt_order[c1][r1], c2, wt_order[c2][r2])
            for c1, r1, c2, r2 in pairs}


native_order_pairs = sorted(order_pairs(native_res_pairs))


def q_residuepair(sysname):
    path = os.path.join(OPENMM, sysname, "minimised.pdb")
    if not os.path.exists(path):
        return None, None
    s = parse(path, heavy_only=True)
    ids = res_ids(s)
    if {c: len(v) for c, v in ids.items()} != {c: len(v)
                                              for c, v in res_ids(WT_MIN).items()}:
        return "chain-length-mismatch", len(native_order_pairs)
    coords = {}
    for c, d in s.items():
        m = {}
        for r, xyz in zip(d["resseq"], d["xyz"]):
            m.setdefault(r, []).append(xyz)
        coords[c] = m
    ids_ord = {c: {r: i for i, r in enumerate(ids[c])} for c in ids}
    rinv = {c: {i: r for r, i in ids_ord[c].items()} for c in ids}
    kept = 0
    for c1, i1, c2, i2 in native_order_pairs:
        r1, r2 = rinv[c1][i1], rinv[c2][i2]
        a, b = np.array(coords[c1][r1]), np.array(coords[c2][r2])
        dd = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
        if dd.min() < 12.0:
            kept += 1
    return kept / len(native_order_pairs), len(native_order_pairs)


def iface_contacts(path):
    ch = parse(path, heavy_only=True)
    keys = [k for k in ch if k in ("A", "B", "D")]
    v = [k for k in keys if len(ch[k]["xyz"]) > 100]
    if len(v) < 2:
        return None
    return len(residue_contacts(ch, 4.5, (v[0], v[1])))


rows = []
reported = json.load(open(os.path.join(OPENMM, "openmm_results.json"),
                          encoding="utf-8"))["results"]
for sysname in systems:
    mp = os.path.join(OPENMM, sysname, "minimised.pdb")
    if not os.path.exists(mp):
        rows.append({"system": sysname, "metric": "ca_rmsd_77_90_vs_WT_A",
                     "value": "not_recomputable", "basis_structure": "-",
                     "note": "no minimised.pdb (run failed: "
                             + json.load(open(os.path.join(
                                 OPENMM, "openmm_results.json"),
                                 encoding="utf-8"))["errors"].get(
                                     sysname, "unknown")[:60] + ")"})
        continue
    s = parse(mp, heavy_only=True)
    s_ca = np.array([res_ca(s, "A", c) for c in range(62, 76)])
    if any(x is None for x in s_ca):
        rmsd = "no CA for 77-90"
        rmsd_g = "no CA for 77-90"
    else:
        rmsd = round(kabsch_rmsd(s_ca, wt_ca), 4)
        # global superposition on every shared CA, then RMSD of 77-90
        all_wt = np.array([res_ca(WT_MIN, "A", c) for c in range(1, 217)])
        all_s = np.array([res_ca(s, "A", c) for c in range(1, 217)])
        ok = np.array([a is not None and b is not None
                       for a, b in zip(all_wt, all_s)])
        P = np.array([all_s[i] for i in range(len(ok)) if ok[i]])
        Q = np.array([all_wt[i] for i in range(len(ok)) if ok[i]])
        Pc, Qc = P - P.mean(0), Q - Q.mean(0)
        V, S, Wt_ = np.linalg.svd(Pc.T @ Qc)
        d = np.sign(np.linalg.det(V @ Wt_))
        R = V @ np.diag([1.0, 1.0, d]) @ Wt_
        Pr = (R @ Pc.T).T
        rmsd_g = round(float(np.sqrt(((Pr - Qc) ** 2).sum(1).mean())), 4)
    q, n = q_residuepair(sysname)
    pref = ""
    if sysname in reported:
        pref = (f"openmm_results interaction_kcal_mol="
                f"{reported[sysname]['interaction_kcal_mol']}")
    rows.append({"system": sysname, "metric": "ca_rmsd_77_90_vs_WT_minimised",
                 "value": rmsd, "basis_structure": "minimised.pdb chain A "
                 "(construct 62-75 = EBOV 77-90)",
                 "note": "static single-structure Kabsch RMSD over 77-90"})
    rows.append({"system": sysname, "metric": "ca_rmsd_77_90_vs_WT_globalfit",
                 "value": rmsd_g, "basis_structure": "minimised.pdb chain A, "
                 "superposed on all shared CA",
                 "note": "static; global-fit then RMSD of 77-90"})
    rows.append({"system": sysname, "metric": "native_contact_fraction_residuepair",
                 "value": (round(q, 4) if isinstance(q, float) else q),
                 "basis_structure": f"minimised.pdb vs WT ({n} native residue pairs "
                 "within 0.8 nm, 1.2 nm allowance)",
                 "note": "static; the pilot's Q_fraction is trajectory-based"})
    rows.append({"system": sysname,
                 "metric": "interface_contacts_lt4.5A_before",
                 "value": iface_contacts(os.path.join(OPENMM, sysname,
                                                      "mutant.pdb")),
                 "basis_structure": "mutant.pdb chains A-D",
                 "note": ("openmm_results reports "
                          f"{reported.get(sysname, {}).get('contacts_lt4.5A_before', 'NA')}"
                          if sysname in reported else "not in openmm_results")})
    rows.append({"system": sysname,
                 "metric": "interface_contacts_lt4.5A_after",
                 "value": iface_contacts(mp),
                 "basis_structure": "minimised.pdb chains A-B",
                 "note": ("openmm_results reports "
                          f"{reported.get(sysname, {}).get('contacts_lt4.5A_after', 'NA')}"
                          if sysname in reported else "not in openmm_results")})

# not-recomputable rows for the trajectory-derived quantities
for sysname in systems:
    for m, src in (("trajectory_local_rmsd_nm", "pilot 1 ns trajectories"),
                   ("trajectory_Q_fraction", "pilot 1 ns trajectories"),
                   ("minimised_energy_kcal_mol", "OpenMM force field"),
                   ("void_volume_nm3_within_0.6nm_of_P83_CA",
                    "trajectory frame")):
        rows.append({"system": sysname, "metric": m, "value":
                     "not_recomputable", "basis_structure": "-",
                     "note": f"requires {src}; only static coordinates are "
                             "deposited in this workspace"})

met = pd.DataFrame(rows, columns=["system", "metric", "value",
                                  "basis_structure", "note"])
met.to_csv(os.path.join(HERE, "X4_md_recomputed_metrics.tsv"), sep="\t",
           index=False)
print("=== X4 recomputed metrics (openmm workspace) ===")
print(met.to_string(index=False))

# ------------------------------------------------- p83 stability pilot inventory
print()
print("=== t4_p83_stability_20260925 static structures ===")
for f in ("vp24_WT.pdb", "vp24_AP83.pdb"):
    ch = parse(os.path.join(P83, f))
    print(" ", f, "atoms", sum(len(d["xyz"]) for d in ch.values()),
          "chains", {c: (min(d["resseq"]), max(d["resseq"]),
                         len(set(d["resseq"]))) for c, d in ch.items()})
