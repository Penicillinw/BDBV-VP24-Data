#!/usr/bin/env python
"""VP24-KPNA1 paralog substitution + a structural audit of every transferred complex.

Part 1  substitution
    There is no experimental VP24-KPNA1 complex, so KPNA1 is placed by
    superposing it onto the KPNA5 chain of 4U2X inside the already-transferred
    VP24-KPNA5 model.  Feasibility is judged by the KPNA1/KPNA5 superposition
    quality (anchors + RMSD) and by the resulting clashes.

Part 2  audit
    Before any transferred model goes into MD it has to pass the same checks a
    crystallographer would run: complete backbone, no chain breaks, no severe
    non-bonded clashes, and an interface that is actually buried.  This prints
    one row per model so the GPU plan can be based on evidence rather than hope.
"""

from __future__ import annotations

import csv
import glob
import os
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import Chain, MMCIFParser, Model, PDBIO, PDBParser, Structure, Superimposer
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
DL = Path("C:/Users/Administrator/Downloads")
TT = ROOT / "data" / "template_transfer_20260919"
ANALYSIS = ROOT / "analysis"

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}
BACKBONE = ("N", "CA", "C", "O")


def load(path, tag):
    if str(path).lower().endswith(".cif"):
        st = MMCIFParser(QUIET=True).get_structure(tag, str(path))
    else:
        st = PDBParser(QUIET=True).get_structure(tag, str(path))
    return next(st.get_models())


def residues(chain):
    return [r for r in chain if r.id[0] == " " and
            r.get_resname().strip().upper() in AA3]


def seq_of(chain):
    return "".join(AA3[r.get_resname().strip().upper()] for r in residues(chain))


def pairs_align(q, r):
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 2, -1
    al.open_gap_score, al.extend_gap_score = -9, -0.5
    a = al.align(r, q)[0]
    out = []
    for (rs, re_), (qs, qe) in zip(a.aligned[0], a.aligned[1]):
        for k in range(min(re_ - rs, qe - qs)):
            out.append((qs + k + 1, rs + k + 1))
    return out


def superpose_onto(fixed_lookup, moving_lookup, fixed_keys, moving_keys):
    """fixed = template, moving = our protein.  Returns (Superimposer, n_anchors)."""
    fx, mv = [], []
    for fk, mk in zip(fixed_keys, moving_keys):
        if fk in fixed_lookup and mk in moving_lookup:
            if "CA" in fixed_lookup[fk] and "CA" in moving_lookup[mk]:
                fx.append(fixed_lookup[fk]["CA"])
                mv.append(moving_lookup[mk]["CA"])
    sup = Superimposer()
    sup.set_atoms(fx, mv)
    return sup, len(fx)


def keys_from_alignment(q_seq, r_seq, q_order, r_order):
    """(ref keys, query keys) from a global alignment of two sequences."""
    p = pairs_align(q_seq, r_seq)          # (query idx, ref idx)
    return ([r_order[i - 1] for _, i in p], [q_order[j - 1] for j, _ in p])


def pick_pair(ids):
    """(viral, host) chain ids.  GP files carry G/P/N, the others V/H."""
    if "V" in ids and "H" in ids:
        return "V", "H"
    if "G" in ids and "N" in ids:          # GP1 / GP2 / NPC1
        return "G", "N"
    return ids[0], ids[1]


def audit(path: Path):
    """Geometry + interface checks for one two-chain complex."""
    m = next(PDBParser(QUIET=True).get_structure("x", str(path)).get_models())
    ids = [c.id for c in m]
    out = {"file": path.name, "chains": "/".join(ids)}
    if len(ids) < 2:
        out["error"] = "fewer than 2 chains"
        return out
    vid, hid = pick_pair(ids)
    out["viral_chain"], out["host_chain"] = vid, hid
    a, b = residues(m[vid]), residues(m[hid])
    out["n_res_V"], out["n_res_H"] = len(a), len(b)

    # incomplete backbone
    miss = sum(1 for r in a + b for at in BACKBONE if at not in r)
    out["missing_backbone_atoms"] = miss

    # chain breaks (CA-CA of consecutive residues > 4.5 A)
    def breaks(res):
        n = 0
        for i in range(1, len(res)):
            if "CA" in res[i] and "CA" in res[i - 1]:
                d = np.linalg.norm(res[i]["CA"].coord - res[i - 1]["CA"].coord)
                if d > 4.5:
                    n += 1
        return n
    out["chain_breaks_V"], out["chain_breaks_H"] = breaks(a), breaks(b)

    # clashes and interface
    ca = np.array([at.coord for r in a for at in r if at.element != "H"])
    cb = np.array([at.coord for r in b for at in r if at.element != "H"])
    d, _ = cKDTree(cb).query(ca, k=1)
    out["clash_atoms_lt2A"] = int((d < 2.0).sum())
    out["clash_atoms_lt2.5A"] = int((d < 2.5).sum())
    out["iface_atoms_lt4A"] = int((d < 4.0).sum())
    out["iface_atoms_lt8A"] = int((d < 8.0).sum())

    # interface residues
    owner = np.repeat(np.arange(len(a)),
                      [len([x for x in r if x.element != "H"]) for r in a])
    out["iface_res_V"] = len(set(owner[d < 8.0].tolist()))
    ownerb = np.repeat(np.arange(len(b)),
                       [len([x for x in r if x.element != "H"]) for r in b])
    db, _ = cKDTree(ca).query(cb, k=1)
    out["iface_res_H"] = len(set(ownerb[db < 8.0].tolist()))

    # buried check: does the partner actually cover the interface?
    out["min_dist_A"] = round(float(d.min()), 2)
    return out


def main() -> None:
    print("=" * 78)
    print("Part 1  VP24-KPNA1 paralog substitution")
    print("=" * 78)
    tpl = load(DL / "4U2X.cif", "4u2x")
    vp24_tpl = residues(tpl["A"])
    kpna5 = residues(tpl["D"])
    kpna5_res = {r.id[1]: r for r in kpna5}
    tpl_vp24_res = {r.id[1]: r for r in vp24_tpl}
    kpna5_seq = seq_of(tpl["D"])

    made = []
    for species in ("bdbv", "ebov"):
        g = load(ROOT / f"data/gate0_harmonized/VP24_KPNA1/{species}_complex_min.pdb",
                 species)
        vp24 = residues(g["A"])
        kpna1 = residues(g["B"])
        vp24_res = {r.id[1]: r for r in vp24}
        kpna1_res = {r.id[1]: r for r in kpna1}
        vp24_seq, kpna1_seq = seq_of(g["A"]), seq_of(g["B"])

        # (1) superpose our VP24 onto the template VP24, transfer KPNA5
        tpl_vp24_seq = seq_of(tpl["A"])
        tpl_order = [r.id[1] for r in vp24_tpl]
        our_order = [r.id[1] for r in vp24]
        key_fixed, key_moving = keys_from_alignment(
            vp24_seq, tpl_vp24_seq, our_order, tpl_order)
        sup1, n1 = superpose_onto(tpl_vp24_res, vp24_res,
                                  key_fixed, key_moving)
        rot, tran = sup1.rotran
        rot_inv = np.asarray(rot).T
        tran_inv = -np.asarray(tran) @ rot_inv
        kpna5_moved = tpl["D"].copy()
        kpna5_moved.id = "H"
        kpna5_moved.transform(rot_inv, tran_inv)
        print(f"\n{species}: VP24 onto 4U2X VP24 -> {n1} anchors, "
              f"RMSD {sup1.rms:.2f} A")

        # (2) superpose KPNA1 onto the transferred KPNA5
        moved5 = {r.id[1]: r for r in residues(kpna5_moved)}
        k5_order = [r.id[1] for r in residues(kpna5_moved)]
        k1_order = [r.id[1] for r in kpna1]
        key5, key1 = keys_from_alignment(kpna1_seq, kpna5_seq,
                                         k1_order, k5_order)
        sup2, n2 = superpose_onto(moved5, kpna1_res, key5, key1)
        pk = pairs_align(kpna1_seq, kpna5_seq)
        ident = sum(1 for qi, ri in pk
                    if kpna1_seq[qi - 1] == kpna5_seq[ri - 1])
        print(f"  KPNA1 onto KPNA5 -> {n2} anchors, RMSD {sup2.rms:.2f} A, "
              f"identity {ident / len(pk):.0%} over {len(pk)} aligned residues")
        r2, t2 = sup2.rotran
        kpna1_moved = g["B"].copy()
        kpna1_moved.id = "H"
        kpna1_moved.transform(np.asarray(r2), np.asarray(t2))

        # (3) emit
        st = Structure.Structure(f"{species}_vp24_kpna1_paralog")
        mo = Model.Model(0)
        st.add(mo)
        cv = Chain.Chain("V")
        for r in vp24:
            cv.add(r.copy())
        mo.add(cv)
        ch = Chain.Chain("H")
        for r in residues(kpna1_moved):
            ch.add(r.copy())
        mo.add(ch)
        out = TT / f"{species}_vp24_kpna1_from4u2x_paralog.pdb"
        io = PDBIO()
        io.set_structure(st)
        io.save(str(out))
        made.append(out)
        print(f"  -> {out.name}")

    print()
    print("=" * 78)
    print("Part 2  structural audit of every transferred complex")
    print("=" * 78)
    rows = []
    for p in sorted(glob.glob(str(TT / "*.pdb"))):
        try:
            rows.append(audit(Path(p)))
        except Exception as exc:  # noqa: BLE001
            rows.append({"file": os.path.basename(p), "error": str(exc)})

    hdr = ["file", "n_res_V", "n_res_H", "missing_backbone_atoms",
           "chain_breaks_V", "chain_breaks_H", "clash_atoms_lt2A",
           "iface_res_V", "iface_res_H", "min_dist_A"]
    print(f"{'file':48s} {'V':>4s} {'H':>5s} {'missBB':>7s} {'brk':>7s} "
          f"{'clash<2A':>9s} {'ifaceV':>7s} {'ifaceH':>7s} {'dmin':>6s}")
    print("-" * 110)
    for r in rows:
        if "error" in r:
            print(f"{r['file']:48s}  ERROR {r['error']}")
            continue
        print(f"{r['file']:48s} {r['n_res_V']:4d} {r['n_res_H']:5d} "
              f"{r['missing_backbone_atoms']:7d} "
              f"{r['chain_breaks_V']:3d}/{r['chain_breaks_H']:<3d} "
              f"{r['clash_atoms_lt2A']:9d} {r['iface_res_V']:7d} "
              f"{r['iface_res_H']:7d} {r['min_dist_A']:6.2f}")

    out = ANALYSIS / "template_transfer_audit_20260919.tsv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr +
                           ["clash_atoms_lt2.5A", "iface_atoms_lt4A",
                            "iface_atoms_lt8A", "chains"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
