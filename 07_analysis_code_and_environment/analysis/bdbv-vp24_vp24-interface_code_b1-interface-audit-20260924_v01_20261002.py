"""B1 positive control - audit: protein-only interfaces and NPC1 construct identity.

Motivation: the first-pass interface table counted *all* residues within 4.0 A,
including covalently attached N-glycans (NAG etc.) on the EBOV 5F1B structure.
That inflates the apparent EBOV interface and must be corrected before any
"the interfaces differ in size" statement.

Outputs:
  analysis/b1_positive_control_20260924/b1_interface_audit.json
"""

from __future__ import annotations

import json
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

RAW = Path(r"G:\本迪布焦研究\analysis\b1_positive_control_20260924\raw")
OUT = Path(r"G:\本迪布焦研究\analysis\b1_positive_control_20260924")
CUT = 4.0

PAIRS = {
    "EBOV_5F1B_GP1A_NPC1C": ("5F1B", "A", "C"),
    "SUDV_9DZ2_GP1I_NPC1C": ("9DZ2", "I", "C"),
    "SUDV_9DZ2_GP1E_NPC1D": ("9DZ2", "E", "D"),
}


def seq_of(chain) -> str:
    return "".join(
        protein_letters_3to1_extended.get(r.get_resname(), "X")
        for r in chain
        if is_aa(r, standard=False)
    )


def interface(gp1, npc1):
    """Return (gp1 protein residues, npc1 protein residues, non-protein contacts)."""
    all_gp1 = [a for a in gp1.get_atoms() if a.element != "H"]
    ns = NeighborSearch(all_gp1)
    res_gp1, res_npc1, other = set(), set(), set()
    for atom in npc1.get_atoms():
        if atom.element == "H":
            continue
        hits = list(ns.search(atom.coord, CUT))
        if not hits:
            continue
        parent = atom.get_parent()
        if is_aa(parent, standard=False):
            res_npc1.add(parent.id[1])
        else:
            other.add((parent.get_resname(), parent.id[1]))
        for hit in hits:
            hp = hit.get_parent()
            if is_aa(hp, standard=False):
                res_gp1.add(hp.id[1])
            else:
                other.add((hp.get_resname(), hp.id[1]))
    return sorted(res_gp1), sorted(res_npc1), sorted(other)


def main() -> None:
    structs = {}
    for stem in {v[0] for v in PAIRS.values()}:
        structs[stem] = MMCIFParser(QUIET=True).get_structure(stem, str(RAW / f"{stem}.cif"))

    report = {}
    for label, (stem, gp1_id, npc1_id) in PAIRS.items():
        model = next(iter(structs[stem]))
        gp1_res, npc1_res, other = interface(model[gp1_id], model[npc1_id])
        report[label] = {
            "structure": stem,
            "gp1_chain": gp1_id,
            "npc1_chain": npc1_id,
            "gp1_protein_interface_residues": gp1_res,
            "npc1_protein_interface_residues": npc1_res,
            "n_gp1": len(gp1_res),
            "n_npc1": len(npc1_res),
            "non_protein_contacts": [list(x) for x in other],
        }

    # are the two NPC1 domain-C chains the same construct?
    n1 = seq_of(next(iter(structs["5F1B"]))["C"])
    n2 = seq_of(next(iter(structs["9DZ2"]))["C"])
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aln = aligner.align(n1, n2)[0]
    ident = sum(1 for a, b in zip(aln[0], aln[1]) if a == b and a != "-")
    report["npc1_construct_check"] = {
        "ebov_5f1b_npc1_len": len(n1),
        "sudv_9dz2_npc1_len": len(n2),
        "alignment_length": len(aln[0]),
        "identity": ident,
        "identity_pct": round(100 * ident / len(aln[0]), 2),
        "ebov_prefix": n1[:30],
        "sudv_prefix": n2[:30],
    }

    (OUT / "b1_interface_audit.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    for label, d in report.items():
        if label == "npc1_construct_check":
            print("NPC1 construct:", d)
            continue
        print(f"=== {label}")
        print(f"    GP1  protein-only: n={d['n_gp1']} {d['gp1_protein_interface_residues']}")
        print(f"    NPC1 protein-only: n={d['n_npc1']} {d['npc1_protein_interface_residues']}")
        print(f"    non-protein contacts: {d['non_protein_contacts']}")


if __name__ == "__main__":
    main()
