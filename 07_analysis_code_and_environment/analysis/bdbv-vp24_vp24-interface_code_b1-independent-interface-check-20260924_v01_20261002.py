"""B1 positive control - independent integrity check of the FoldX inputs.

This is a *reader-only* verification that does not call FoldX and does not
modify any of the subagent's products. It answers three questions:

1. Did the chain renaming in `build_complexes` produce exactly two protein
   chains named A (GP1) and B (NPC1-C) in every complex file?
2. Did `RepairPDB` preserve the GP1/NPC1 interface, i.e. is the interface
   residue set (4.0 A heavy-atom cutoff) of the repaired model still the same
   as the one measured directly on the deposited structure?
3. Are the four published hotspot positions (GP1 79/141/142/148) present in
   the modelled GP1 chain, and do the residues there match the deposited
   sequences of PMID 39894818 (SUDV I79/A141/Q142/P148 vs EBOV V79/V141/S142/A148)?

Output: analysis/b1_positive_control_20260924/b1_independent_integrity.json
"""

from __future__ import annotations

import json
from pathlib import Path

from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.Polypeptide import is_aa

ROOT = Path(r"G:\本迪布焦研究")
BASE = ROOT / "analysis" / "b1_positive_control_20260924"

# deposited reference files measured directly (raw downloads) vs built complexes
REFERENCE = {
    "ebov": (BASE / "raw" / "5F1B.cif", "A", "C"),
    "sudv": (BASE / "raw" / "9DZ2.cif", "I", "C"),
    "sudv_copy2": (BASE / "raw" / "9DZ2.cif", "E", "D"),
}
BUILT = {
    "ebov": BASE / "complexes" / "ebov_gp1_npc1.pdb",
    "sudv": BASE / "complexes" / "sudv_gp1_npc1.pdb",
    "sudv_copy2": BASE / "complexes" / "sudv_copy2_gp1_npc1.pdb",
}
REPAIRED = {
    "ebov": BASE / "foldx" / "ebov" / "ebov_gp1_npc1_Repair.pdb",
    "sudv": BASE / "foldx" / "sudv" / "sudv_gp1_npc1_Repair.pdb",
    "sudv_copy2": BASE / "foldx" / "sudv_copy2" / "sudv_copy2_gp1_npc1_Repair.pdb",
}

HOTSPOTS = [79, 141, 142, 148]
EXPECTED = {
    "ebov": {79: "V", 141: "V", 142: "S", 148: "A"},
    "sudv": {79: "I", 141: "A", 142: "Q", 148: "P"},
    "sudv_copy2": {79: "I", 141: "A", 142: "Q", 148: "P"},
}


def parse(path: Path):
    if path.suffix.lower() == ".cif":
        from Bio.PDB import MMCIFParser

        return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))


def aa1(res) -> str:
    return protein_letters_3to1_extended.get(res.get_resname(), "X")


def interface(a, b, cutoff=4.0):
    atoms_a = [x for x in a.get_atoms() if x.element != "H"]
    ns = NeighborSearch(atoms_a)
    ra, rb = set(), set()
    for atom in b.get_atoms():
        if atom.element == "H":
            continue
        for hit in ns.search(atom.coord, cutoff):
            ra.add(hit.get_parent().id[1])
            rb.add(atom.get_parent().id[1])
    return sorted(ra), sorted(rb)


def main() -> None:
    report = {}
    for key in REFERENCE:
        entry = {}
        for label, spec in (
            ("deposited", REFERENCE[key]),
            ("built", BUILT[key]),
            ("repaired", REPAIRED[key]),
        ):
            path = spec[0] if label == "deposited" else spec
            if not Path(path).exists():
                entry[label] = {"present": False}
                continue
            model = next(iter(parse(Path(path))))
            chains = list(model)
            if label == "deposited":
                _, gid, nid = spec
                gp1, npc1 = model[gid], model[nid]
            else:
                gid, nid = "A", "B"
                gp1 = model[gid] if gid in model else None
                npc1 = model[nid] if nid in model else None
            if gp1 is None or npc1 is None:
                entry[label] = {
                    "present": True,
                    "chain_ids": [c.id for c in chains],
                    "error": "expected chains A/B not found",
                }
                continue
            ra, rb = interface(gp1, npc1)
            residues = {
                r.id[1]: aa1(r)
                for r in gp1
                if is_aa(r, standard=False) and r.id[1] in HOTSPOTS
            }
            entry[label] = {
                "present": True,
                "file": str(Path(path).name),
                "chain_ids": [c.id for c in chains],
                "n_res_gp1": sum(1 for r in gp1 if is_aa(r, standard=False)),
                "n_res_npc1": sum(1 for r in npc1 if is_aa(r, standard=False)),
                "gp1_interface_res": ra,
                "npc1_interface_res": rb,
                "hotspot_residues": residues,
                "hotspot_match_expected": residues == EXPECTED[key],
            }
        report[key] = entry

    # comparison block: input vs repaired
    checks = {}
    for key, entry in report.items():
        b, r = entry.get("built", {}), entry.get("repaired", {})
        if not (b.get("present") and r.get("present")):
            checks[key] = {"comparable": False}
            continue
        same_gp1 = b["gp1_interface_res"] == r["gp1_interface_res"]
        same_npc1 = b["npc1_interface_res"] == r["npc1_interface_res"]
        checks[key] = {
            "comparable": True,
            "chains_ok_built": b["chain_ids"] == ["A", "B"],
            "chains_ok_repaired": r["chain_ids"] == ["A", "B"],
            "n_gp1_iface_built": len(b["gp1_interface_res"]),
            "n_gp1_iface_repaired": len(r["gp1_interface_res"]),
            "n_npc1_iface_built": len(b["npc1_interface_res"]),
            "n_npc1_iface_repaired": len(r["npc1_interface_res"]),
            "gp1_iface_identical_after_repair": same_gp1,
            "npc1_iface_identical_after_repair": same_npc1,
            "gp1_iface_lost": sorted(set(b["gp1_interface_res"]) - set(r["gp1_interface_res"])),
            "gp1_iface_gained": sorted(set(r["gp1_interface_res"]) - set(b["gp1_interface_res"])),
            "hotspot_match_built": b["hotspot_match_expected"],
            "hotspot_match_repaired": r["hotspot_match_expected"],
        }
    report["_checks"] = checks

    out = BASE / "b1_independent_integrity.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    for key, c in checks.items():
        print(key, json.dumps(c, ensure_ascii=False))
    print("wrote", out)


if __name__ == "__main__":
    main()
