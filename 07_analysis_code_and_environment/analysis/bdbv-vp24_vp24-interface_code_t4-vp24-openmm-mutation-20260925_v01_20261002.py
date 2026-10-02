"""T4 structural positive control via OpenMM: point mutations inside ONE structure.

FoldX BuildModel refuses to resolve any residue of 4U2X (verified for the raw
RCSB file, the FoldX-repaired file, a single-chain extract and a renumbered
extract; format `AN135A;` works on other project files, so this is specific to
this structure).  This script therefore uses the *independent* method that the
pre-registration already requires.

Protocol (identical for wild type and every mutant; no cross-structure step):
  1. take chains A (EBOV VP24) and D (human KPNA5) from PDB 4U2X
  2. mutate side chains by truncation at CB (pdbfixer then rebuilds the mutant
     side chain and all hydrogens)
  3. build amber14 + GBn2 implicit solvent, minimise the complex
  4. single-point interaction energy = E(complex) - E(chain A) - E(chain D)
     on the minimised complex coordinates
  5. Delta-Delta-G = interaction(mutant) - interaction(WT)

Positive controls that must behave as published:
  N135A and R140A must destabilise the interface (PMID 37243162, 27974555)
Species ordering to reproduce: EBOV ~ RESTV > BDBV (PMID 27974555)
"""

import csv
import json
import os
import sys
from copy import deepcopy

from openmm import app, unit, LangevinIntegrator, Platform
from openmm.app import ForceField, Modeller, PDBFile, Simulation
from pdbfixer import PDBFixer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "analysis", "t4_vp24_struct_seqmap_20260925")
WORK = os.path.join(ROOT, "analysis", "t4_vp24_openmm_20260925")

STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def extract_complex():
    path = os.path.join(WORK, "complex_wt.pdb")
    seen = set()
    by_chain = {"A": [], "D": []}
    for line in open(os.path.join(SRC, "4U2X.pdb"), encoding="utf-8"):
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        if chain not in ("A", "D"):
            continue
        if line[17:20].strip() not in STANDARD:
            continue
        if line[16] != " ":
            key = (chain, line[22:27], line[12:16])
            if key in seen:
                continue
            seen.add(key)
        by_chain[chain].append(line.rstrip("\r\n"))
    lines = []
    for chain in ("A", "D"):
        lines.extend(by_chain[chain])
        lines.append("TER")
    lines.append("END")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def apply_mutations(pdb_text, mutations):
    """mutations: list of (chain, wt_one_letter, resnum, mut_one_letter)."""
    wanted = {(c, num): (wt, mut) for c, wt, num, mut in mutations}
    out = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "TER", "END")):
            continue
        if line.startswith("ATOM"):
            chain, resnum = line[21], int(line[22:26])
            if (chain, resnum) in wanted:
                wt, mut = wanted[(chain, resnum)]
                atom = line[12:16].strip()
                if atom not in ("N", "CA", "C", "O", "CB"):
                    continue  # drop old side chain, pdbfixer rebuilds it
                line = line[:17] + f"{THREE[mut]:>3s}" + line[20:]
        out.append(line)
    return "\n".join(out) + "\n"


def build_simulation(pdb_path):
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,
    )
    sim = Simulation(fixer.topology, system, LangevinIntegrator(300 * unit.kelvin,
                                                                 1 / unit.picosecond,
                                                                 0.002 * unit.picoseconds))
    sim.context.setPositions(fixer.positions)
    return fixer, sim


def minimise(sim, max_iterations=400):
    sim.minimizeEnergy(maxIterations=max_iterations)
    return sim


def subsystem_energy(topology, positions, keep_chain):
    """Single-point energy of the subset of atoms whose chain id is keep_chain."""
    atoms = [a for a in topology.atoms() if a.residue.chain.id == keep_chain]
    atom_ids = {a.index for a in atoms}
    modeller = Modeller(topology, positions)
    modeller.delete([a for a in topology.atoms() if a.index not in atom_ids])
    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,
    )
    sim = Simulation(modeller.topology, system,
                     LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                        0.002 * unit.picoseconds))
    sim.context.setPositions(modeller.positions)
    return sim.context.getState(energy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole)


def interaction_energy(pdb_path):
    fixer, sim = build_simulation(pdb_path)
    minimise(sim)
    state = sim.context.getState(getPositions=True, energy=True)
    positions = state.getPositions()
    e_complex = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    e_a = subsystem_energy(fixer.topology, positions, "A")
    e_d = subsystem_energy(fixer.topology, positions, "D")
    return e_complex, e_a, e_d, e_complex - e_a - e_d


def interface_substitutions():
    rows = list(csv.DictReader(
        open(os.path.join(SRC, "interface_species.tsv"), encoding="utf-8"), delimiter="\t"))
    out = {}
    for sp in ["BDBV", "SUDV", "TAFV", "RESTV"]:
        muts = []
        for r in rows:
            wt, mut = r["EBOV"], r[sp]
            if wt != mut and wt in THREE and mut in THREE:
                muts.append(("A", wt, int(r["ebov_full_pos"]), mut))
        out[sp] = sorted(muts, key=lambda m: m[2])
    return out


def main():
    os.makedirs(WORK, exist_ok=True)
    wt_path = extract_complex()
    wt_text = open(wt_path, encoding="utf-8").read()

    queries = {
        "WT": [],
        "EBOV_N135A": [("A", "N", 135, "A")],
        "EBOV_R140A": [("A", "R", 140, "A")],
    }
    subs = interface_substitutions()
    for sp, muts in subs.items():
        queries[f"{sp}_interface"] = muts

    results = {}
    for tag, muts in queries.items():
        path = os.path.join(WORK, f"{tag}.pdb")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(apply_mutations(wt_text, muts))
        e_complex, e_a, e_d, e_int = interaction_energy(path)
        results[tag] = {
            "n_mutations": len(muts),
            "mutations": [f"{c}{w}{n}{m}" for c, w, n, m in muts],
            "E_complex": round(e_complex, 2),
            "E_chainA": round(e_a, 2),
            "E_chainD": round(e_d, 2),
            "interaction_kcal_mol": round(e_int, 2),
        }
        print(f"{tag:18s} n_mut={len(muts):2d}  E_int={e_int:9.2f} kcal/mol")

    wt = results["WT"]["interaction_kcal_mol"]
    for tag, rec in results.items():
        rec["ddG_vs_WT"] = round(rec["interaction_kcal_mol"] - wt, 2)

    with open(os.path.join(WORK, "openmm_results.json"), "w", encoding="utf-8") as fh:
        json.dump(
            {
                "structure": "4U2X chains A (EBOV VP24) + D (human KPNA5)",
                "method": "amber14-all + GBn2 implicit solvent, minimised complex, "
                          "single-point interaction energy (E_complex - E_A - E_D)",
                "design_note": "all queries are point mutations inside ONE experimental "
                               "structure; no cross-structure comparison",
                "positive_controls": {
                    "EBOV_N135A": "must be less stable than WT",
                    "EBOV_R140A": "must be less stable than WT",
                },
                "published_ordering_to_reproduce": "eVP24 ~= rVP24 > bVP24",
                "results": results,
            },
            fh, ensure_ascii=False, indent=1,
        )
    print()
    print("written", os.path.join(WORK, "openmm_results.json"))


if __name__ == "__main__":
    main()
