"""T4 structural gate, version 2 (robust): point mutations inside ONE structure.

Why this version exists: FoldX BuildModel cannot resolve any residue of 4U2X
(verified on the raw RCSB file, the FoldX-repaired file, single-chain extract and
a renumbered extract), and the v1 OpenMM script crashed inside Modeller.delete.

Protocol, identical for wild type and every mutant:
  1. chains A (EBOV VP24) + D (human KPNA5) from PDB 4U2X
  2. mutate by truncating the side chain at CB; pdbfixer rebuilds it + hydrogens
  3. amber14-all + GBn2 implicit solvent, minimisation of the complex
  4. write the minimised complex, then split it into chain-A-only and chain-D-only
     FILES (no atom deletion in the OpenMM topology) and evaluate single-point
     energies with the same force field
  5. interaction energy = E(complex) - E(chain A) - E(chain D); ddG vs WT

Positive controls that must behave as published:
  N135A and R140A must destabilise the interface (PMID 37243162 / 27974555).
Species ordering to reproduce: EBOV ~ RESTV > BDBV (PMID 27974555).
"""

import csv
import json
import os
import traceback

from openmm import LangevinIntegrator, Platform, unit
from openmm.app import ForceField, PDBFile, Simulation
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
PLATFORM = Platform.getPlatformByName("CPU")
FF_ARGS = ("amber14-all.xml", "implicit/gbn2.xml")


def extract_complex():
    path = os.path.join(WORK, "complex_wt.pdb")
    seen = set()
    by_chain = {"A": [], "D": []}
    for line in open(os.path.join(SRC, "4U2X.pdb"), encoding="utf-8"):
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        if chain not in ("A", "D") or line[17:20].strip() not in STANDARD:
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


def mutate_text(pdb_text, mutations):
    wanted = {(c, num): (wt, mut) for c, wt, num, mut in mutations}
    out = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            key = (line[21], int(line[22:26]))
            if key in wanted:
                atom = line[12:16].strip()
                if atom not in ("N", "CA", "C", "O", "CB"):
                    continue
                line = line[:17] + f"{THREE[wanted[key][1]]:>3s}" + line[20:]
        out.append(line)
    return "\n".join(out) + "\n"


def load_fixer(path):
    fixer = PDBFixer(filename=path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer


def minimised_complex(path, out_path, steps=300):
    fixer = load_fixer(path)
    system = ForceField(*FF_ARGS).createSystem(
        fixer.topology, nonbondedMethod=app_NoCutoff(), constraints=None)
    sim = Simulation(fixer.topology, system,
                     LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                        0.002 * unit.picoseconds),
                     PLATFORM)
    sim.context.setPositions(fixer.positions)
    sim.minimizeEnergy(maxIterations=steps)
    state = sim.context.getState(getPositions=True, energy=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        PDBFile.writeFile(fixer.topology, state.getPositions(), fh)
    return state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def app_NoCutoff():
    from openmm import app
    return app.NoCutoff


def split_chain(path, chain, out_path):
    lines = []
    for line in open(path, encoding="utf-8"):
        if line.startswith("ATOM") and line[21] == chain:
            lines.append(line.rstrip("\r\n"))
    lines.append("TER")
    lines.append("END")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def single_point(path):
    fixer = load_fixer(path)
    system = ForceField(*FF_ARGS).createSystem(
        fixer.topology, nonbondedMethod=app_NoCutoff(), constraints=None)
    sim = Simulation(fixer.topology, system,
                     LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                        0.002 * unit.picoseconds),
                     PLATFORM)
    sim.context.setPositions(fixer.positions)
    return sim.context.getState(energy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole)


def interface_substitutions():
    rows = list(csv.DictReader(
        open(os.path.join(SRC, "interface_species.tsv"), encoding="utf-8"), delimiter="\t"))
    out = {}
    for sp in ["BDBV", "RESTV", "SUDV", "TAFV"]:
        muts = [("A", r["EBOV"], int(r["ebov_full_pos"]), r[sp])
                for r in rows
                if r["EBOV"] != r[sp] and r["EBOV"] in THREE and r[sp] in THREE]
        out[sp] = sorted(muts, key=lambda m: m[2])
    return out


def run_query(tag, muts, wt_text):
    d = os.path.join(WORK, tag)
    os.makedirs(d, exist_ok=True)
    raw = os.path.join(d, "mutant.pdb")
    with open(raw, "w", encoding="utf-8") as fh:
        fh.write(mutate_text(wt_text, muts))
    minimised = os.path.join(d, "minimised.pdb")
    e_complex = minimised_complex(raw, minimised)
    a_path, d_path = os.path.join(d, "chainA.pdb"), os.path.join(d, "chainD.pdb")
    split_chain(minimised, "A", a_path)
    split_chain(minimised, "D", d_path)
    e_a, e_d = single_point(a_path), single_point(d_path)
    return {
        "n_mutations": len(muts),
        "mutations": [f"{c}{w}{n}{m}" for c, w, n, m in muts],
        "E_complex": round(e_complex, 2),
        "E_chainA": round(e_a, 2),
        "E_chainD": round(e_d, 2),
        "interaction_kcal_mol": round(e_complex - e_a - e_d, 2),
    }


def main():
    os.makedirs(WORK, exist_ok=True)
    wt_path = extract_complex()
    wt_text = open(wt_path, encoding="utf-8").read()

    queries = {
        "WT": [],
        "EBOV_N135A": [("A", "N", 135, "A")],
        "EBOV_R140A": [("A", "R", 140, "A")],
    }
    for sp, muts in interface_substitutions().items():
        queries[f"{sp}_interface"] = muts

    results, errors = {}, {}
    for tag, muts in queries.items():
        try:
            results[tag] = run_query(tag, muts, wt_text)
            print(f"{tag:18s} n_mut={len(muts):2d}  "
                  f"E_int={results[tag]['interaction_kcal_mol']:10.2f} kcal/mol")
        except Exception as exc:  # noqa: BLE001
            errors[tag] = f"{type(exc).__name__}: {exc}"
            print(f"{tag:18s} FAILED  {errors[tag]}")
            traceback.print_exc()

    if "WT" in results:
        wt = results["WT"]["interaction_kcal_mol"]
        for rec in results.values():
            rec["ddG_vs_WT"] = round(rec["interaction_kcal_mol"] - wt, 2)

    with open(os.path.join(WORK, "openmm_results.json"), "w", encoding="utf-8") as fh:
        json.dump(
            {
                "structure": "4U2X chains A (EBOV VP24) + D (human KPNA5)",
                "method": ("amber14-all + GBn2 implicit solvent; complex minimised 300 steps; "
                           "single-point interaction energy E_complex - E_A - E_D"),
                "design_note": ("every query is a point mutation inside ONE experimental "
                                "structure; no cross-structure comparison"),
                "positive_controls": {
                    "EBOV_N135A": "must be less stable than WT",
                    "EBOV_R140A": "must be less stable than WT",
                },
                "published_ordering_to_reproduce": "eVP24 ~= rVP24 > bVP24",
                "results": results,
                "errors": errors,
            },
            fh, ensure_ascii=False, indent=1,
        )
    print("\nwritten", os.path.join(WORK, "openmm_results.json"))


if __name__ == "__main__":
    main()
