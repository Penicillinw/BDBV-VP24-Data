"""T4 structural positive control, done the ONLY way this project allows:
point mutations applied *inside one experimental structure* (PDB 4U2X).

No cross-structure comparison is made, so the failure mode that invalidated the
B1 positive control (different crystal/cryo-EM structures of the same partner
chain differing by 66 kcal/mol in internal stability) cannot apply here.

Pipeline
  1. extract chains A (EBOV VP24) + D (human KPNA5) from 4U2X
  2. FoldX RepairPDB on the complex
  3. FoldX AnalyseComplex on the repaired wild type   -> reference interaction energy
  4. FoldX BuildModel with one mutation list per query (species variants and
     alanine controls), then AnalyseComplex on each mutant
  5. report interaction energy and DeltadeltaG vs wild type

Positive controls that MUST behave as published:
  - AN135A and AR140A must destabilise the interface (PMID 37243162, 27974555)
  - species ordering must be EBOV ~ RESTV > BDBV (PMID 27974555)
"""

import csv
import json
import os
import re
import shutil
import subprocess
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "analysis", "t4_vp24_struct_seqmap_20260925")
WORK = os.path.join(ROOT, "analysis", "t4_vp24_foldx_20260925")
FOLDX = r"F:\FoldX5\FoldX.cmd"

STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def extract_complex():
    """Chains A (VP24) and D (KPNA5), ATOM records only, first altloc."""
    out = os.path.join(WORK, "complex_wt_raw.pdb")
    seen_alt = set()
    by_chain = {"A": [], "D": []}
    for line in open(os.path.join(SRC, "4U2X.pdb"), encoding="utf-8"):
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        if chain not in ("A", "D"):
            continue
        if line[17:20].strip() not in STANDARD:
            continue
        altloc = line[16]
        if altloc != " ":
            key = (chain, line[22:27], line[12:16])
            if key in seen_alt:
                continue
            seen_alt.add(key)
            line = line[:16] + " " + line[17:]
        by_chain[chain].append(line.rstrip("\n"))
    lines = []
    for chain in ("A", "D"):
        lines.extend(by_chain[chain])
        lines.append("TER")
    lines.append("END")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return out


def run_foldx(args, cwd):
    cmd = [FOLDX] + args
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def parse_interaction(fxout_path):
    """Return the FoldX interaction energy (kcal/mol) from an AnalyseComplex output."""
    if not os.path.exists(fxout_path):
        return None
    with open(fxout_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if parts[0].lower().startswith("interaction"):
                m = re.search(r"(-?\d+\.\d+)", line)
                return float(m.group(1)) if m else None
    return None


def interface_substitutions():
    """species -> list of (chain, wt, resnum, mut) using only modelled residues."""
    rows = list(csv.DictReader(
        open(os.path.join(SRC, "interface_species.tsv"), encoding="utf-8"), delimiter="\t"))
    modelled = set()
    for line in open(os.path.join(WORK, "complex_wt.pdb"), encoding="utf-8"):
        if line.startswith("ATOM") and line[21] == "A":
            modelled.add(int(line[22:26]))
    out = {}
    for sp in ["BDBV", "SUDV", "TAFV", "RESTV"]:
        muts = []
        for r in rows:
            pos = int(r["ebov_full_pos"])
            if pos not in modelled:
                continue
            wt, mut = r["EBOV"], r[sp]
            if mut in STANDARD and wt != mut and wt in STANDARD:
                muts.append((wt, pos, mut))
        out[sp] = sorted(muts, key=lambda m: m[1])
    return out


def write_mutant_file(path, chain, muts):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(",".join(f"{chain}{wt}{pos}{mut}" for wt, pos, mut in muts) + ";\n")


def analyse(tag, pdb_path):
    """Run AnalyseComplex in a dedicated subdir and return the interaction energy."""
    d = os.path.join(WORK, tag)
    os.makedirs(d, exist_ok=True)
    shutil.copy(pdb_path, os.path.join(d, os.path.basename(pdb_path)))
    code, log = run_foldx(
        ["--command=AnalyseComplex", f"--pdb={os.path.basename(pdb_path)}"], d)
    with open(os.path.join(d, "analyse.log"), "w", encoding="utf-8") as fh:
        fh.write(log)
    base = os.path.basename(pdb_path)[:-4]
    energy = parse_interaction(os.path.join(d, f"Interaction_{base}_AC.fxout"))
    return energy, d


def build_mutant(tag, muts):
    d = os.path.join(WORK, tag)
    os.makedirs(d, exist_ok=True)
    shutil.copy(os.path.join(WORK, "complex_wt_repaired.pdb"), os.path.join(d, "complex.pdb"))
    write_mutant_file(os.path.join(d, "individual_list.txt"), "A", muts)
    code, log = run_foldx(
        ["--command=BuildModel", "--pdb=complex.pdb",
         "--mutant-file=individual_list.txt", "--numberOfRuns=5",
         "--out-pdb=true"], d)
    with open(os.path.join(d, "build.log"), "w", encoding="utf-8") as fh:
        fh.write(log)
    candidates = [f for f in os.listdir(d) if re.match(r"complex_\d+\.pdb$", f)]
    if not candidates:
        return None, d, "no BuildModel output"
    energies = []
    for cand in sorted(candidates):
        e, _ = analyse(f"{tag}__{cand[:-4]}", os.path.join(d, cand))
        if e is not None:
            energies.append(e)
    if not energies:
        return None, d, "AnalyseComplex returned no interaction term"
    return sum(energies) / len(energies), d, f"n={len(energies)} runs"


def main():
    os.makedirs(WORK, exist_ok=True)
    raw = extract_complex()
    print("extracted complex:", raw)

    code, log = run_foldx(["--command=RepairPDB", "--pdb=complex_wt_raw.pdb"], WORK)
    with open(os.path.join(WORK, "repair.log"), "w", encoding="utf-8") as fh:
        fh.write(log)
    repaired = os.path.join(WORK, "complex_wt_repaired.pdb")
    produced = os.path.join(WORK, "complex_wt_raw_Repair.pdb")
    if not os.path.exists(repaired) and os.path.exists(produced):
        shutil.copy(produced, repaired)
    print("repair exit", code, "| repaired exists:", os.path.exists(repaired))
    if not os.path.exists(repaired):
        print(log[-2000:])
        return

    wt_energy, _ = analyse("WT", repaired)
    print(f"WT interaction energy = {wt_energy} kcal/mol")

    results = OrderedDict()
    results["EBOV_WT"] = {"interaction": wt_energy, "ddg": 0.0, "muts": []}

    controls = {
        "EBOV_N135A": [("N", 135, "A")],
        "EBOV_R140A": [("R", 140, "A")],
    }
    for tag, muts in controls.items():
        e, d, note = build_mutant(tag, muts)
        results[tag] = {"interaction": e, "muts": muts, "note": note,
                        "ddg": None if e is None or wt_energy is None else e - wt_energy}
        print(f"{tag:16s} interaction={e} ddg={results[tag]['ddg']} ({note})")

    subs = interface_substitutions()
    for sp, muts in subs.items():
        tag = f"{sp}_interface"
        e, d, note = build_mutant(tag, muts)
        results[tag] = {"interaction": e, "muts": muts, "note": note,
                        "ddg": None if e is None or wt_energy is None else e - wt_energy}
        print(f"{tag:16s} n_mut={len(muts)} interaction={e} ddg={results[tag]['ddg']} ({note})")

    with open(os.path.join(WORK, "foldx_results.json"), "w", encoding="utf-8") as fh:
        json.dump(
            {
                "structure": "4U2X chains A (EBOV VP24) + D (human KPNA5)",
                "design_note": ("all queries are point mutations inside ONE experimental "
                                "structure; no cross-structure comparison is made"),
                "wt_interaction_kcal_mol": wt_energy,
                "results": results,
                "published_ordering_to_reproduce": "eVP24 ~= rVP24 > bVP24 (~10x lower)",
            },
            fh, ensure_ascii=False, indent=1,
        )
    print()
    print("written:", os.path.join(WORK, "foldx_results.json"))


if __name__ == "__main__":
    main()
