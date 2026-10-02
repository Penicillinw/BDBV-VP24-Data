"""T4 structural gate, version 3: fast, cutoff-based, same-structure point mutations.

Changes vs v2: GBn2 with a 1.0 nm non-bonded cutoff instead of NoCutoff
(the NoCutoff GB run was ~O(N^2) on 6000 atoms and did not finish in 15 min),
an orthogonal contact-count readout, and the functional determinant set from
PMID 37243162 (P83S/N135Q/R140H/V141A) added explicitly.

Every query is a point mutation inside the ONE experimental structure 4U2X
(chains A = EBOV VP24, D = human KPNA5).  No cross-structure comparison.
"""

import csv
import json
import os
import time

from openmm import LangevinIntegrator, Platform, app, unit
from openmm.app import ForceField, PDBFile, Simulation
from pdbfixer import PDBFixer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "analysis", "t4_vp24_struct_seqmap_20260925")
WORK = os.path.join(ROOT, "analysis", "t4_vp24_openmm_20260925")

STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL"}
THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
         "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
         "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
         "Y": "TYR", "V": "VAL"}
PLATFORM = Platform.getPlatformByName("CPU")
CUTOFF_NM = 1.0


def extract_complex():
    path = os.path.join(WORK, "complex_wt.pdb")
    seen, by_chain = set(), {"A": [], "D": []}
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


def mutate_text(text, mutations):
    wanted = {(c, n): (w, m) for c, w, n, m in mutations}
    out = []
    for line in text.splitlines():
        if line.startswith("ATOM"):
            key = (line[21], int(line[22:26]))
            if key in wanted:
                if line[12:16].strip() not in ("N", "CA", "C", "O", "CB"):
                    continue
                line = line[:17] + f"{THREE[wanted[key][1]]:>3s}" + line[20:]
        out.append(line)
    return "\n".join(out) + "\n"


def fixer_for(path):
    f = PDBFixer(filename=path)
    f.findMissingResidues()
    f.findMissingAtoms()
    f.addMissingAtoms()
    f.addMissingHydrogens(7.0)
    return f


def make_system(topology):
    return ForceField("amber14-all.xml", "implicit/gbn2.xml").createSystem(
        topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=CUTOFF_NM * unit.nanometer,
        constraints=None,
        removeCMMotion=True,
    )


def minimise(path, out_path, steps=300):
    f = fixer_for(path)
    sim = Simulation(f.topology, make_system(f.topology),
                     LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                        0.002 * unit.picoseconds), PLATFORM)
    sim.context.setPositions(f.positions)
    sim.minimizeEnergy(maxIterations=steps)
    st = sim.context.getState(getPositions=True, energy=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        PDBFile.writeFile(f.topology, st.getPositions(), fh)
    return st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def split_chain(path, chain, out_path):
    lines = [l.rstrip("\r\n") for l in open(path, encoding="utf-8")
             if l.startswith("ATOM") and l[21] == chain]
    lines += ["TER", "END"]
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def chain_ids(path):
    """Chain ids as actually written (pdbfixer renames D -> B, so never hard-code)."""
    import collections
    c = collections.Counter(l[21] for l in open(path, encoding="utf-8")
                            if l.startswith("ATOM"))
    return [ch for ch, _ in c.most_common()]


def single_point(path):
    f = fixer_for(path)
    sim = Simulation(f.topology, make_system(f.topology),
                     LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                        0.002 * unit.picoseconds), PLATFORM)
    sim.context.setPositions(f.positions)
    return sim.context.getState(energy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole)


def interchain_contacts(path, cutoff=4.5):
    """Orthogonal readout: number of heavy-atom contacts between chains A and D."""
    a, d = [], []
    for line in open(path, encoding="utf-8"):
        if not line.startswith("ATOM"):
            continue
        if line[76:78].strip().upper() == "H":
            continue
        rec = (line[21], int(line[22:26]), float(line[30:38]),
               float(line[38:46]), float(line[46:54]))
        (a if line[21] == "A" else d).append(rec)
    pairs = set()
    c2 = cutoff ** 2
    for ca, ra, xa, ya, za in a:
        for cb, rb, xb, yb, zb in d:
            if (xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2 <= c2:
                pairs.add((ra, rb))
    return len(pairs)


def interface_subs():
    rows = list(csv.DictReader(open(os.path.join(SRC, "interface_species.tsv"),
                                    encoding="utf-8"), delimiter="\t"))
    out = {}
    for sp in ["BDBV", "RESTV", "SUDV", "TAFV"]:
        out[sp] = sorted([("A", r["EBOV"], int(r["ebov_full_pos"]), r[sp])
                          for r in rows
                          if r["EBOV"] != r[sp] and r["EBOV"] in THREE and r[sp] in THREE],
                         key=lambda m: m[2])
    return out


def main():
    os.makedirs(WORK, exist_ok=True)
    wt_text = open(extract_complex(), encoding="utf-8").read()
    subs = interface_subs()

    queries = {
        "WT": [],
        "EBOV_N135A": [("A", "N", 135, "A")],
        "EBOV_R140A": [("A", "R", 140, "A")],
        # the published mutant that lost IFN-lambda1 / IFN-beta / ISG15 inhibition
        "EBOV_4x_BDBVres": [("A", "P", 83, "S"), ("A", "N", 135, "Q"),
                            ("A", "R", 140, "H"), ("A", "V", 141, "A")],
        "BDBV_interface5": subs["BDBV"],
        "BDBV_iface5_plus83": sorted(subs["BDBV"] + [("A", "P", 83, "S")],
                                     key=lambda m: m[2]),
    }
    for sp in ["RESTV", "SUDV", "TAFV"]:
        queries[f"{sp}_interface"] = subs[sp]

    results, errors = {}, {}
    for tag, muts in queries.items():
        t0 = time.time()
        try:
            d = os.path.join(WORK, tag)
            os.makedirs(d, exist_ok=True)
            raw = os.path.join(d, "mutant.pdb")
            open(raw, "w", encoding="utf-8").write(mutate_text(wt_text, muts))
            contacts_in = interchain_contacts(raw)
            e_c = minimise(raw, os.path.join(d, "minimised.pdb"))
            mini = os.path.join(d, "minimised.pdb")
            contacts_out = interchain_contacts(mini)
            chains = chain_ids(mini)
            if len(chains) < 2:
                raise RuntimeError(f"expected 2 chains in {mini}, found {chains}")
            split_chain(mini, chains[0], os.path.join(d, "A.pdb"))
            split_chain(mini, chains[1], os.path.join(d, "B.pdb"))
            e_a = single_point(os.path.join(d, "A.pdb"))
            e_d = single_point(os.path.join(d, "B.pdb"))
            results[tag] = {
                "n_mutations": len(muts),
                "mutations": [f"{c}{w}{n}{m}" for c, w, n, m in muts],
                "E_complex": round(e_c, 1),
                "E_chainA": round(e_a, 1),
                "E_chainD": round(e_d, 1),
                "interaction_kcal_mol": round(e_c - e_a - e_d, 2),
                "contacts_lt4.5A_before": contacts_in,
                "contacts_lt4.5A_after": contacts_out,
                "seconds": round(time.time() - t0, 1),
            }
            print(f"{tag:20s} n={len(muts):2d}  E_int="
                  f"{results[tag]['interaction_kcal_mol']:9.2f}  contacts "
                  f"{contacts_in:3d}->{contacts_out:3d}  "
                  f"({results[tag]['seconds']}s)", flush=True)
        except Exception as exc:  # noqa: BLE001
            errors[tag] = f"{type(exc).__name__}: {exc}"
            print(f"{tag:20s} FAILED {errors[tag]}", flush=True)

    if "WT" in results:
        wt = results["WT"]["interaction_kcal_mol"]
        wtc = results["WT"]["contacts_lt4.5A_after"]
        for rec in results.values():
            rec["ddG_vs_WT"] = round(rec["interaction_kcal_mol"] - wt, 2)
            rec["d_contacts_vs_WT"] = rec["contacts_lt4.5A_after"] - wtc

    json.dump(
        {
            "structure": "4U2X chains A (EBOV VP24) + D (human KPNA5)",
            "method": (f"amber14-all + GBn2 implicit solvent, {CUTOFF_NM} nm cutoff; "
                       "300-step minimisation; single-point "
                       "E_complex - E_A - E_D"),
            "design_note": "all queries are point mutations inside ONE experimental "
                           "structure; no cross-structure comparison",
            "controls": {"EBOV_N135A": "published attenuating mutant",
                         "EBOV_R140A": "published attenuating mutant",
                         "EBOV_4x_BDBVres": "published mutant that lost IFN-lambda1, "
                                            "IFN-beta and ISG15 inhibition (PMID 37243162)"},
            "published_ordering_to_reproduce": "eVP24 ~= rVP24 > bVP24",
            "results": results, "errors": errors,
        },
        open(os.path.join(WORK, "openmm_results.json"), "w", encoding="utf-8"),
        ensure_ascii=False, indent=1,
    )
    print("\nwritten", os.path.join(WORK, "openmm_results.json"))


if __name__ == "__main__":
    main()
