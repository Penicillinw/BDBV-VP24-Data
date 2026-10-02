"""T4 pilot: is P83S a VP24 *stability* substitution?  (local feasibility run)

Hypothesis under test
  BDBV carries S83 where EBOV carries P83.  In PDB 4U2X, P83 sits in a tight
  turn (N82 at 1.35 A, N84 at 1.33 A) and is 11.4 A away from KPNA5, i.e. it is
  NOT a direct interface contact.  Proline is the residue that best tolerates a
  turn backbone, so P83S is expected to destabilise the local fold rather than
  the interface.  PMID 27974555 independently reported that weaker KPNA binding
  goes together with a shorter VP24 half-life, i.e. a stability component.

What this script does (pilot scale, chosen to confirm feasibility on this box)
  A. deterministic local packing analysis: void volume in a sphere around the
     P83 CA, wild type vs P83S, after identical minimisation
  B. short MD: 2 replicas x N ps at 300 K + 1 heat ramp, implicit solvent,
     OpenCL platform; metrics = local RMSD (77-90), native-contact fraction,
     radius of gyration
  C. throughput report (ns/day) so the scale-up on a GTX 5090 can be planned

Usage
  python tools/t4_p83_stability_pilot_20260925.py bench
  python tools/t4_p83_stability_pilot_20260925.py run [ps_per_replica] [replicas]
"""

import json
import math
import os
import sys
import time

from openmm import LangevinIntegrator, Platform, unit
from openmm.app import ForceField, PDBFile, Simulation
from pdbfixer import PDBFixer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "analysis", "t4_vp24_struct_seqmap_20260925", "4U2X.pdb")
WORK = os.path.join(ROOT, "analysis", "t4_p83_stability_20260925")

STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL"}
THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
         "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
         "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
         "Y": "TYR", "V": "VAL"}
PLATFORM = Platform.getPlatformByName("OpenCL")
LOOP = {str(i) for i in range(77, 91)}  # keys come from residue.id, i.e. strings


def write_monomer(mutation):
    """Chain A (EBOV VP24) only; optional single point mutation."""
    os.makedirs(WORK, exist_ok=True)
    tag = "WT" if mutation is None else f"{mutation[0]}{mutation[1]}{mutation[2]}"
    path = os.path.join(WORK, f"vp24_{tag}.pdb")
    seen, lines = set(), []
    for line in open(SRC, encoding="utf-8"):
        if not line.startswith("ATOM") or line[21] != "A":
            continue
        if line[17:20].strip() not in STANDARD:
            continue
        if line[16] != " ":
            key = (line[22:27], line[12:16])
            if key in seen:
                continue
            seen.add(key)
        if mutation is not None:
            _, wt, resnum, mut = mutation
            if int(line[22:26]) == resnum:
                if line[12:16].strip() not in ("N", "CA", "C", "O", "CB"):
                    continue
                line = line[:17] + f"{THREE[mut]:>3s}" + line[20:]
        lines.append(line.rstrip("\r\n"))
    lines += ["TER", "END"]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


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
        nonbondedMethod=__import__("openmm").app.CutoffNonPeriodic,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=__import__("openmm").app.HBonds,
        hydrogenMass=1.5 * unit.amu,
        removeCMMotion=True,
    )


def atoms_by_residue(topology, positions):
    out = {}
    for atom in topology.atoms():
        if atom.residue.name == "HOH":
            continue
        out.setdefault(atom.residue.id, []).append(
            (atom.element.symbol, positions[atom.index].value_in_unit(unit.nanometer)))
    return out


def void_volume_near(residues, center_resnum, radius=0.6, grid=0.25, probe=0.14):
    """Grid-based void probe in a sphere around the CA of `center_resnum` (nm)."""
    center = None
    for sym, xyz in residues.get(center_resnum, []):
        if sym == "C":
            center = xyz
            break
    if center is None:
        return None
    vdw = {"C": 0.17, "N": 0.155, "O": 0.152, "S": 0.18, "H": 0.12}
    pts = []
    steps = int(radius / grid)
    for i in range(-steps, steps + 1):
        for j in range(-steps, steps + 1):
            for k in range(-steps, steps + 1):
                p = (center[0] + i * grid, center[1] + j * grid, center[2] + k * grid)
                if math.dist(p, center) > radius:
                    continue
                pts.append(p)
    free = 0
    for p in pts:
        if all(
            math.dist(p, xyz) > vdw.get(sym, 0.17) + probe
            for res in residues.values() for sym, xyz in res
        ):
            free += 1
    return free * grid ** 3


def native_contacts(residues):
    keys = {}
    ids = sorted(residues)
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            ca, cb = ids[a], ids[b]
            if abs(int(ca) - int(cb)) <= 3:
                continue
            d = min((math.dist(x, y)
                     for _, x in residues[ca] for _, y in residues[cb]), default=99)
            if d < 0.8:
                keys[(ca, cb)] = d
    return keys


def q_fraction(reference, current):
    ok = 0
    for (a, b), d0 in reference.items():
        try:
            d = min((math.dist(x, y) for _, x in current.get(a, [])
                     for _, y in current.get(b, [])), default=99)
        except ValueError:
            d = 99
        if d < d0 * 1.5:
            ok += 1
    return ok / max(1, len(reference))


def local_rmsd(reference, current):
    num, den = 0.0, 0
    for res in sorted(LOOP):
        ca0 = next((x for s, x in reference.get(res, []) if s == "C"), None)
        ca1 = next((x for s, x in current.get(res, []) if s == "C"), None)
        if ca0 and ca1:
            num += math.dist(ca0, ca1) ** 2
            den += 1
    return math.sqrt(num / den) if den else None


def snapshot(residues):
    return {k: [(s, tuple(float(c) for c in v)) for s, v in vals]
            for k, vals in residues.items()}


def prepare(mutation):
    path = write_monomer(mutation)
    f = fixer_for(path)
    sim = Simulation(f.topology, make_system(f.topology),
                     LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                        0.004 * unit.picoseconds), PLATFORM)
    sim.context.setPositions(f.positions)
    sim.minimizeEnergy(maxIterations=500)
    st = sim.context.getState(getPositions=True, energy=True)
    pos = st.getPositions()
    residues = atoms_by_residue(f.topology, pos)
    return sim, f.topology, pos, snapshot(residues), st


def run_traj(sim, tag, ps, temp_profile, sample_ps=20.0):
    """temp_profile: list of (temperature_K, ps) segments."""
    integ = sim.context.getIntegrator()
    integ.setStepSize(0.004 * unit.picoseconds)
    total_steps = int(ps / 0.004)
    interval = max(1, int(sample_ps / 0.004))
    samples = []
    t0 = time.time()
    done = 0
    for temp, seg_ps in temp_profile:
        integ.setTemperature(temp * unit.kelvin)
        steps = int(seg_ps / 0.004)
        sim.step(steps)
        done += steps
        st = sim.context.getState(getPositions=True)
        samples.append({"temp_K": temp,
                        "residues": atoms_by_residue(sim.topology, st.getPositions())})
    wall = time.time() - t0
    return samples, wall, done * 0.004  # ps simulated


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "bench"
    os.makedirs(WORK, exist_ok=True)
    summary = {"platform": PLATFORM.getName(), "mode": mode, "variants": {}}

    for label, mutation in [("WT", None), ("P83S", ("A", "P", 83, "S"))]:
        sim, topology, pos, ref, st = prepare(mutation)
        e_min = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        void = void_volume_near(ref, "83")
        contacts = native_contacts(ref)
        summary["variants"][label] = {
            "minimised_energy_kcal_mol": round(e_min, 1),
            "void_volume_nm3_within_0.6nm_of_P83_CA": round(void, 4),
            "n_native_contacts": len(contacts),
        }
        print(f"[{label}] E_min={e_min:9.1f} kcal/mol  void={void:.4f} nm^3  "
              f"native contacts={len(contacts)}", flush=True)

        if mode == "bench":
            samples, wall, sim_ps = run_traj(sim, label, 20.0, [(300, 20.0)])
            ns_per_day = (sim_ps / 1000.0) / (wall / 86400.0)
            summary["variants"][label]["bench_ns_per_day"] = round(ns_per_day, 3)
            print(f"[{label}] benchmark: {sim_ps:.0f} ps in {wall:.0f} s "
                  f"-> {ns_per_day:.2f} ns/day", flush=True)
        else:
            reps = int(sys.argv[3]) if len(sys.argv) > 3 else 2
            ps = float(sys.argv[2]) if len(sys.argv) > 2 else 200.0
            for r in range(reps):
                samples, wall, sim_ps = run_traj(sim, label, ps, [(300, ps)])
                last = samples[-1]["residues"]
                snap = snapshot(last)
                summary["variants"][label].setdefault("replicas", []).append({
                    "replica": r,
                    "ps": sim_ps,
                    "wall_s": round(wall, 1),
                    "local_rmsd_nm": round(local_rmsd(ref, snap) or 0, 4),
                    "Q_fraction": round(q_fraction(contacts, snap), 3),
                })
                print(f"[{label}] replica {r}: RMSD(77-90)="
                      f"{local_rmsd(ref, snap):.3f} nm  Q={q_fraction(contacts, snap):.3f}  "
                      f"({wall:.0f}s)", flush=True)

    with open(os.path.join(WORK, f"pilot_{mode}.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print("\nwritten", os.path.join(WORK, f"pilot_{mode}.json"))


if __name__ == "__main__":
    main()
