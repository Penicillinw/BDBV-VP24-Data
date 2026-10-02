"""Map the EBOV VP24 - KPNA5 interface (PDB 4U2X) onto the orthoebolavirus MSA.

4U2X contains three cognate copies of the complex in the asymmetric unit:
  entity 1 = EBOV VP24 (Mayinga 1976), chains A/B/C, 216 resolved residues
  entity 2 = human importin subunit alpha-6 (= KPNA5), chains D/E/F, 175 aa
so cognate pairs are A-D, B-E, C-F.

Outputs (analysis/t4_vp24_struct_seqmap_20260925/):
  4U2X.pdb                 - complex structure
  chain_summary.tsv        - chains, lengths
  interface_contacts.tsv   - cognate-pair residue contacts within the cutoff
  interface_map.tsv        - per EBOV VP24 interface residue: MSA column and the
                             residue in BDBV / SUDV / TAFV / RESTV
"""

import csv
import os
import urllib.request
from collections import OrderedDict, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "t4_vp24_struct_seqmap_20260925")
MSA = os.path.join(ROOT, "data", "orthologs", "msa_clustalw", "VP24.aln.fasta")
PDB_ID = "4U2X"
CUTOFF = 4.5  # A, heavy atom


def fetch(url, path):
    if not os.path.exists(path):
        with urllib.request.urlopen(url, timeout=120) as resp:
            data = resp.read()
        with open(path, "wb") as fh:
            fh.write(data)
    return open(path, "rb").read().decode("utf-8", "replace")


def parse_pdb(text):
    chains = defaultdict(list)
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if element == "H":
                continue
            chains[line[21]].append(
                {
                    "resname": line[17:20].strip(),
                    "resseq": int(line[22:26]),
                    "icode": line[26],
                    "atom": line[12:16].strip(),
                    "xyz": (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ),
                }
            )
    return chains


def read_fasta(path):
    seqs = OrderedDict()
    name, chunks = None, []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(chunks)
            name, chunks = line[1:].strip(), []
        elif line:
            chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def main():
    os.makedirs(OUT, exist_ok=True)
    text = fetch(f"https://files.rcsb.org/download/{PDB_ID}.pdb", os.path.join(OUT, f"{PDB_ID}.pdb"))
    chains = parse_pdb(text)

    order = sorted(chains)
    with open(os.path.join(OUT, "chain_summary.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["chain", "n_atoms", "n_residues", "first_resseq", "last_resseq"])
        for ch in order:
            residues = {a["resseq"] for a in chains[ch]}
            writer.writerow([ch, len(chains[ch]), len(residues), min(residues), max(residues)])
    print("chains in", PDB_ID, {ch: len({a['resseq'] for a in chains[ch]}) for ch in order})

    cognate = [("A", "D"), ("B", "E"), ("C", "F")]
    rows = []
    best = {}
    for ca, cb in cognate:
        for a in chains[ca]:
            for b in chains[cb]:
                dx = a["xyz"][0] - b["xyz"][0]
                dy = a["xyz"][1] - b["xyz"][1]
                dz = a["xyz"][2] - b["xyz"][2]
                d2 = dx * dx + dy * dy + dz * dz
                if d2 <= CUTOFF * CUTOFF:
                    rows.append(
                        [ca, a["resname"], a["resseq"], cb, b["resname"], b["resseq"],
                         a["atom"], b["atom"], round(d2 ** 0.5, 3)]
                    )
                    key = (a["resseq"], b["resseq"])
                    d = d2 ** 0.5
                    if key not in best or d < best[key]:
                        best[key] = d
    with open(os.path.join(OUT, "interface_contacts.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["chain_vp24", "res_vp24", "seq_vp24", "chain_kpna",
                         "res_kpna", "seq_kpna", "atom_vp24", "atom_kpna", "dist_A"])
        writer.writerows(rows)

    # VP24 interface residues, with per-copy reproducibility
    from collections import Counter
    vp24_res = Counter()
    for r in rows:
        vp24_res[(r[1], r[2])] += 1
    per_copy = defaultdict(set)
    for r in rows:
        per_copy[r[0]].add(r[2])
    print(f"cognate contact atom pairs: {len(rows)}; unique residue pairs: {len(best)}")
    print(f"VP24 interface residues (union over copies): {len({k[1] for k in vp24_res})}")
    print("per-copy interface residue counts:", {k: len(v) for k, v in sorted(per_copy.items())})
    shared = set.intersection(*per_copy.values())
    print(f"interface residues present in ALL three copies: {len(shared)}")
    with open(os.path.join(OUT, "vp24_interface_residues.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["vp24_resseq", "vp24_resname", "n_contact_atoms_total",
                         "in_copy_A", "in_copy_B", "in_copy_C"])
        for (resname, resseq), n in sorted(vp24_res.items(), key=lambda kv: kv[0][1]):
            writer.writerow([resseq, resname, n,
                             resseq in per_copy["A"], resseq in per_copy["B"], resseq in per_copy["C"]])
    return chains, shared


if __name__ == "__main__":
    main()
