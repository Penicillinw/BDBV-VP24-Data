"""Map the PDB 4U2X VP24 interface residues onto full-length VP24 numbering
and then onto the five-species orthoebolavirus MSA.

Inputs : analysis/t4_vp24_struct_seqmap_20260925/vp24_interface_residues.tsv
         analysis/t4_vp24_struct_seqmap_20260925/4U2X.pdb
         data/orthologs/msa_clustalw/VP24.aln.fasta
Outputs: interface_numbering.tsv  - structure numbering -> EBOV full-length numbering
         interface_species.tsv    - per interface residue: EBOV vs BDBV/SUDV/TAFV/RESTV
         summary.json
"""

import csv
import json
import os
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT_DIR = os.path.join(ROOT, "analysis", "t4_vp24_struct_seqmap_20260925")
MSA_PATH = os.path.join(ROOT, "data", "orthologs", "msa_clustalw", "VP24.aln.fasta")

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


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


def chain_sequence(pdb_path, chain):
    residues = OrderedDict()
    for line in open(pdb_path, encoding="utf-8"):
        if line.startswith("ATOM") and line[21] == chain:
            resname = line[17:20].strip()
            resseq = int(line[22:26])
            residues.setdefault(resseq, THREE2ONE.get(resname, "X"))
    return list(residues.items())


def best_offset(struct_seq, ref_seq):
    """Slide struct_seq over ref_seq (ungapped) and return (offset, matches, coverage)."""
    best = (-1, -1)
    best_off = None
    for off in range(0, len(ref_seq) - len(struct_seq) + 1):
        matches = sum(1 for i, ch in enumerate(struct_seq) if ref_seq[off + i] == ch)
        if matches > best[1]:
            best = (off, matches)
            best_off = off
    return best_off, best[1]


def main():
    msa = read_fasta(MSA_PATH)
    species = {}
    for header, seq in msa.items():
        parts = header.split("|")
        species[parts[1]] = seq.upper()
    ebov_aln = species["EBOV"]
    ebov_ungapped = ebov_aln.replace("-", "")

    res = list(csv.DictReader(
        open(os.path.join(STRUCT_DIR, "vp24_interface_residues.tsv"), encoding="utf-8"),
        delimiter="\t",
    ))
    chain = chain_sequence(os.path.join(STRUCT_DIR, "4U2X.pdb"), "A")
    struct_positions = [p for p, _ in chain]
    struct_seq = "".join(aa for _, aa in chain)

    offset, matches = best_offset(struct_seq, ebov_ungapped)
    print(f"structure chain A: {len(struct_seq)} aa, resseq {struct_positions[0]}..{struct_positions[-1]}")
    print(f"best placement in full-length EBOV VP24 (251 aa): offset={offset}, "
          f"identity={matches}/{len(struct_seq)} = {matches/len(struct_seq):.3f}")

    # MSA column lookup: ungapped EBOV index (1-based) -> alignment column
    col_of = {}
    idx = 0
    for col, ch in enumerate(ebov_aln):
        if ch != "-":
            idx += 1
            col_of[idx] = col

    rows = []
    for r in res:
        s_pos = int(r["vp24_resseq"])
        if s_pos not in struct_positions:
            continue
        full_pos = offset + struct_positions.index(s_pos) + 1
        col = col_of.get(full_pos)
        row = {
            "struct_resseq": s_pos,
            "struct_resname": struct_seq[struct_positions.index(s_pos)],
            "ebov_full_pos": full_pos,
            "n_contact_atoms": int(r["n_contact_atoms_total"]),
            "in_all_copies": (r["in_copy_A"] == "True" and r["in_copy_B"] == "True"
                              and r["in_copy_C"] == "True"),
        }
        for sp in ["EBOV", "BDBV", "SUDV", "TAFV", "RESTV"]:
            seq = species.get(sp)
            row[sp] = seq[col] if (seq and col is not None) else "NA"
        row["bdbv_differs"] = row["BDBV"] not in (row["EBOV"], "-", "NA")
        rows.append(row)

    rows.sort(key=lambda r: (-r["n_contact_atoms"], r["ebov_full_pos"]))
    with open(os.path.join(STRUCT_DIR, "interface_species.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    diffs = [r for r in rows if r["bdbv_differs"]]
    summary = {
        "pdb": "4U2X",
        "structure_chain": "A (EBOV VP24)",
        "structure_length": len(struct_seq),
        "placement_offset": offset,
        "placement_identity": round(matches / len(struct_seq), 4),
        "n_interface_residues": len(rows),
        "n_in_all_three_copies": sum(1 for r in rows if r["in_all_copies"]),
        "bdbv_differing_interface_residues": [
            {k: r[k] for k in ("ebov_full_pos", "EBOV", "BDBV", "SUDV", "TAFV", "RESTV",
                               "n_contact_atoms", "in_all_copies")}
            for r in diffs
        ],
    }
    with open(os.path.join(STRUCT_DIR, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    print()
    print(f"{'EBOV#':>6s} {'ebov':>5s} {'bdbv':>5s} {'sudv':>5s} {'tafv':>5s} {'restv':>5s} "
          f"{'contacts':>9s}  all3  BDBV-diff")
    for r in rows:
        print(f"{r['ebov_full_pos']:>6d} {r['EBOV']:>5s} {r['BDBV']:>5s} {r['SUDV']:>5s} "
              f"{r['TAFV']:>5s} {r['RESTV']:>5s} {r['n_contact_atoms']:>9d}  "
              f"{str(r['in_all_copies']):>4s}  {'***' if r['bdbv_differs'] else ''}")
    print()
    print(f"interface residues: {len(rows)}; BDBV differs at {len(diffs)} of them:",
          [(r["ebov_full_pos"], r["EBOV"], r["BDBV"]) for r in diffs])


if __name__ == "__main__":
    main()
