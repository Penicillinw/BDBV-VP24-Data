"""B1 positive control - step 1: recon of the two experimental GPcl-NPC1 complexes.

Inputs (downloaded from RCSB):
  5F1B  EBOV GP1+GP2 + human NPC1 domain C, X-ray 2.3 A
  9DZ2  SUDV GP + human NPC1-C, cryo-EM 3.31 A

Outputs:
  analysis/b1_positive_control_20260924/recon_chains.tsv
  analysis/b1_positive_control_20260924/recon_contacts.tsv

Only reading is done here; no structure is modified.
"""

from __future__ import annotations

import json
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

ROOT = Path(r"G:\本迪布焦研究")
RAW = ROOT / "analysis" / "b1_positive_control_20260924" / "raw"
OUT = ROOT / "analysis" / "b1_positive_control_20260924"
OUT.mkdir(parents=True, exist_ok=True)

CONTACT_CUTOFF = 4.0


def load(stem: str):
    cif = RAW / f"{stem}.cif"
    pdb = RAW / f"{stem}.pdb"
    if cif.exists():
        return MMCIFParser(QUIET=True).get_structure(stem, str(cif))
    return PDBParser(QUIET=True).get_structure(stem, str(pdb))


def chain_seq(chain) -> str:
    from Bio.Data.PDBData import protein_letters_3to1_extended

    out = []
    for res in chain:
        if is_aa(res, standard=False):
            out.append(protein_letters_3to1_extended.get(res.get_resname(), "X"))
    return "".join(out)


def main() -> None:
    chain_rows = []
    contact_rows = []
    summary = {}

    for stem in ("5F1B", "9DZ2"):
        structure = load(stem)
        model = next(iter(structure))
        chains = [c for c in model if c.id.strip()]
        summary[stem] = {"n_chains": len(chains), "chains": [c.id for c in chains]}

        for chain in chains:
            seq = chain_seq(chain)
            heavy = [a for a in chain.get_atoms() if a.element != "H"]
            chain_rows.append(
                {
                    "structure": stem,
                    "chain": chain.id,
                    "n_residues": len(seq),
                    "n_heavy_atoms": len(heavy),
                    "first_res": chain.child_list[0].id[1] if len(chain) else "",
                    "last_res": chain.child_list[-1].id[1] if len(chain) else "",
                    "sequence": seq,
                }
            )

        # pairwise interface contacts between protein chains
        for i, ci in enumerate(chains):
            atoms_i = [a for a in ci.get_atoms() if a.element != "H"]
            if not atoms_i:
                continue
            ns = NeighborSearch(atoms_i)
            for cj in chains[i + 1 :]:
                atoms_j = [a for a in cj.get_atoms() if a.element != "H"]
                if not atoms_j:
                    continue
                pairs = set()
                for atom in atoms_j:
                    for hit in ns.search(atom.coord, CONTACT_CUTOFF):
                        pairs.add((hit.get_parent().id[1], atom.get_parent().id[1]))
                if not pairs:
                    continue
                ri = sorted({p[0] for p in pairs})
                rj = sorted({p[1] for p in pairs})
                contact_rows.append(
                    {
                        "structure": stem,
                        "chain_a": ci.id,
                        "chain_b": cj.id,
                        "n_atom_pairs": len(pairs),
                        "n_res_a": len(ri),
                        "n_res_b": len(rj),
                        "res_a": ",".join(str(x) for x in ri),
                        "res_b": ",".join(str(x) for x in rj),
                    }
                )

    with (OUT / "recon_chains.tsv").open("w", encoding="utf-8", newline="") as fh:
        import csv

        w = csv.DictWriter(fh, fieldnames=list(chain_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(chain_rows)

    with (OUT / "recon_contacts.tsv").open("w", encoding="utf-8", newline="") as fh:
        import csv

        w = csv.DictWriter(fh, fieldnames=list(contact_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(contact_rows)

    (OUT / "recon_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
