"""G7 orthogonal cross-check: HPA PROTEIN layer for the KPNA axis genes.

Why: the transcript-level G7 result puts gut epithelium (enterocytes,
colonocytes, ...) into Q2 ("high Delta, low KPNA").  Before that can be read as
a biological statement, it has to survive a non-transcript layer.  HPA provides
per-gene protein annotations (immunohistochemistry reliability, tissue/cell-type
intensity) from raw JSON already fetched into
analysis/g7_kpna_axis_20260926/hpa_raw/.

POST-HOC: not part of the G7 pre-registration.
Outputs:
    analysis/g7_kpna_axis_20260926/kpna_protein_layer.tsv
    analysis/g7_kpna_axis_20260926/kpna_protein_layer.json
"""

import csv
import glob
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "analysis", "g7_kpna_axis_20260926", "hpa_raw")
OUT = os.path.join(ROOT, "analysis", "g7_kpna_axis_20260926")

CELLS = ["Enterocytes", "Colonocytes", "Goblet cells", "Paneth cells",
         "Hepatocytes", "Kupffer cells", "Macrophages", "Neutrophils",
         "Monocytes", "Late primary spermatocytes", "Rod photoreceptor cells",
         "Sertoli cells", "Leydig cells", "Myonuclei", "Adrenal cortex cells"]


def as_map(value):
    return value if isinstance(value, dict) else {}


def main():
    rows, payload = [], {}
    for path in sorted(glob.glob(os.path.join(RAW, "*.json"))):
        d = json.load(open(path, encoding="utf-8"))
        gene = d["Gene"]
        sc = as_map(d.get("RNA single cell type specific nCPM"))
        tissue = as_map(d.get("RNA tissue specific nTPM"))
        rec = {
            "gene": gene,
            "ensembl": d.get("Ensembl"),
            "uniprot": ";".join(d.get("Uniprot") or []),
            "evidence": d.get("Evidence"),
            "rna_tissue_specificity": d.get("RNA tissue specificity"),
            "rna_tissue_distribution": d.get("RNA tissue distribution"),
            "rna_tissue_top": ";".join(tissue),
            "rna_sc_specificity": d.get("RNA single cell type specificity"),
            "rna_sc_distribution": d.get("RNA single cell type distribution"),
            "rna_sc_top": ";".join(sc),
            "protein_tissue_specificity": d.get("Protein tissue specificity"),
            "protein_tissue_distribution": d.get("Protein tissue distribution"),
            "protein_tissue_intensity": ";".join(
                as_map(d.get("Protein tissue specific Intensity"))),
            "protein_celltype_specificity": d.get("Protein cell type specificity"),
            "protein_celltype_intensity": ";".join(
                as_map(d.get("Protein cell type specific Intensity"))),
            "reliability_IH": d.get("Reliability (IH)"),
            "reliability_IF": d.get("Reliability (IF)"),
            "antibody": ";".join(d.get("Antibody") or []),
            "subcellular": ";".join(d.get("Subcellular location") or []),
        }
        rows.append(rec)
        payload[gene] = {"annotations": rec,
                         "rna_single_cell_nCPM_selected": {c: sc.get(c) for c in CELLS}}

    fields = list(rows[0].keys())
    with open(os.path.join(OUT, "kpna_protein_layer.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["gene"]))
    with open(os.path.join(OUT, "kpna_protein_layer.json"), "w", encoding="utf-8") as fh:
        json.dump({"note": "POST-HOC cross-check; HPA per-gene JSON, fetched 2026-09-26",
                   "source_dir": "analysis/g7_kpna_axis_20260926/hpa_raw",
                   "genes": payload}, fh, ensure_ascii=False, indent=1)

    hdr = ["gene", "protein_tissue_specificity", "protein_tissue_intensity",
           "protein_celltype_specificity", "reliability_IH", "rna_tissue_top",
           "rna_sc_top"]
    print(" | ".join(hdr))
    for r in sorted(rows, key=lambda r: r["gene"]):
        print(" | ".join(str(r[h]) for h in hdr))
    print("")
    for g in sorted(payload):
        sel = payload[g]["rna_single_cell_nCPM_selected"]
        print("%-7s " % g + "  ".join("%s=%s" % (k, v) for k, v in sel.items()
                                      if v is not None))
    print("")
    print("wrote kpna_protein_layer.tsv / .json")


if __name__ == "__main__":
    main()
