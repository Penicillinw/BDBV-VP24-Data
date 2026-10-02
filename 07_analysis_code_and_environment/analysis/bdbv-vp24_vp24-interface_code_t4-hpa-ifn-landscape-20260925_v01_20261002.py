"""Extract the IFN-system landscape from the LOCAL HPA single-cell-type table.

Input : data/hpa/rna_single_cell_type.tsv  (Gene, Gene name, Cell type, nCPM)
Output: analysis/t4_ifn_landscape_20260925/hpa_ifn_landscape_{long,wide}.tsv
        analysis/t4_ifn_landscape_20260925/hpa_ifn_summary.json

The HPA single-cell table is derived from the Human Protein Atlas
(Consensus dataset, scRNA-seq of human tissues).  Values are nCPM
(normalised counts per million).
"""

import csv
import json
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")
OUT = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")

GENES = [
    # type III (lambda) system
    "IFNLR1", "IL10RB", "IFNL1", "IFNL2", "IFNL3", "IFNL4",
    # type I system
    "IFNAR1", "IFNAR2", "IFNB1", "IFNA1", "IFNA2",
    # signalling
    "STAT1", "STAT2", "JAK1", "JAK2", "IRF3", "IRF7",
    # ISGs
    "ISG15", "MX1", "MX2", "OAS1", "IFITM1", "IFITM2", "IFITM3", "BST2",
    # VP24 host targets
    "KPNA1", "KPNA5", "KPNA6",
    # reference
    "NPC1",
]


def main():
    os.makedirs(OUT, exist_ok=True)
    want = set(GENES)
    values = defaultdict(dict)  # gene -> cell type -> nCPM
    cell_types = set()
    n_rows = 0
    with open(SRC, encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        gi, ci, vi = idx["Gene name"], idx["Cell type"], idx["nCPM"]
        for row in reader:
            n_rows += 1
            gene = row[gi]
            if gene not in want:
                continue
            ct = row[ci]
            try:
                value = float(row[vi])
            except ValueError:
                continue
            values[gene][ct] = value
            cell_types.add(ct)

    cell_types = sorted(cell_types)
    long_path = os.path.join(OUT, "hpa_ifn_landscape_long.tsv")
    with open(long_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene", "cell_type", "nCPM"])
        for gene in GENES:
            for ct in sorted(values.get(gene, {})):
                writer.writerow([gene, ct, values[gene][ct]])

    wide_path = os.path.join(OUT, "hpa_ifn_landscape_wide.tsv")
    present = [g for g in GENES if g in values and values[g]]
    with open(wide_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["cell_type"] + present)
        for ct in cell_types:
            writer.writerow([ct] + [values[g].get(ct, "") for g in present])

    summary = {
        "source_file": os.path.relpath(SRC, ROOT),
        "source_rows_scanned": n_rows,
        "unit": "nCPM (HPA single-cell-type RNA, Consensus)",
        "genes_requested": GENES,
        "genes_present": present,
        "genes_missing": [g for g in GENES if g not in present],
        "n_cell_types": len(cell_types),
        "cell_types": cell_types,
        "top_cell_types_by_gene": {},
    }
    for g in present:
        top = sorted(values[g].items(), key=lambda kv: -kv[1])[:15]
        summary["top_cell_types_by_gene"][g] = [
            {"cell_type": ct, "nCPM": round(v, 1)} for ct, v in top
        ]
    with open(os.path.join(OUT, "hpa_ifn_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    print(f"scanned {n_rows} rows; cell types with any target gene: {len(cell_types)}")
    print(f"genes present: {len(present)}/{len(GENES)}; missing: {summary['genes_missing']}")
    print()
    for g in ["IFNLR1", "IL10RB", "IFNAR1", "IFNAR2", "STAT1", "ISG15", "MX1", "IFITM3"]:
        if g not in summary["top_cell_types_by_gene"]:
            continue
        tops = summary["top_cell_types_by_gene"][g][:6]
        print(f"{g:8s} top: " + ", ".join(f"{t['cell_type']}={t['nCPM']}" for t in tops))


if __name__ == "__main__":
    main()
