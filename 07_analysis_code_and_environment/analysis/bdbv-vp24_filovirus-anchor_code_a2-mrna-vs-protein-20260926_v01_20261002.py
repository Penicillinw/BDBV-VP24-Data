"""A2: do the genes that make up Δ agree between the mRNA layer and the protein layer
in the same reference resource?

Δ is built from HPA single-cell-type mRNA (nCPM) for IFNAR1/2, IFNLR1, IL10RB, ISG15,
MX1 and the importin-alpha genes. HPA also annotates protein-level tissue and cell-type
distribution from immunohistochemistry. G7 produced one contradiction (KPNA5 mRNA is
low in enterocytes while the protein annotation calls the intestine enriched), so this
script audits every gene in the score the same way and records antibody reliability.
"""

from __future__ import annotations

import csv
import io
import json
import os
from urllib.request import Request, urlopen

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "a2_layer_consistency_20260926")
RAW = os.path.join(OUT, "raw")
HPA_TABLE = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")
GRADIENT = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925",
                        "hpa_ifn_landscape_wide.tsv")

GENES = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1",
         "KPNA1", "KPNA5", "KPNA6", "STAT1", "STAT2", "IRF9"]

# the cell types the manuscript's prediction actually names
PREDICTION_CELLS = [
    "enterocytes", "colonocytes", "urothelial cells", "kupffer cells", "hepatocytes",
    "monocytes", "macrophages", "microglia", "pdcs", "respiratory ciliated cells",
    "ocular epithelial cells", "retinal pigment epithelial cells", "sertoli cells",
    "leydig cells", "rod photoreceptor cells", "choroid plexus epithelial cells",
]


def fetch(url: str, dest: str) -> dict | None:
    if os.path.exists(dest) and os.path.getsize(dest) > 100:
        try:
            return json.load(open(dest, encoding="utf-8"))
        except Exception:  # noqa: BLE001
            pass
    for attempt in range(4):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=120) as r:  # noqa: S310
                data = r.read()
            with open(dest, "wb") as fh:
                fh.write(data)
            return json.loads(data.decode("utf-8", "replace"))
        except Exception:  # noqa: BLE001
            import time

            time.sleep(2 * (attempt + 1))
    return None


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    # ensembl ids + mRNA values
    ens: dict[str, str] = {}
    vals: dict[str, dict[str, float]] = {}
    with open(HPA_TABLE, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        gi, ni, ci, vi = (header.index("Gene"), header.index("Gene name"),
                          header.index("Cell type"), header.index("nCPM"))
        for row in rd:
            g = row[ni]
            if g in GENES:
                ens.setdefault(g, row[gi])
                try:
                    vals.setdefault(g, {})[row[ci]] = float(row[vi])
                except ValueError:
                    pass

    rows = []
    for g in GENES:
        eid = ens.get(g, "")
        js = fetch(f"https://www.proteinatlas.org/{eid}.json",
                   os.path.join(RAW, f"{g}_{eid}.json")) if eid else None
        rec: dict[str, object] = {"gene": g, "ensembl": eid}
        if js:
            rec.update({
                "rna_sc_specificity": js.get("RNA single cell type specificity"),
                "rna_sc_top": json.dumps(js.get("RNA single cell type specific nCPM"),
                                         ensure_ascii=False)[:160],
                "protein_celltype_specificity": js.get("Protein cell type specificity"),
                "protein_celltype_intensity":
                    json.dumps(js.get("Protein cell type specific Intensity"),
                               ensure_ascii=False)[:200],
                "protein_tissue_specificity": js.get("Protein tissue specificity"),
                "protein_tissue_intensity":
                    json.dumps(js.get("Protein tissue specific Intensity"),
                               ensure_ascii=False)[:200],
                "reliability_IH": js.get("Reliability (IH)"),
                "reliability_IF": js.get("Reliability (IF)"),
                "evidence": js.get("Evidence"),
            })
        for ct in PREDICTION_CELLS:
            rec[f"nCPM_{ct}"] = round(vals.get(g, {}).get(ct, float("nan")), 1) \
                if ct in vals.get(g, {}) else None
        rows.append(rec)

    with open(os.path.join(OUT, "layer_consistency.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"{'gene':8s} {'RNA-sc top':26s} {'protein tissue (intensity)':30s} "
          f"{'IH':10s} {'IF':10s}")
    for r in rows:
        print(f"{r['gene']:8s} {str(r.get('rna_sc_top'))[:25]:26s} "
              f"{str(r.get('protein_tissue_intensity'))[:29]:30s} "
              f"{str(r.get('reliability_IH'))[:9]:10s} {str(r.get('reliability_IF'))[:9]:10s}")

    print("\n-- mRNA in the cell types the prediction names --")
    for r in rows:
        got = {k[5:]: v for k, v in r.items() if k.startswith("nCPM_") and v is not None}
        top = sorted(got.items(), key=lambda kv: -kv[1])[:5]
        print(f"{r['gene']:8s} " + ", ".join(f"{k}={v}" for k, v in top))

    json.dump(rows, open(os.path.join(OUT, "layer_consistency.json"), "w",
                         encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\nwrote", os.path.join(OUT, "layer_consistency.tsv"))


if __name__ == "__main__":
    main()
