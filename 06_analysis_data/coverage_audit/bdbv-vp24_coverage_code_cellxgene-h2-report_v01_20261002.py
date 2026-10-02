"""H2 line: print the candidate atlas table (cell-type hits) sorted by size."""
from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader((HERE / "cellxgene_dataset_details.tsv").open(encoding="utf-8"),
                           delimiter="\t"))


def size_gb(r):
    try:
        return float(r.get("h5ad_size") or 0) / 1e9
    except (TypeError, ValueError):
        return 0.0


hits = [r for r in rows if r["celltype_hits"]]
hits.sort(key=lambda r: (r["celltype_hits"], size_gb(r)))
print(f"{'celltype_hits':36s} {'cells':>9s} {'size_GB':>8s}  primary        ds_id      title")
for r in hits:
    print(f"{r['celltype_hits'][:34]:36s} {r['cell_count']:>9s} {size_gb(r):8.2f} "
          f" {r['is_primary_data'][:14]:14s} {r['dataset_id'][:8]} {r['dataset_title'][:56]}")

print("\n--- spleen by tissue label (top 20 by cells) ---")
sp = [r for r in rows if "spleen" in (r.get("tissue") or "").lower()]
sp.sort(key=lambda r: -(int(r["cell_count"] or 0)))
for r in sp[:20]:
    print(f"{r['cell_count']:>9s} {size_gb(r):8.2f}GB prim={r['is_primary_data'][:12]:12s} "
          f"{r['dataset_id'][:8]} {r['dataset_title'][:50]:52s} | nCT={r['n_cell_types']:>3s} "
          f"| {r['cell_type_labels'][:90]}")
