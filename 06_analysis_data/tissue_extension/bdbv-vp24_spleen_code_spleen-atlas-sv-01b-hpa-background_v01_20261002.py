"""SRB-VERIFY step 1b: what exactly is the HPA gene background used by the
expression-matched null?

SRB's `gene_background()` (srb_06_nulls.py) kept every gene with *exactly 154*
rows in data/hpa/rna_single_cell_type.tsv and reported 18,394 such genes.  This
probe checks whether "154 rows" really coincides with the frozen 154-type
reference panel, because the null is only interpretable if it does.

Read-only; writes raw/sv_hpa_background_probe.json
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
RAW = os.path.join(ROOT, "analysis", "srb_verify_20260926", "raw")
T4 = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")
HPA_SC = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")


def main() -> None:
    with open(os.path.join(T4, "hpa_ifn_landscape_wide.tsv"), encoding="utf-8") as fh:
        frozen_types = [r["cell_type"] for r in csv.DictReader(fh, delimiter="\t")]
    frozen_set = set(frozen_types)

    per_gene = defaultdict(list)
    with open(HPA_SC, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                per_gene[r["Gene name"]].append((r["Cell type"], math.log10(float(r["nCPM"]) + 1.0)))
            except (ValueError, TypeError, KeyError):
                continue

    all_types = sorted({ct for v in per_gene.values() for ct, _ in v})
    widths = Counter(len(v) for v in per_gene.values())
    by_width = defaultdict(list)
    for g, v in per_gene.items():
        by_width[len(v)].append(g)

    # how well do the 154-row genes correspond to the frozen panel?
    exact154 = [g for g, v in per_gene.items() if len(v) == 154]
    set154 = {g for g, v in per_gene.items()
              if {ct for ct, _ in v} == frozen_set}
    covers154 = {g for g, v in per_gene.items()
                 if frozen_set <= {ct for ct, _ in v}}
    rows154 = set(by_width[154])

    out = {
        "n_genes_total": len(per_gene),
        "n_cell_types_in_file": len(all_types),
        "frozen_panel_n": len(frozen_types),
        "frozen_panel_inside_file": len(frozen_set & set(all_types)),
        "n_genes_with_exactly_154_rows": len(exact154),
        "n_genes_whose_rowset_equals_frozen_panel": len(set154),
        "n_genes_covering_all_frozen_types": len(covers154),
        "is_exact154_same_as_rowset_equals": rows154 == set154,
        "width_histogram_top20": widths.most_common(20),
        "width_values": sorted(widths)[:12] + ["..."] + sorted(widths)[-12:],
    }
    with open(os.path.join(RAW, "sv_hpa_background_probe.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
