"""Step 4: genome-wide response size per interferon arm (independent check).

Frozen criterion: a gene counts as differentially expressed if the paired
 t-test on log2(TPM+1) across the 5 replicate-matched control/treated pairs of
 the same (time, route) cell reaches BH-FDR < 0.05 **and** |mean log2 fold
 change| >= 1. With n = 5 pairs, this is a permissive but symmetric rule; it is
 applied identically to every arm, so the comparison between arms is fair.

Purpose: the source study states that IFN-beta1 activates thousands of genes
 while type III IFNs activate hundreds, and that IFN-lambda4 is inactive. If our
pipeline reproduces that ordering from the deposited TPM matrix alone, the
pipeline is validated against the authors' own claim.

Outputs: deg_counts.tsv, g342_deg.json (this directory only).
"""

import gzip
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
MATRIX = HERE / "raw" / "GSE342661_counts_tpm.matrix.gz"


def main() -> None:
    design = pd.read_csv(HERE / "design.tsv", sep="\t")
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        tpm = pd.read_csv(fh, sep="\t", quotechar='"', index_col=0, low_memory=False)
    tpm = tpm.astype("float64")
    log2tpm = np.log2(tpm + 1.0)
    log2tpm.columns = [c.strip() for c in log2tpm.columns]
    design = design.copy()
    design["column"] = design["column"].str.strip()
    missing = [c for c in design["column"] if c not in log2tpm.columns]
    assert not missing, missing[:5]

    rows = []
    for route in ["apical-basal", "apical-only"]:
        for time in sorted(design["time"].unique()):
            cell = design[(design["time"] == time) & (design["route"] == route)]
            ctrl_cols = cell[cell["treatment"] == "Ctrl"].sort_values("rep")["column"]
            for trt in ["IFNB1", "IFNL1", "IFNL2", "IFNL3", "IFNL4"]:
                trt_cols = cell[cell["treatment"] == trt].sort_values("rep")["column"]
                if len(trt_cols) != len(ctrl_cols) or len(trt_cols) == 0:
                    continue
                a = log2tpm[list(trt_cols)].values
                b = log2tpm[list(ctrl_cols)].values
                diff = a - b
                mean_diff = diff.mean(axis=1)
                t = stats.ttest_rel(a, b, axis=1)
                p = np.asarray(t.pvalue, dtype=float)
                p = np.where(np.isnan(p), 1.0, p)
                order = np.argsort(p)
                ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
                ranked = np.minimum.accumulate(ranked[::-1])[::-1]
                q = np.empty_like(ranked)
                q[order] = np.clip(ranked, 0, 1)
                sig = (q < 0.05) & (np.abs(mean_diff) >= 1.0)
                rows.append(
                    {
                        "route": route,
                        "time": time,
                        "treatment": trt,
                        "n_pairs": int(len(trt_cols)),
                        "n_tested": int(len(p)),
                        "n_sig_up": int((sig & (mean_diff > 0)).sum()),
                        "n_sig_down": int((sig & (mean_diff < 0)).sum()),
                        "max_abs_mean_log2_fc": float(np.abs(mean_diff).max()),
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "deg_counts.tsv", sep="\t", index=False)
    summary = {
        "rule": "paired t on log2(TPM+1), BH-FDR < 0.05 and |mean log2FC| >= 1",
        "rows": out.to_dict("records"),
    }
    (HERE / "g342_deg.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    pd.set_option("display.width", 200)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
