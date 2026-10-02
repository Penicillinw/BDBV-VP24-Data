"""G6c: does baseline receptor/ISG expression proxy interferon responsiveness in
PRIMARY human immune cells? Primary test on GSE327707.

WHY THIS DATASET
  GSE327707 profiles FACS-sorted immune populations from 6 healthy donors after
  48 h of in-vitro stimulation with four stimuli (BCS, CytoStim, LPS, IFN-alpha)
  plus unstimulated baseline controls. Baseline and response therefore come from
  the same donor, the same cell type and the same experiment - the design the
  question actually needs, in the cell types the manuscript ranks (CD14 myeloid
  plus lymphoid comparators).

FROZEN BEFORE RUNNING (recorded here, not adjusted after seeing values)
  * Module definitions are copied verbatim from the cell-line arm
    (tools/g6_reliability_20260926.py): 14-gene ISG module, 4-gene receptor
    module, 6-gene housekeeping control, and the manuscript's own components
    (IFNAR1+IFNAR2; IFNLR1+IL10RB; ISG15+MX1) and delta_like composite.
  * Primary endpoint: Spearman rho between BASELINE ISG module score and the
    IFN-alpha induction magnitude (IFN-alpha minus baseline) of the same module,
    computed WITHIN each sorted cell type (n = 6 donors).
  * Primary pooled statistic: rank within cell type, then pool the ranks
    (n = 30 donor x cell-type units), so a between-cell-type mixture cannot create
    the correlation. The naive pooled value is reported as sensitivity only.
  * Secondary: baseline type-I receptor module, the manuscript's 2-gene priming
    term, the delta_like composite, and the cross-cell-type version (n = 5 cell
    types) that mirrors how the manuscript uses delta.
  * Negative controls: housekeeping module, and a 1000-draw random 14-gene-set
    null for |rho|.
  * Bimodality screen is run on the induction axis first, because a bimodal
    responder split is what invalidated the pooled statistic in the cell-line arm.

Outputs (analysis/g6c_primary_immune_20260926/gse327707/):
  sample_scores.tsv, g6c_gse327707_results.json
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import random
import re
import sys
from collections import defaultdict
from urllib.request import Request, urlopen

ROOT = r"G:\本迪布焦研究"
BASE = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926")
RAW = os.path.join(BASE, "raw")
OUT = os.path.join(BASE, "gse327707")

ISG = ["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3",
       "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6"]
RECEPTORS = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB"]
HK = ["ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "HPRT1"]
GENES = list(dict.fromkeys(ISG + RECEPTORS + HK))

CELL_TYPES = ["Whole PBMC", "CD14", "CD4", "CD8", "B cell"]
TREATMENTS = ["Baseline", "IFNa", "LPS", "Cytostim", "HEK BCS (CD79B)"]


def fetch_ensembl_symbols() -> dict[str, str]:
    """symbol -> Ensembl gene id, cached as raw JSON (project convention)."""
    cache = os.path.join(RAW, "ensembl_symbols_20260926.json")
    if os.path.exists(cache):
        return json.load(open(cache, encoding="utf-8"))
    payload = json.dumps({"symbols": GENES}).encode()
    req = Request(
        "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "CodexResearch/1.0",
        },
    )
    try:
        data = json.loads(urlopen(req, timeout=120).read().decode())  # noqa: S310
    except Exception as exc:  # noqa: BLE001 - fall back to one-by-one lookups
        print(f"ensembl POST failed ({exc}); falling back to per-symbol GET")
        data = {}
        for gene in GENES:
            try:
                url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json"
                rec = json.loads(
                    urlopen(Request(url, headers={"User-Agent": "CodexResearch/1.0"}), timeout=60).read().decode()
                )
                data[gene] = {"id": rec["id"]}
            except Exception as exc2:  # noqa: BLE001
                print(f"  {gene}: {exc2}")
    mapping = {sym: rec["id"].split(".")[0] for sym, rec in data.items() if rec and rec.get("id")}
    json.dump(mapping, open(cache, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return mapping


def load_design() -> dict[str, dict[str, str]]:
    """column name (library name) -> {donor, cell_type, treatment}."""
    text = open(os.path.join(RAW, "GSE327707_gsm_brief.txt"), encoding="utf-8", errors="replace").read()
    design: dict[str, dict[str, str]] = {}
    for block in re.split(r"(?m)^\^SAMPLE = ", text)[1:]:
        title = re.search(r"(?m)^!Sample_title = (.*)$", block)
        desc = re.search(r"(?m)^!Sample_description = Library name: (.*)$", block)
        if not title or not desc:
            continue
        parts = [p.strip() for p in title.group(1).split(",")]
        donor, cell_type, treatment = parts[0], parts[1], parts[2]
        design[desc.group(1).strip()] = {
            "gsm": block.split("\n", 1)[0].strip(),
            "donor": donor,
            "cell_type": cell_type,
            "treatment": treatment,
        }
    return design


def load_counts() -> tuple[list[str], dict[str, list[float]]]:
    path = os.path.join(RAW, "GSE327707_raw_counts.csv.gz")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)[1:]
        rows: dict[str, list[float]] = {}
        for row in reader:
            if not row:
                continue
            rows[row[0].split(".")[0]] = [float(x) for x in row[1:]]
    return header, rows


def log2_cpm(rows: dict[str, list[float]], n_cols: int) -> dict[str, list[float]]:
    totals = [sum(v[i] for v in rows.values()) for i in range(n_cols)]
    out: dict[str, list[float]] = {}
    for gene, vec in rows.items():
        out[gene] = [
            math.log2(v / totals[i] * 1e6 + 1) if totals[i] else 0.0 for i, v in enumerate(vec)
        ]
    return out


def spearman(a: list[float], b: list[float]) -> float:
    def rank(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    da = math.sqrt(sum((v - ma) ** 2 for v in ra))
    db = math.sqrt(sum((v - mb) ** 2 for v in rb))
    return num / (da * db) if da and db else float("nan")


def permutation_p(a: list[float], b: list[float], observed: float, draws: int = 20000) -> float:
    rng = random.Random(20260926)
    x = list(a)
    hits = 0
    for _ in range(draws):
        rng.shuffle(x)
        if abs(spearman(x, b)) >= abs(observed) - 1e-12:
            hits += 1
    return (hits + 1) / (draws + 1)


def rank_within(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    r = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(OUT, exist_ok=True)

    symbols = fetch_ensembl_symbols()
    design = load_design()
    columns, counts = load_counts()
    expr = log2_cpm(counts, len(columns))

    resolved: dict[str, str] = {}
    missing: list[str] = []
    for gene in GENES:
        eid = symbols.get(gene)
        if eid and eid in expr:
            resolved[gene] = eid
        else:
            missing.append(gene)
    print(f"module genes resolved: {len(resolved)}/{len(GENES)}")
    if missing:
        print(f"  missing from matrix/mapping: {', '.join(missing)}")

    # sample -> design row, and guard against silent mismatches
    units = []
    for i, col in enumerate(columns):
        meta = design.get(col)
        if not meta:
            raise SystemExit(f"no design metadata for count-matrix column {col}")
        units.append({"index": i, "column": col, **meta})
    print(f"samples: {len(units)}  donors: {len({u['donor'] for u in units})}  "
          f"cell types: {sorted({u['cell_type'] for u in units})}  "
          f"treatments: {sorted({u['treatment'] for u in units})}")

    def module_score(genes: list[str], idx: int) -> float:
        vals = [expr[resolved[g]][idx] for g in genes if g in resolved]
        return sum(vals) / len(vals) if vals else float("nan")

    # per-sample scores
    for u in units:
        u["isg"] = module_score(ISG, u["index"])
        u["isg_2gene"] = module_score(["ISG15", "MX1"], u["index"])
        u["typeI"] = module_score(["IFNAR1", "IFNAR2"], u["index"])
        u["typeIII"] = module_score(["IFNLR1", "IL10RB"], u["index"])
        u["hk"] = module_score(HK, u["index"])
        u["delta_like"] = (u["typeI"] + u["typeIII"] + u["isg_2gene"]) / 3

    def pick(cell_type: str, treatment: str) -> list[dict]:
        return [u for u in units if u["cell_type"] == cell_type and u["treatment"] == treatment]

    # ---- per cell type: is the induction phenotype graded or binary? --------
    bimodality = {}
    for ct in CELL_TYPES:
        base = {u["donor"]: u for u in pick(ct, "Baseline")}
        ifn = {u["donor"]: u for u in pick(ct, "IFNa")}
        donors = sorted(set(base) & set(ifn))
        ind = [ifn[d]["isg"] - base[d]["isg"] for d in donors]
        srt = sorted(ind)
        gap = srt[-1] - srt[0]
        max_gap = max((srt[i + 1] - srt[i] for i in range(len(srt) - 1)), default=0.0)
        bimodality[ct] = {
            "n_donors": len(donors),
            "induced_log_units": [round(v, 3) for v in ind],
            "min": round(srt[0], 3),
            "max": round(srt[-1], 3),
            "largest_single_gap": round(max_gap, 3),
            "largest_gap_fraction_of_range": round(max_gap / gap, 3) if gap else None,
            "looks_binary": bool(gap and max_gap / gap > 0.5),
        }

    # ---- primary: within cell type, baseline vs induction -------------------
    per_cell_type = {}
    for ct in CELL_TYPES:
        base = {u["donor"]: u for u in pick(ct, "Baseline")}
        ifn = {u["donor"]: u for u in pick(ct, "IFNa")}
        donors = sorted(set(base) & set(ifn))
        b_isg = [base[d]["isg"] for d in donors]
        i_isg = [ifn[d]["isg"] - base[d]["isg"] for d in donors]
        b_2g = [base[d]["isg_2gene"] for d in donors]
        i_2g = [ifn[d]["isg_2gene"] - base[d]["isg_2gene"] for d in donors]
        per_cell_type[ct] = {
            "donors": donors,
            "baseline_isg": [round(v, 3) for v in b_isg],
            "induction_isg": [round(v, 3) for v in i_isg],
            "induced_mean": round(sum(i_isg) / len(i_isg), 3),
            "rho_baseline_isg_vs_induction": round(spearman(b_isg, i_isg), 3),
            "p_permutation": round(permutation_p(b_isg, i_isg, spearman(b_isg, i_isg)), 4),
            "rho_baseline_2gene_vs_own_induction": round(spearman(b_2g, i_2g), 3),
            "rho_baseline_typeI_vs_isg_induction":
                round(spearman([base[d]["typeI"] for d in donors], i_isg), 3),
            "rho_baseline_typeIII_vs_isg_induction":
                round(spearman([base[d]["typeIII"] for d in donors], i_isg), 3),
            "rho_delta_like_vs_isg_induction":
                round(spearman([base[d]["delta_like"] for d in donors], i_isg), 3),
            "rho_housekeeping_vs_isg_induction":
                round(spearman([base[d]["hk"] for d in donors], i_isg), 3),
        }

    # ---- pooled: rank within cell type (mixture-free) -----------------------
    rank_b: list[float] = []
    rank_i: list[float] = []
    raw_b: list[float] = []
    raw_i: list[float] = []
    for ct in CELL_TYPES:
        base = {u["donor"]: u for u in pick(ct, "Baseline")}
        ifn = {u["donor"]: u for u in pick(ct, "IFNa")}
        donors = sorted(set(base) & set(ifn))
        bs = [base[d]["isg"] for d in donors]
        ii = [ifn[d]["isg"] - base[d]["isg"] for d in donors]
        rank_b += rank_within(bs)
        rank_i += rank_within(ii)
        raw_b += bs
        raw_i += ii
    pooled_within = spearman(rank_b, rank_i)
    pooled_within_p = permutation_p(rank_b, rank_i, pooled_within)
    rows = [
        {"resource": "GSE327707", "donor": u["donor"], "cell_type": u["cell_type"],
         "treatment": u["treatment"], "isg": round(u["isg"], 4), "isg_2gene": round(u["isg_2gene"], 4),
         "typeI": round(u["typeI"], 4), "typeIII": round(u["typeIII"], 4),
         "housekeeping": round(u["hk"], 4), "delta_like": round(u["delta_like"], 4)}
        for u in units
    ]

    # ---- cross-cell-type axis (mirrors how the manuscript uses delta) -------
    cross = []
    for ct in CELL_TYPES:
        base = [u["isg"] for u in pick(ct, "Baseline")]
        ind = [
            u["isg"] - b["isg"]
            for u, b in zip(sorted(pick(ct, "IFNa"), key=lambda x: x["donor"]),
                            sorted(pick(ct, "Baseline"), key=lambda x: x["donor"]))
        ]
        cross.append({
            "cell_type": ct,
            "mean_baseline_isg": round(sum(base) / len(base), 3),
            "mean_induction_isg": round(sum(ind) / len(ind), 3),
        })
    cross_rho = round(
        spearman([c["mean_baseline_isg"] for c in cross], [c["mean_induction_isg"] for c in cross]), 3
    )

    # ---- random 14-gene-set null for |rho| ---------------------------------
    # background is restricted to genes with detectable expression, otherwise the
    # null is dominated by constant zero vectors and |rho| is undefined
    n_cols = len(columns)
    background = [
        g for g, vec in expr.items() if sum(vec) / n_cols >= 3.0
    ]
    rng = random.Random(20260926)
    null_rhos = []
    for _ in range(1000):
        sample_genes = rng.sample(background, 14)
        rb: list[float] = []
        ri: list[float] = []
        for ct in CELL_TYPES:
            base = {u["donor"]: u for u in pick(ct, "Baseline")}
            ifn = {u["donor"]: u for u in pick(ct, "IFNa")}
            donors = sorted(set(base) & set(ifn))
            vb = [sum(expr[g][base[d]["index"]] for g in sample_genes) / 14 for d in donors]
            vi = [sum(expr[g][ifn[d]["index"]] - expr[g][base[d]["index"]] for g in sample_genes) / 14
                  for d in donors]
            rb += rank_within(vb)
            ri += rank_within(vi)
        null_rhos.append(abs(spearman(rb, ri)))
    null_rhos.sort()
    observed_abs = abs(pooled_within)
    percentile = 100 * sum(1 for v in null_rhos if v <= observed_abs) / len(null_rhos)

    result = {
        "dataset": "GSE327707",
        "design": "6 healthy donors; FACS-sorted Whole PBMC / CD14 / CD4 / CD8 / B cell; "
                  "48 h BCS, CytoStim, LPS, IFN-alpha, and unstimulated baseline; bulk RNA-seq",
        "modules_frozen": {"ISG": ISG, "receptors": RECEPTORS, "housekeeping": HK},
        "genes_missing": missing,
        "bimodality_screen": bimodality,
        "per_cell_type": per_cell_type,
        "pooled_within_cell_type": {
            "n_units": len(rank_b),
            "spearman_baseline_rank_vs_induction_rank": round(pooled_within, 3),
            "p_permutation": round(pooled_within_p, 4),
        },
        "pooled_naive": {
            "spearman_raw_baseline_vs_induction": round(spearman(raw_b, raw_i), 3),
        },
        "cross_cell_type": {"rows": cross, "spearman_mean_baseline_vs_mean_induction": cross_rho},
        "random_gene_set_null": {
            "background_genes": len(background),
            "draws": len(null_rhos),
            "median_abs_rho": round(null_rhos[len(null_rhos) // 2], 3),
            "p95_abs_rho": round(null_rhos[int(0.95 * len(null_rhos))], 3),
            "observed_abs_rho_percentile": round(percentile, 1),
        },
    }
    json.dump(result, open(os.path.join(OUT, "g6c_gse327707_results.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open(os.path.join(OUT, "sample_scores.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps({k: result[k] for k in
                      ("bimodality_screen", "per_cell_type", "pooled_within_cell_type",
                       "pooled_naive", "cross_cell_type", "random_gene_set_null")},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
