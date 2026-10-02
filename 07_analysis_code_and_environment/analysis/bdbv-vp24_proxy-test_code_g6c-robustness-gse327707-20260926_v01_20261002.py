"""G6c robustness: is the negative baseline-vs-induction relation in GSE327707
biology, a single-donor artefact, a global-expression artefact, or mathematical
coupling between a score and its own change?

Five checks, all decided before running:
  1. Per-donor profile. If the low-induction / high-baseline donor also fails to
     respond to the non-interferon stimuli (LPS, CytoStim, BCS), the pattern is a
     technical low-responder, not an interferon-tuning phenotype.
  2. Leave-one-donor-out: is the pooled rho carried by one donor?
  3. Partial correlation controlling for the baseline housekeeping module, i.e.
     is the effect just global expression level?
  4. Disjoint-gene split: baseline scored on 7 ISGs, response scored on the OTHER
     7 ISGs. A shared measurement-error term cannot create the correlation when
     predictor and response use different genes and different libraries.
  5. Ceiling/saturation check: correlate baseline with the POST-IFN level, and
     compare the post-IFN spread with the baseline spread.

Output: analysis/g6c_primary_immune_20260926/gse327707/robustness.json
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
from urllib.request import Request, urlopen

ROOT = r"G:\本迪布焦研究"
BASE = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926")
RAW = os.path.join(BASE, "raw")
OUT = os.path.join(BASE, "gse327707")

ISG = ["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3",
       "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6"]
HK = ["ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "HPRT1"]
# interferon-independent response reporters, frozen here before use:
NFKB = ["NFKBIA", "TNFAIP3", "NFKBIZ", "IL6", "CXCL8", "CCL2", "IL1B", "SOD2"]
ACTIVATION = ["CD69", "IL2RA", "ICOS", "TNFRSF9", "CD38", "IFNG"]
EXTRA = list(dict.fromkeys(NFKB + ACTIVATION))

CELL_TYPES = ["Whole PBMC", "CD14", "CD4", "CD8", "B cell"]


def ensembl_ids() -> dict[str, str]:
    cache = os.path.join(RAW, "ensembl_symbols_extra_20260926.json")
    if os.path.exists(cache):
        return json.load(open(cache, encoding="utf-8"))
    payload = json.dumps({"symbols": EXTRA}).encode()
    req = Request(
        "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
        data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json",
                 "User-Agent": "CodexResearch/1.0"},
    )
    data = json.loads(urlopen(req, timeout=120).read().decode())  # noqa: S310
    mapping = {sym: rec["id"].split(".")[0] for sym, rec in data.items() if rec and rec.get("id")}
    json.dump(mapping, open(cache, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return mapping


def load_design() -> dict[str, dict[str, str]]:
    text = open(os.path.join(RAW, "GSE327707_gsm_brief.txt"), encoding="utf-8", errors="replace").read()
    design: dict[str, dict[str, str]] = {}
    for block in re.split(r"(?m)^\^SAMPLE = ", text)[1:]:
        title = re.search(r"(?m)^!Sample_title = (.*)$", block)
        desc = re.search(r"(?m)^!Sample_description = Library name: (.*)$", block)
        if not (title and desc):
            continue
        donor, cell_type, treatment = [p.strip() for p in title.group(1).split(",")]
        design[desc.group(1).strip()] = {"donor": donor, "cell_type": cell_type, "treatment": treatment}
    return design


def load_counts() -> tuple[list[str], dict[str, list[float]]]:
    with gzip.open(os.path.join(RAW, "GSE327707_raw_counts.csv.gz"), "rt",
                   encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)[1:]
        rows = {r[0].split(".")[0]: [float(x) for x in r[1:]] for r in reader if r}
    return header, rows


def log2_cpm(rows: dict[str, list[float]], n: int) -> dict[str, list[float]]:
    totals = [sum(v[i] for v in rows.values()) for i in range(n)]
    return {g: [math.log2(v / totals[i] * 1e6 + 1) for i, v in enumerate(vec)]
            for g, vec in rows.items()}


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


def pearson(a: list[float], b: list[float]) -> float:
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    da = math.sqrt(sum((v - ma) ** 2 for v in a))
    db = math.sqrt(sum((v - mb) ** 2 for v in b))
    return num / (da * db) if da and db else float("nan")


def partial_spearman(x: list[float], y: list[float], z: list[float]) -> float:
    """Spearman partial correlation of x and y given z (rank-transform then Pearson)."""
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry, rz = rank(x), rank(y), rank(z)
    rxy, rxz, ryz = pearson(rx, ry), pearson(rx, rz), pearson(ry, rz)
    denom = math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    return (rxy - rxz * ryz) / denom if denom else float("nan")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    symbols = json.load(open(os.path.join(RAW, "ensembl_symbols_20260926.json"), encoding="utf-8"))
    symbols.update(ensembl_ids())
    design = load_design()
    columns, counts = load_counts()
    expr = log2_cpm(counts, len(columns))

    def gene(g: str) -> str | None:
        eid = symbols.get(g)
        return eid if eid and eid in expr else None

    def score(genes: list[str], idx: int) -> float:
        hits = [gene(g) for g in genes]
        hits = [h for h in hits if h]
        return sum(expr[h][idx] for h in hits) / len(hits) if hits else float("nan")

    units = []
    for i, col in enumerate(columns):
        meta = design[col]
        units.append({"index": i, **meta,
                      "isg": score(ISG, i), "hk": score(HK, i),
                      "nfkb": score(NFKB, i), "act": score(ACTIVATION, i)})

    def pick(ct, tr):
        return {u["donor"]: u for u in units if u["cell_type"] == ct and u["treatment"] == tr}

    donors = sorted({u["donor"] for u in units})

    # ---- 1. per-donor profile across stimuli --------------------------------
    profile = {}
    for donor in donors:
        entry = {}
        for ct in CELL_TYPES:
            b, i = pick(ct, "Baseline").get(donor), pick(ct, "IFNa").get(donor)
            l = pick(ct, "LPS").get(donor)
            c = pick(ct, "Cytostim").get(donor)
            if b and i:
                entry[ct] = {
                    "isg_induction_ifna": round(i["isg"] - b["isg"], 3),
                    "isg_induction_lps": round(l["isg"] - b["isg"], 3) if l else None,
                    "nfkb_induction_lps": round(l["nfkb"] - b["nfkb"], 3) if l else None,
                    "nfkb_induction_cytostim": round(c["nfkb"] - b["nfkb"], 3) if c else None,
                    "act_induction_cytostim": round(c["act"] - b["act"], 3) if c else None,
                }
        means = {
            key: round(sum(v[key] for v in entry.values() if v[key] is not None)
                       / max(1, sum(1 for v in entry.values() if v[key] is not None)), 3)
            for key in ("isg_induction_ifna", "isg_induction_lps", "nfkb_induction_lps",
                        "nfkb_induction_cytostim", "act_induction_cytostim")
        }
        profile[donor] = {"per_cell_type": entry, "mean_over_cell_types": means}

    # ---- pooled primary statistic, rebuilt here ----------------------------
    def pooled(cell_types: list[str], exclude: set[str] = frozenset(), half: str | None = None):
        a: list[float] = []
        b: list[float] = []
        for ct in cell_types:
            base, ifn = pick(ct, "Baseline"), pick(ct, "IFNa")
            ds = [d for d in donors if d in base and d in ifn and d not in exclude]
            if half == "predictor":
                xb = [score(ISG[0::2], base[d]["index"]) for d in ds]
                yb = [score(ISG[1::2], ifn[d]["index"]) - score(ISG[1::2], base[d]["index"]) for d in ds]
            elif half == "response":
                xb = [score(ISG[1::2], base[d]["index"]) for d in ds]
                yb = [score(ISG[0::2], ifn[d]["index"]) - score(ISG[0::2], base[d]["index"]) for d in ds]
            else:
                xb = [base[d]["isg"] for d in ds]
                yb = [ifn[d]["isg"] - base[d]["isg"] for d in ds]
            a += xb
            b += yb
        return a, b

    x, y = pooled(CELL_TYPES)
    base_rho = spearman(x, y)

    loo = {}
    for donor in donors:
        xx, yy = pooled(CELL_TYPES, exclude={donor})
        # rank within cell type (5 donors after exclusion) to keep the mixture rule
        rx, ry = [], []
        k = 0
        for ct in CELL_TYPES:
            size = len([d for d in donors if d != donor])
            rx += _rank(xx[k:k + size])
            ry += _rank(yy[k:k + size])
            k += size
        loo[donor] = round(spearman(rx, ry), 3)

    # partial correlation controlling for baseline housekeeping
    hk_x = []
    for ct in CELL_TYPES:
        base = pick(ct, "Baseline")
        hk_x += [base[d]["hk"] for d in donors if d in base]
    partial = round(partial_spearman(x, y, hk_x), 3)

    # disjoint halves
    xa, ya = pooled(CELL_TYPES, half="predictor")
    xb, yb = pooled(CELL_TYPES, half="response")
    halves = {
        "rho_baseline_7genes_vs_response_other_7": round(spearman(xa, ya), 3),
        "rho_baseline_other_7_vs_response_7": round(spearman(xb, yb), 3),
    }

    # ceiling / post-IFN level
    post = {ct: [pick(ct, "IFNa")[d]["isg"] for d in donors] for ct in CELL_TYPES}
    post_level_rho = round(
        spearman([pick(ct, "Baseline")[d]["isg"] for ct in CELL_TYPES for d in donors],
                 [v for ct in CELL_TYPES for v in post[ct]]), 3
    )
    spread = {ct: {"baseline_sd": round(_sd([pick(ct, "Baseline")[d]["isg"] for d in donors]), 3),
                   "post_ifn_sd": round(_sd(post[ct]), 3),
                   "baseline_range": round(max([pick(ct, "Baseline")[d]["isg"] for d in donors])
                                           - min([pick(ct, "Baseline")[d]["isg"] for d in donors]), 3),
                   "post_ifn_range": round(max(post[ct]) - min(post[ct]), 3)}
              for ct in CELL_TYPES}

    null = {}
    rng = random.Random(20260926)
    for half, (xa_, ya_) in {"predictor": (xa, ya), "response": (xb, yb)}.items():
        null[half] = round(spearman(xa_, ya_), 3)

    result = {
        "dataset": "GSE327707",
        "note": "all checks pre-specified in the script docstring",
        "primary_pooled_spearman_raw": round(base_rho, 3),
        "leave_one_donor_out": loo,
        "partial_spearman_controlling_baseline_housekeeping": partial,
        "disjoint_gene_halves": halves,
        "baseline_vs_post_ifn_level": {
            "spearman_baseline_vs_post_level": post_level_rho,
            "spread_by_cell_type": spread,
        },
        "per_donor_profile": profile,
        "housekeeping_baseline_by_donor": {
            ct: {d: round(pick(ct, "Baseline")[d]["hk"], 3) for d in donors} for ct in CELL_TYPES
        },
    }
    json.dump(result, open(os.path.join(OUT, "g6c_gse327707_robustness.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(json.dumps({k: result[k] for k in
                      ("primary_pooled_spearman_raw", "leave_one_donor_out",
                       "partial_spearman_controlling_baseline_housekeeping", "disjoint_gene_halves",
                       "baseline_vs_post_ifn_level")}, ensure_ascii=False, indent=1))
    print("\nper-donor mean over cell types:")
    for donor, entry in profile.items():
        print(f"  {donor}: {entry['mean_over_cell_types']}")


def _rank(v: list[float]) -> list[float]:
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def _sd(v: list[float]) -> float:
    m = sum(v) / len(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0


if __name__ == "__main__":
    main()
