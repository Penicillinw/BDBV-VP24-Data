"""G6 step 4: does baseline expression proxy interferon responsiveness?

Primary test inside a single platform and a single experiment (GSE21158):
for each of 10 human cancer cell lines, take the baseline (24 h control medium)
expression of the interferon system and the magnitude of the interferon-alpha-2a
response (24 h treatment minus 24 h control), and correlate the two.

Modules are frozen before any correlation is computed. A housekeeping module is
carried through the same pipeline as a negative control.

Outputs: analysis/g6_ifn_responsiveness_20260926/
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
from collections import defaultdict
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "g6_ifn_responsiveness_20260926")
RAW = os.path.join(OUT, "raw")
MATRIX = os.path.join(RAW, "GSE21158_series_matrix.txt.gz")
ANNOT = os.path.join(RAW, "GPL6098.annot.gz")

# frozen before any correlation was computed
MODULES = {
    "ISG_module": ["ISG15", "MX1", "MX2", "OAS1", "IFIT1", "IFIT2", "IFIT3",
                   "BST2", "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6"],
    "receptor_module": ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB"],
    "housekeeping_control": ["ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "HPRT1"],
}
LINEAGE = {
    "ME15": "melanocyte lineage (melanoma)",
    "D10": "melanocyte lineage (melanoma)",
    "JUSO": "melanocyte lineage (melanoma)",
    "LS174T": "colorectal epithelium",
    "HCT116": "colorectal epithelium",
    "AsPC1": "pancreatic epithelium",
    "MIAPACA": "pancreatic epithelium",
    "PANC1": "pancreatic epithelium",
    "A549": "lung epithelium",
    "Calu6": "lung epithelium",
}


def download(url: str, dest: str) -> str:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    last: Exception | None = None
    for attempt in range(5):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=600) as resp, open(dest, "wb") as fh:  # noqa: S310
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            return dest
        except Exception as exc:  # noqa: BLE001 - transient TLS/FTP resets
            last = exc
            if os.path.exists(dest):
                os.remove(dest)
            import time

            time.sleep(3 * (attempt + 1))
    raise SystemExit(f"download failed: {url}\n{last}")


def load_annotation() -> dict[str, str]:
    """Probe -> symbol from the complete GPL6098 platform table.

    GPL6098 (Illumina humanRef-8 v1.0) has no annot/ directory and the SOFT family
    download is unstable, so the platform table is fetched once from the GEO text
    view into raw/GPL6098_platform_table.txt; its ID column matches all 24,385
    series-matrix probe IDs exactly.
    """
    path = os.path.join(RAW, "GPL6098_platform_table.txt")
    if not os.path.exists(path):
        raise SystemExit(
            "missing raw/GPL6098_platform_table.txt - run "
            "tools/g6_fetch_gpl6098_table_20260926.py first"
        )
    probe2gene: dict[str, str] = {}
    header: list[str] | None = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                continue
            if header and line.strip() and not line.startswith("!"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < len(header):
                    parts += [""] * (len(header) - len(parts))
                row = dict(zip(header, parts))
                sym = (row.get("Symbol") or "").strip()
                if sym and sym not in {"---", ""}:
                    probe2gene[row["ID"].strip()] = sym.split("///")[0].strip()
    return probe2gene


def load_matrix() -> tuple[list[str], dict[str, list[float]]]:
    titles: list[str] = []
    rows: dict[str, list[float]] = {}
    in_table = False
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                next(fh)
            elif in_table:
                if line.startswith("!series_matrix_table_end"):
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                probe = parts[0].strip('"')
                try:
                    rows[probe] = [float(x) for x in parts[1:]]
                except ValueError:
                    continue
    return titles, rows


def parse_design(titles: list[str]) -> list[dict[str, str]]:
    out = []
    for t in titles:
        low = t.lower()
        line = next((k for k in LINEAGE if low.startswith(k.lower())), low.split()[0])
        cond = "IFN24" if ("24 hr interferon" in low) else (
            "IFN4" if ("4 hr interferon" in low) else "control")
        out.append({"title": t, "cell_line": line, "condition": cond})
    return out


def spearman(a: list[float], b: list[float]) -> float:
    def rank(x: list[float]) -> list[float]:
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


def permutation_p(a: list[float], b: list[float], n_perm: int = 20000, seed: int = 20260926) -> float:
    import random

    rng = random.Random(seed)
    obs = abs(spearman(a, b))
    hits = 0
    bb = list(b)
    for _ in range(n_perm):
        rng.shuffle(bb)
        if abs(spearman(a, bb)) >= obs:
            hits += 1
    return (hits + 1) / (n_perm + 1)


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    probe2gene = load_annotation()
    titles, rows = load_matrix()
    design = parse_design(titles)
    print(f"annotation probes: {len(probe2gene)}; matrix probes: {len(rows)}; samples: {len(titles)}")
    n_samp = len(titles)
    assert all(len(v) == n_samp for v in rows.values()), "ragged matrix"
    vals = [v for v in rows.values()]
    flat = [x for v in vals[:200] for x in v]
    print(f"value range: {min(flat):.2f} .. {max(flat):.2f} (looks "
          f"{'log-scale' if max(flat) < 30 else 'linear'})")

    gene_series: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for probe, vec in rows.items():
        sym = probe2gene.get(probe)
        if not sym:
            continue
        for i, d in enumerate(design):
            gene_series[sym][d["cell_line"]].append(vec[i]) if False else None
    # per gene, per cell line, per condition
    acc: dict[str, dict[str, dict[str, list[float]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for probe, vec in rows.items():
        sym = probe2gene.get(probe)
        if not sym:
            continue
        for i, d in enumerate(design):
            acc[sym][d["cell_line"]][d["condition"]].append(vec[i])

    lines = sorted({d["cell_line"] for d in design})
    def mean(xs: list[float]) -> float:
        return sum(xs) / len(xs) if xs else float("nan")

    rows_out = []
    for line in lines:
        rec: dict[str, object] = {"cell_line": line, "lineage": LINEAGE.get(line, "")}
        for mod, genes in MODULES.items():
            base_vals, ifn_vals = [], []
            for g in genes:
                if g in acc and line in acc[g]:
                    base_vals.append(mean(acc[g][line].get("control", [])))
                    ifn_vals.append(mean(acc[g][line].get("IFN24", [])))
            rec[f"{mod}_baseline"] = mean(base_vals)
            rec[f"{mod}_IFN24"] = mean(ifn_vals)
            rec[f"{mod}_induction"] = mean(ifn_vals) - mean(base_vals)
            rec[f"{mod}_n_genes"] = len(base_vals)
        rows_out.append(rec)

    with open(os.path.join(OUT, "gse21158_per_cell_line.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    def col(name: str) -> list[float]:
        return [float(r[name]) for r in rows_out]

    tests = {}
    for base_mod in ["ISG_module", "receptor_module", "housekeeping_control"]:
        for target in ["ISG_module"]:
            a = col(f"{base_mod}_baseline")
            b = col(f"{target}_induction")
            rho = spearman(a, b)
            tests[f"{base_mod}_baseline_vs_{target}_induction"] = {
                "n_cell_lines": len(lines),
                "spearman": round(rho, 3),
                "permutation_p": round(permutation_p(a, b), 4),
                "values_baseline": [round(x, 3) for x in a],
                "values_induction": [round(x, 3) for x in b],
            }
    # also: baseline vs absolute induced level (less exposed to regression to the mean)
    a = col("ISG_module_baseline")
    b = col("ISG_module_IFN24")
    tests["ISG_module_baseline_vs_ISG_module_IFN24_level"] = {
        "n_cell_lines": len(lines),
        "spearman": round(spearman(a, b), 3),
        "permutation_p": round(permutation_p(a, b), 4),
    }

    summary = {
        "dataset": "GSE21158",
        "design": "10 human cancer cell lines; 24 h control medium, 4 h and 24 h "
                  "IFN-alpha-2a (100 U/ml); triplicates; Illumina GPL6098",
        "note": "the 4 h time point has no matched 4 h control, so the primary "
                "comparison uses 24 h IFN versus 24 h control",
        "modules_frozen": MODULES,
        "lineage_map": LINEAGE,
        "tests": tests,
        "cell_lines": lines,
    }
    with open(os.path.join(OUT, "g6_results.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    print("\ncell line        ISGbase  ISGind  recbase  HKbase")
    for r in rows_out:
        print(f"{r['cell_line']:10s} {r['ISG_module_baseline']:8.3f} "
              f"{r['ISG_module_induction']:8.3f} {r['receptor_module_baseline']:8.3f} "
              f"{r['housekeeping_control_baseline']:8.3f}")
    print()
    for k, v in tests.items():
        print(f"{k}: rho={v['spearman']:+.3f} p={v['permutation_p']}")


if __name__ == "__main__":
    main()
