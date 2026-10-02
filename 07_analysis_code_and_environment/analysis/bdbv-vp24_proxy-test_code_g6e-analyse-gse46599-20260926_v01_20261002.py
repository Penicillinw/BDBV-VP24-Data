"""Material test: within one experiment, are immortalized lines graded responders while
primary human cells respond coherently?

GSE46599: 9 immortalized lines (CEM, CEM-SS, HT1080, Jurkat, PMA-THP-1, PMA-U937, THP-1,
U87-MG, U937) and 2 primary cell types (macrophages, CD4+ T cells), each untreated and
IFN-treated, matched within experiment and platform (Illumina HumanHT-12 v4).
"""

from __future__ import annotations

import gzip
import json
import math
import os
import re
from collections import defaultdict
from urllib.request import Request, urlopen

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "g6e_gse46599_20260926")
RAW = os.path.join(OUT, "raw")
MATRIX = os.path.join(RAW, "GSE46599_series_matrix.txt.gz")
ANNOT = os.path.join(RAW, "GPL10558.annot.gz")

ISG = ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "IFIT1", "IFIT2", "IFIT3", "BST2",
       "RSAD2", "USP18", "STAT1", "IRF7", "IFI27", "IFI6", "DDX58", "IFIH1"]
PRIMARY = {"primary-macrophages", "primary-CD4+-T-cells"}


def fetch(url: str, dest: str) -> str:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    last = None
    for a in range(5):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=600) as r, open(dest, "wb") as fh:  # noqa: S310
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            return dest
        except Exception as exc:  # noqa: BLE001
            last = exc
            if os.path.exists(dest):
                os.remove(dest)
            import time

            time.sleep(3 * (a + 1))
    raise SystemExit(f"failed {url}: {last}")


def load_annot() -> dict[str, str]:
    path = fetch(
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
        ANNOT,
    )
    out: dict[str, str] = {}
    header: list[str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                continue
            if header is None or line.startswith(("!", "#")) or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sym = parts[2].strip()
            if sym and sym not in {"---", ""}:
                out[parts[0].strip()] = sym.split("///")[0].strip()
    return out


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    p2g = load_annot()
    print("annotation probes:", len(p2g))
    fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46599/matrix/"
          "GSE46599_series_matrix.txt.gz", MATRIX)

    titles: list[str] = []
    probes: dict[str, list[float]] = {}
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        in_t = False
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_t = True
                next(fh)
            elif in_t:
                if line.startswith("!series_matrix_table_end"):
                    break
                p = line.rstrip("\n").split("\t")
                if len(p) > 1:
                    try:
                        probes[p[0].strip('"')] = [float(x) for x in p[1:]]
                    except ValueError:
                        pass
    print("samples:", len(titles), "probes:", len(probes))

    gene_rows: dict[str, list[list[float]]] = defaultdict(list)
    for probe, vec in probes.items():
        sym = p2g.get(probe)
        if sym and sym in ISG:
            gene_rows[sym].append(vec)

    def sample_idx(cond: str, cell: str) -> list[int]:
        return [i for i, t in enumerate(titles) if t.startswith(cell + "_")
                and f"_{cond}_" in t]

    cells = []
    for t in titles:
        cell = t.rsplit("_", 2)[0]
        if cell not in cells:
            cells.append(cell)

    rows = []
    for cell in cells:
        base_i = sample_idx("None", cell)
        ifn_i = sample_idx("IFN", cell)
        # module score per sample = mean over available ISG probes (first probe per gene)
        def score(idxs: list[int]) -> float:
            vals = []
            for g in ISG:
                if gene_rows.get(g):
                    vec = gene_rows[g][0]
                    vals.append(sum(vec[i] for i in idxs) / len(idxs))
            return sum(vals) / len(vals) if vals else float("nan")

        b, f = score(base_i), score(ifn_i)
        rows.append({
            "cell": cell,
            "kind": "primary" if cell in PRIMARY else "immortalized",
            "n_base": len(base_i),
            "n_ifn": len(ifn_i),
            "baseline": round(b, 3),
            "ifn": round(f, 3),
            "induction": round(f - b, 3),
        })
    rows.sort(key=lambda r: r["induction"])
    print(f"\n{'cell':22s} {'kind':13s} {'baseline':>9s} {'IFN':>8s} {'induction':>10s}")
    for r in rows:
        print(f"{r['cell']:22s} {r['kind']:13s} {r['baseline']:9.3f} {r['ifn']:8.3f} "
              f"{r['induction']:10.3f}")

    prim = [r["induction"] for r in rows if r["kind"] == "primary"]
    imm = [r["induction"] for r in rows if r["kind"] == "immortalized"]
    nonresp = [r for r in rows if r["kind"] == "immortalized" and r["induction"] < 0.5]

    def spearman(a, b):
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

    base_all = [r["baseline"] for r in rows]
    ind_all = [r["induction"] for r in rows]
    base_imm = [r["baseline"] for r in rows if r["kind"] == "immortalized"]
    ind_imm = [r["induction"] for r in rows if r["kind"] == "immortalized"]
    spread = max(ind_all) - min(ind_all)
    out = {
        "dataset": "GSE46599",
        "design": "9 immortalized lines and 2 primary human cell types, untreated vs "
                  "type-I IFN treated, matched within experiment (Illumina HT-12 v4)",
        "rows": rows,
        "primary_inductions": prim,
        "immortalized_inductions": imm,
        "immortalized_non_responders_lt_0.5": [r["cell"] for r in nonresp],
        "n_immortalized_non_responders": len(nonresp),
        "primary_all_respond": all(x > 0.5 for x in prim),
        "induction_range": [round(min(ind_all), 3), round(max(ind_all), 3)],
        "induction_spread": round(spread, 3),
        "is_bimodal": bool(spread > 0 and
                           sorted(ind_all)[len(ind_all) // 2] - sorted(ind_all)[len(ind_all) // 2 - 1] > 0.5),
        "spearman_baseline_vs_induction_all11": round(spearman(base_all, ind_all), 3),
        "spearman_baseline_vs_induction_immortalized_only": round(spearman(base_imm, ind_imm), 3),
    }
    with open(os.path.join(OUT, "material_test.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print("\nprimary inductions:", prim)
    print("immortalized non-responders (<0.5):", out["immortalized_non_responders_lt_0.5"])


if __name__ == "__main__":
    main()
