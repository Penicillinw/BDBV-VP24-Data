# -*- coding: utf-8 -*-
"""dv_05_genomewide_same_source -- genome-wide version of the P5 same-source test.

The 29-gene Delta panel is too narrow to settle whether DVP `Matched nCPM` is the same
matrix as the HPA single-cell consensus that Delta was built from.  This script streams the
raw HPA single-cell table (data/hpa/rna_single_cell_type.tsv, 148 MB) and keeps only the
Delta cell types that the crosswalk matched, then compares every shared gene.

Reads only.  Writes only inside analysis/dvp_delta_alignment_20260926/.
"""
import csv
import io
import json
import math
import os
import sys
import zipfile

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "dvp_delta_alignment_20260926")
DVP = os.path.join(ROOT, "analysis", "protein_layer_scan_20260926", "raw", "hpa_dvp")
HPA = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")
LAND = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")


def rankdata(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def pearson(a, b):
    n = len(a)
    if n < 3:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    if da == 0 or db == 0:
        return float("nan")
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (da * db)


def spearman(a, b):
    return pearson(rankdata(a), rankdata(b))


def normalise(name):
    s = name.lower().strip().replace("(", " ").replace(")", " ").replace(",", " ")
    s = s.replace("/", " ").replace("-", " ")
    toks = []
    for t in s.split():
        if len(t) > 3 and t.endswith("s") and not t.endswith("ss"):
            t = t[:-1]
        toks.append(t)
    return " ".join(toks)


def tokset(name):
    return {t for t in normalise(name).split() if t != "and"}


def main():
    # ---- DVP group matrix (24 groups; gene -> group -> matched nCPM)
    grp, types = {}, []
    with zipfile.ZipFile(os.path.join(DVP, "dvp_cell_type_group_data.tsv.zip")) as zf:
        with zf.open("dvp_cell_type_group_data.tsv") as fh:
            rdr = csv.reader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t")
            next(rdr)
            for r in rdr:
                gene, ctype, ncpm = r[1], r[2], r[4]
                if ctype not in types:
                    types.append(ctype)
                grp.setdefault(gene, {})[ctype] = float(ncpm) if ncpm.strip() else None
    with open(os.path.join(LAND, "hpa_ifn_landscape_wide.tsv"), encoding="utf-8") as fh:
        hdr = next(csv.reader(fh, delimiter="\t"))
        delta_types = [r[0] for r in csv.reader(open(
            os.path.join(LAND, "hpa_ifn_landscape_wide.tsv"), encoding="utf-8"), delimiter="\t")][1:]

    # crosswalk = sens_a + family (identical construction to dv_03)
    prim, back = {}, {}
    for d in types:
        hits = [a for a in delta_types if normalise(a) == normalise(d)]
        if len(hits) == 1:
            prim[d] = hits[0]
            back.setdefault(hits[0], []).append(d)
    prim = {d: a for d, a in prim.items() if len(back[a]) == 1}
    pairs = {d: [a] for d, a in prim.items()}
    for d in types:
        if d in pairs:
            continue
        td = tokset(d)
        hits = [a for a in delta_types if td and tokset(a) and
                (td <= tokset(a) or tokset(a) <= td)]
        if len(hits) == 1:
            pairs[d] = [hits[0]]
    for d in types:
        if d in pairs:
            continue
        td = tokset(d)
        hits = [a for a in delta_types if td and td < tokset(a)]
        if 1 <= len(hits) <= 3:
            pairs[d] = sorted(hits)

    wanted = set()
    for lst in pairs.values():
        wanted.update(lst)
    print("matched Delta cell types to pull from the raw HPA table: %d" % len(wanted),
          file=sys.stderr)

    # ---- stream the raw HPA single-cell table
    hpav = {}
    with open(HPA, encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        head = next(rdr)
        assert head[:4] == ["Gene", "Gene name", "Cell type", "nCPM"], head
        for i, r in enumerate(rdr):
            if r[2] in wanted:
                hpav.setdefault(r[1], {})[r[2]] = float(r[3]) if r[3].strip() else None
            if i % 500000 == 0:
                print("  ... %d rows" % i, file=sys.stderr)

    shared = sorted(set(hpav) & set(grp))
    print("genes shared between the two matrices: %d" % len(shared), file=sys.stderr)

    rhos, exact_all, used = [], 0, 0
    rows = ["gene\tn_groups\tspearman\tidentical_all"]
    for g in shared:
        xs, ys = [], []
        for d, lst in pairs.items():
            dv = grp[g].get(d)
            vals = [hpav[g].get(a) for a in lst]
            if dv is None or any(v is None for v in vals):
                continue
            xs.append(dv)
            ys.append(sum(vals) / len(vals))
        if len(xs) < 3:
            continue
        used += 1
        r = spearman(xs, ys)
        ident = all(abs(a - b) < 1e-9 for a, b in zip(xs, ys))
        if ident:
            exact_all += 1
        if not math.isnan(r):
            rhos.append(r)
        rows.append("%s\t%d\t%.6f\t%s" % (g, len(xs), r, "yes" if ident else "no"))
    rhos.sort()
    n = len(rhos)
    report = {
        "source": "data/hpa/rna_single_cell_type.tsv (raw, streamed)",
        "n_genes_shared": len(shared),
        "n_genes_used": used,
        "n_matched_groups": len(pairs),
        "median_rho": rhos[n // 2] if n else None,
        "q05_rho": rhos[int(0.05 * n)] if n else None,
        "q25_rho": rhos[int(0.25 * n)] if n else None,
        "q75_rho": rhos[int(0.75 * n)] if n else None,
        "q95_rho": rhos[int(0.95 * n)] if n else None,
        "frac_rho_ge_099": (sum(1 for r in rhos if r >= 0.99) / n) if n else None,
        "frac_rho_ge_095": (sum(1 for r in rhos if r >= 0.95) / n) if n else None,
        "frac_identical_all_pairs": exact_all / used if used else None,
    }
    json.dump(report, open(os.path.join(OUT, "p5b_genomewide_same_source.json"), "w",
                           encoding="utf-8"), indent=2, ensure_ascii=False)
    open(os.path.join(OUT, "p5b_genomewide_pergene.tsv"), "w", encoding="utf-8").write(
        "\n".join(rows) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
