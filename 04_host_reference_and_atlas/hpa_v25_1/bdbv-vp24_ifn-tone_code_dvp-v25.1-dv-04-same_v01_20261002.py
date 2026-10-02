# -*- coding: utf-8 -*-
"""dv_04_same_source -- decisive test of whether DVP `Matched nCPM` is the SAME matrix as
the HPA single-cell consensus matrix that Delta was built from.

Method (frozen in PREREG v2): for EVERY gene shared by the two matrices, compute the
Spearman rank correlation across the matched DVP groups.  If the two columns come from one
and the same underlying matrix, per-gene rho must be ~1 for essentially all genes and a
large fraction must agree digit-for-digit.  If they are distinct matrices, rho is high but
materially below 1 and exact agreement is rare.

Also emits min-detectable-rho for every n used anywhere in this task.

Writes only inside analysis/dvp_delta_alignment_20260926/.
"""
import csv
import io
import json
import math
import os
import zipfile

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "dvp_delta_alignment_20260926")
DVP = os.path.join(ROOT, "analysis", "protein_layer_scan_20260926", "raw", "hpa_dvp")
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


def min_detectable_rho(n, power=0.80, alpha=0.05):
    if n < 4:
        return float("nan")
    return math.tanh((1.959963985 + 0.8416212336) / math.sqrt(n - 3))


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
    # ---- DVP group matrix
    grp = {}
    types = []
    with zipfile.ZipFile(os.path.join(DVP, "dvp_cell_type_group_data.tsv.zip")) as zf:
        with zf.open("dvp_cell_type_group_data.tsv") as fh:
            rdr = csv.reader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t")
            next(rdr)
            for r in rdr:
                gene, ctype, ncpm = r[1], r[2], r[4]
                if ctype not in types:
                    types.append(ctype)
                grp.setdefault(gene, {})[ctype] = float(ncpm) if ncpm.strip() else None
    # ---- Delta-side matrix
    with open(os.path.join(LAND, "hpa_ifn_landscape_wide.tsv"), encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        genes = hdr[1:]
        delta = {}
        for r in rdr:
            delta[r[0]] = {g: (float(v) if v.strip() else None) for g, v in zip(genes, r[1:])}
    delta_types = list(delta.keys())

    # ---- crosswalk (sens_a variant, matching dv_03)
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
    # families
    for d in types:
        if d in pairs:
            continue
        td = tokset(d)
        hits = [a for a in delta_types if td and td < tokset(a)]
        if 1 <= len(hits) <= 3:
            pairs[d] = sorted(hits)

    shared = [g for g in genes if g in grp]
    rhos, exact_all, n_used = [], 0, 0
    per_gene = []
    for g in shared:
        xs, ys = [], []
        for d, lst in pairs.items():
            dv = grp[g].get(d)
            vals = [delta[a].get(g) for a in lst]
            if dv is None or any(v is None for v in vals):
                continue
            xs.append(dv)
            ys.append(sum(vals) / len(vals))
        if len(xs) < 3:
            continue
        n_used += 1
        r = spearman(xs, ys)
        if not math.isnan(r):
            rhos.append(r)
        if all(abs(a - b) < 1e-9 for a, b in zip(xs, ys)):
            exact_all += 1
        per_gene.append((g, len(xs), r))
    rhos.sort()
    n = len(rhos)
    report = {
        "n_shared_genes": len(shared),
        "n_genes_used_ge3": n_used,
        "n_rho_finite": n,
        "median_rho": rhos[n // 2] if n else None,
        "q05_rho": rhos[int(0.05 * n)] if n else None,
        "q25_rho": rhos[int(0.25 * n)] if n else None,
        "q75_rho": rhos[int(0.75 * n)] if n else None,
        "q95_rho": rhos[int(0.95 * n)] if n else None,
        "frac_rho_ge_099": (sum(1 for r in rhos if r >= 0.99) / n) if n else None,
        "frac_rho_ge_095": (sum(1 for r in rhos if r >= 0.95) / n) if n else None,
        "frac_identical_all_pairs": exact_all / n_used if n_used else None,
        "n_matched_groups": len(pairs),
    }
    json.dump(report, open(os.path.join(OUT, "p5_same_source.json"), "w", encoding="utf-8"),
              indent=2, ensure_ascii=False)
    with open(os.path.join(OUT, "p5_same_source_pergene.tsv"), "w", encoding="utf-8") as fh:
        fh.write("gene\tn_groups\tspearman\n")
        for g, k, r in per_gene:
            fh.write("%s\t%d\t%.6f\n" % (g, k, r))
    print(json.dumps(report, indent=2))

    powers = {}
    for k in [5, 6, 8, 9, 10, 13, 14, 15, 17, 19, 21, 24, 39]:
        powers[k] = min_detectable_rho(k)
    json.dump(powers, open(os.path.join(OUT, "p6_min_detectable_rho.json"), "w",
                           encoding="utf-8"), indent=2)
    print("\nmin detectable |rho| at 80% power, alpha=0.05 (two-sided):")
    for k, v in powers.items():
        print("   n=%2d -> |rho| >= %.4f" % (k, v))


if __name__ == "__main__":
    main()
