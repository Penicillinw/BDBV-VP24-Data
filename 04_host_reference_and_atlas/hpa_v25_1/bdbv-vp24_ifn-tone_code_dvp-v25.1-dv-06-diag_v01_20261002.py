# -*- coding: utf-8 -*-
"""dv_06_diag -- print the paired values behind P2 (which pairs are digit-identical and why).

Reads only; writes only inside this task directory.
"""
import csv
import io
import os
import zipfile

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "dvp_delta_alignment_20260926")
DVP = os.path.join(ROOT, "analysis", "protein_layer_scan_20260926", "raw", "hpa_dvp")
LAND = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")
FOCAL = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1", "KPNA5"]


def normalise(x):
    s = x.lower().strip().replace("(", " ").replace(")", " ").replace(",", " ")
    s = s.replace("/", " ").replace("-", " ")
    out = []
    for t in s.split():
        if len(t) > 3 and t.endswith("s") and not t.endswith("ss"):
            t = t[:-1]
        out.append(t)
    return " ".join(out)


def tokset(x):
    return {t for t in normalise(x).split() if t != "and"}


def main():
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
    hpa = {}
    with open(os.path.join(LAND, "hpa_ifn_landscape_wide.tsv"), encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        genes = next(rdr)[1:]
        for r in rdr:
            hpa[r[0]] = {g: (float(v) if v.strip() else None) for g, v in zip(genes, r[1:])}
    dts = list(hpa.keys())
    prim, back = {}, {}
    for d in types:
        h = [a for a in dts if normalise(a) == normalise(d)]
        if len(h) == 1:
            prim[d] = h[0]
            back.setdefault(h[0], []).append(d)
    prim = {d: a for d, a in prim.items() if len(back[a]) == 1}
    pairs = {d: [a] for d, a in prim.items()}
    for d in types:
        if d in pairs:
            continue
        td = tokset(d)
        h = [a for a in dts if td and tokset(a) and (td <= tokset(a) or tokset(a) <= td)]
        if len(h) == 1:
            pairs[d] = [h[0]]
    for d in types:
        if d in pairs:
            continue
        td = tokset(d)
        h = [a for a in dts if td and td < tokset(a)]
        if 1 <= len(h) <= 3:
            pairs[d] = sorted(h)
    lines = []
    for g in FOCAL:
        lines.append("=== %s ===" % g)
        lines.append("  %-26s %10s %10s  %s" % ("DVP group", "DVP nCPM", "HPA nCPM", "delta target"))
        for d, lst in pairs.items():
            dv = grp[g].get(d)
            vals = [hpa[a].get(g) for a in lst]
            if dv is None or any(v is None for v in vals):
                continue
            hv = sum(vals) / len(vals)
            flag = "   <-- identical" if abs(dv - hv) < 1e-9 else ""
            lines.append("  %-26s %10.1f %10.1f  %s%s" % (d, dv, hv, ";".join(lst), flag))
        lines.append("")
    text = "\n".join(lines)
    print(text)
    open(os.path.join(OUT, "p2_paired_values.txt"), "w", encoding="utf-8").write(text + "\n")


if __name__ == "__main__":
    main()
