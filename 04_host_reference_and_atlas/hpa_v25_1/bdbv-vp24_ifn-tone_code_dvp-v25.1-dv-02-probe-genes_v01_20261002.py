# -*- coding: utf-8 -*-
"""dv_02_probe_genes: does the DVP group file actually carry the focal genes, and does
`Matched nCPM` have values for them?  Read-only; writes only inside this task dir."""
import csv
import io
import os
import zipfile

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "dvp_delta_alignment_20260926")
DVP = os.path.join(ROOT, "analysis", "protein_layer_scan_20260926", "raw", "hpa_dvp")
FOCAL = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1", "KPNA5"]


def rows_of(zpath, member):
    with zipfile.ZipFile(zpath) as zf:
        with zf.open(member) as fh:
            txt = io.TextIOWrapper(fh, encoding="utf-8")
            rdr = csv.reader(txt, delimiter="\t")
            for row in rdr:
                yield row


def main():
    os.makedirs(OUT, exist_ok=True)
    out = []
    gz = os.path.join(DVP, "dvp_cell_type_group_data.tsv.zip")
    genes_g = {}
    header = None
    for i, r in enumerate(rows_of(gz, "dvp_cell_type_group_data.tsv")):
        if i == 0:
            header = r
            continue
        gene, name, ctype, inten, ncpm = r[1], r[1], r[2], r[3], r[4]
        d = genes_g.setdefault(gene, {"name": name, "n": 0, "inten_pos": 0, "ncpm_pos": 0})
        d["n"] += 1
        if inten.strip():
            d["inten_pos"] += 1
        if ncpm.strip():
            d["ncpm_pos"] += 1
    out.append("group header: %s" % header)
    out.append("group genes: %d ; rows: %d" % (len(genes_g), sum(v["n"] for v in genes_g.values())))
    out.append("genes with ANY non-empty Matched nCPM: %d" % sum(1 for v in genes_g.values() if v["ncpm_pos"]))
    out.append("")
    out.append("--- focal genes in the GROUP file ---")
    for g in FOCAL:
        out.append("  %-8s present=%s  %s" % (g, bool(genes_g.get(g)), genes_g.get(g)))
    out.append("")
    pop = sorted(g for g, v in genes_g.items() if v["ncpm_pos"] > 0)
    out.append("--- genes with populated Matched nCPM (first 30 of %d) ---" % len(pop))
    for g in pop[:30]:
        out.append("  %-12s %s  ncpm_pos=%d/%d" % (g, genes_g[g]["name"], genes_g[g]["ncpm_pos"], genes_g[g]["n"]))
    cz = os.path.join(DVP, "dvp_cell_type.tsv.zip")
    genes_c = {}
    for i, r in enumerate(rows_of(cz, "dvp_cell_type.tsv")):
        if i == 0:
            out.append("")
            out.append("cell_type header: %s" % r)
            continue
        gene, name, ctype, inten = r[1], r[1], r[2], r[3]
        d = genes_c.setdefault(gene, {"name": name, "n": 0, "inten_pos": 0})
        d["n"] += 1
        if inten.strip():
            d["inten_pos"] += 1
    out.append("cell_type genes: %d ; rows: %d" % (len(genes_c), sum(v["n"] for v in genes_c.values())))
    out.append("--- focal genes in the CELL_TYPE file ---")
    for g in FOCAL:
        out.append("  %-8s present=%s  %s" % (g, bool(genes_c.get(g)), genes_c.get(g)))
    gt, ct = [], []
    for i, r in enumerate(rows_of(gz, "dvp_cell_type_group_data.tsv")):
        if i and r[2] not in gt:
            gt.append(r[2])
    for i, r in enumerate(rows_of(cz, "dvp_cell_type.tsv")):
        if i and r[2] not in ct:
            ct.append(r[2])
    out.append("")
    out.append("GROUP cell types (%d): %s" % (len(gt), gt))
    out.append("CELL_TYPE cell types (%d): %s" % (len(ct), ct))
    for g in ["IL10RB", "IFNAR1", "IFNLR1", "KPNA5"]:
        out.append("")
        out.append("--- group-file rows for %s ---" % g)
        for i, r in enumerate(rows_of(gz, "dvp_cell_type_group_data.tsv")):
            if i and r[1] == g:
                out.append("  %s" % "\t".join(r))
    text = "\n".join(out)
    print(text)
    with open(os.path.join(OUT, "dv_02_probe_genes.txt"), "w", encoding="utf-8") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
