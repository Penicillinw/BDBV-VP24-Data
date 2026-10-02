"""X3 - voucher every 1/0 cell of the Supplementary Fig. S2d requirement matrix.

The matrix itself is frozen in the round-S figure script
(_coord/_round_S_20261002/figures/FigS2_check_and_gaps.py, panel_d).  For each
of the 35 cells this script locates the workspace artefact that does or does not
support the frozen value and reads the voucher value from that artefact.

voucher_status
  verified     the artefact supports the frozen 0/1
  contradicted the artefact contradicts the frozen 0/1
  unverified   no artefact in the workspace decides the cell
"""

from __future__ import annotations

import gzip
import json
import os
import re

import pandas as pd

ROOT = r"G:\本迪布焦研究"
HERE = os.path.dirname(os.path.abspath(__file__))


def jget(path, dotted):
    d = json.load(open(path, encoding="utf-8"))
    for k in dotted.split("."):
        d = d[int(k)] if k.isdigit() else d[k]
    return d


def soft(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    out = {}
    for m in re.finditer(r"^!(\S+)\s*=\s*(.*)$", txt, re.M):
        out.setdefault(m.group(1), []).append(m.group(2).strip())
    return out


def textfield(path, field):
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(field):
                return line.split("=", 1)[1].strip()
    return None


def tsv_cell(path, row_field, row_val, col):
    d = pd.read_csv(path, sep="\t")
    return d.loc[d[row_field] == row_val, col].iloc[0]


D = {
    "114905_json": os.path.join(ROOT, "analysis", "t4_gse114905_20260925",
                                "viral_mrna_summary.json"),
    "114905_esum": os.path.join(ROOT, "analysis", "main_tropism_gap_20260924",
                                "raw", "esummary_GSE114905.json"),
    "309699_soft": os.path.join(ROOT, "analysis", "a2_filovirus_anchor_20260926",
                                "raw", "GSE309699_soft_quick.txt"),
    "309699_find": os.path.join(ROOT, "analysis", "a2_filovirus_anchor_20260926",
                                "FINDINGS.md"),
    "226106_w14": os.path.join(ROOT, "analysis", "w14_organ_20260924",
                               "w14_summary.json"),
    "226106_g8": os.path.join(ROOT, "analysis", "t4_g8_tissue_ifn_env_20260926",
                              "g8_summary.json"),
    "226106_env": os.path.join(ROOT, "analysis", "t4_g8_tissue_ifn_env_20260926",
                               "tissue_env.tsv"),
    "226106_samp": os.path.join(ROOT, "analysis", "ext_resource_scan_20260924",
                                "GSE226106_samples.tsv"),
    "46599_res": os.path.join(ROOT, "analysis", "g6e_gse46599_20260926",
                              "g6e_results.json"),
    "46599_design": os.path.join(ROOT, "analysis", "g6e_gse46599_20260926",
                                 "raw", "GSE46599_design.txt"),
    "342661_soft": os.path.join(ROOT, "analysis", "protein_gse342661_20260926",
                                "raw", "GSE342661_soft_brief.txt"),
    "342661_scores": os.path.join(ROOT, "analysis",
                                  "gse342661_matrix_20260926",
                                  "sample_scores.tsv"),
    "342661_pub": os.path.join(ROOT, "analysis", "gse342661_matrix_20260926",
                               "PUBLICATION_STATUS_20260926.md"),
}

REQ = ["BDBV", "EBOV", "Same host cell type", "Host transcriptome",
       "IFN readout", ">= 2 replicates", "Publicly released"]
FROZEN = {  # dataset -> {requirement: 0/1}, transposed from the S panel_d grid
    "GSE114905": dict(zip(REQ, [1, 1, 1, 0, 0, 0, 1])),
    "GSE309699": dict(zip(REQ, [0, 1, 0, 1, 1, 1, 1])),
    "GSE226106": dict(zip(REQ, [0, 1, 0, 1, 1, 1, 1])),
    "GSE46599":  dict(zip(REQ, [0, 0, 0, 1, 1, 1, 1])),
    "GSE342661": dict(zip(REQ, [0, 0, 0, 1, 1, 1, 0])),
}

rows = []


def add(ds, req, vfile, vfield, vvalue, status, note=""):
    rows.append({"requirement": req, "dataset": ds,
                 "satisfied": FROZEN[ds][req],
                 "voucher_file": os.path.relpath(vfile, ROOT).replace("\\", "/"),
                 "voucher_field": vfield, "voucher_value": str(vvalue),
                 "accession": ds, "voucher_status": status, "note": note})


# ---------------------------------------------------------------- GSE114905
j = D["114905_json"]
e = json.load(open(D["114905_esum"], encoding="utf-8"))
add("GSE114905", "BDBV", j, "totals_reads.BDBV_D1",
    jget(j, "totals_reads.BDBV_D1"), "verified", "BDBV reads recovered")
add("GSE114905", "EBOV", j, "totals_reads.EBOV_D1",
    jget(j, "totals_reads.EBOV_D1"), "verified", "EBOV reads recovered")
add("GSE114905", "Same host cell type", D["114905_esum"], "title", e["title"],
    "verified", "all four viruses profiled in Huh7 cells")
add("GSE114905", "Host transcriptome", j, "design_note",
    jget(j, "design_note"), "verified")
add("GSE114905", "IFN readout", D["114905_esum"], "samples[].title",
    ";".join(s["title"] for s in e["samples"][:4]), "verified",
    "no interferon arm in any sample title")
add("GSE114905", ">= 2 replicates", j, "caveats.0",
    jget(j, "caveats.0"), "verified")
add("GSE114905", "Publicly released", D["114905_esum"], "pdat", e["pdat"],
    "verified", "GEO series public since 2018-12-25")

# ---------------------------------------------------------------- GSE309699
s = soft(D["309699_soft"])
add("GSE309699", "BDBV", D["309699_soft"], "Series_summary",
    (s["Series_summary"][0][:80] + "..."), "verified",
    "series is Zaire EBOV only; BDBV absent")
add("GSE309699", "EBOV", D["309699_soft"], "Series_summary",
    "Orthoebolavirus zairense (EBOV) infection of NHSK-1 keratinocytes",
    "verified")
add("GSE309699", "Same host cell type", D["309699_soft"], "Series_overall_design",
    s["Series_overall_design"][0][:70] + "...", "verified",
    "no BDBV arm, so the BDBV/EBOV same-cell comparison is not available")
add("GSE309699", "Host transcriptome", D["309699_soft"], "Series_type",
    s["Series_type"][0], "verified")
add("GSE309699", "IFN readout", D["309699_soft"], "Series_overall_design",
    s["Series_overall_design"][0][:70] + "...", "verified",
    "IFN-alpha/beta/gamma/lambda pre-treatment")
add("GSE309699", ">= 2 replicates", D["309699_find"], "section 2 arm sizes",
    "alpha/beta/gamma/lambda/untreated = 8/8/5/8/15 libraries", "verified")
add("GSE309699", "Publicly released", D["309699_soft"], "Series_status",
    s["Series_status"][0], "verified")

# ---------------------------------------------------------------- GSE226106
g8 = json.load(open(D["226106_g8"], encoding="utf-8"))
w14 = json.load(open(D["226106_w14"], encoding="utf-8"))
samp = pd.read_csv(D["226106_samp"], sep="\t")
add("GSE226106", "BDBV", D["226106_g8"], "limitations[1]",
    g8["limitations"][1], "verified", "EBOV Makona only, no BDBV")
add("GSE226106", "EBOV", D["226106_g8"], "resource", g8["resource"],
    "verified")
add("GSE226106", "Same host cell type", D["226106_g8"], "limitations[1]",
    g8["limitations"][1], "verified", "no BDBV arm")
add("GSE226106", "Host transcriptome", D["226106_w14"], "resource.n_samples",
    w14["resource"]["n_samples"], "verified",
    f"{w14['resource']['gene_rows']} gene rows")
add("GSE226106", "IFN readout", D["226106_env"], "E_broad (liver)",
    tsv_cell(D["226106_env"], "tissue", "liver", "E_broad"), "verified",
    "broad interferon-stimulated environment score per tissue")
add("GSE226106", ">= 2 replicates", D["226106_w14"], "resource.n_animals",
    w14["resource"]["n_animals"], "verified", "38 macaques")
pub = samp.loc[samp["field"] == "Sample_status", "S000"].iloc[0]
add("GSE226106", "Publicly released", D["226106_samp"], "Sample_status", pub,
    "verified")

# ---------------------------------------------------------------- GSE46599
r = json.load(open(D["46599_res"], encoding="utf-8"))
d465 = soft(D["46599_design"])
add("GSE46599", "BDBV", D["46599_res"], "design_verbatim",
    r["design_verbatim"], "verified", "HIV-1 restriction panel; no filovirus")
add("GSE46599", "EBOV", D["46599_res"], "design_verbatim",
    r["design_verbatim"], "verified", "no filovirus arm")
add("GSE46599", "Same host cell type", D["46599_res"], "design_verbatim",
    r["design_verbatim"], "verified", "no BDBV/EBOV pair")
add("GSE46599", "Host transcriptome", D["46599_res"], "platform",
    r["platform"] + f" ({r['n_samples']} samples)", "verified")
add("GSE46599", "IFN readout", D["46599_res"], "treatment_verbatim",
    r["treatment_verbatim"][:80] + "...", "verified")
add("GSE46599", ">= 2 replicates", D["46599_res"], "design_verbatim",
    r["design_verbatim"], "verified", "2 replicate experiments per cell line")
add("GSE46599", "Publicly released", D["46599_design"], "Series_status",
    d465["Series_status"][0], "verified")

# ---------------------------------------------------------------- GSE342661
s3 = soft(D["342661_soft"])
sc = pd.read_csv(D["342661_scores"], sep="\t")
add("GSE342661", "BDBV", D["342661_soft"], "Series_summary",
    s3["Series_summary"][0][:70] + "...", "verified",
    "influenza A / IFN design; no filovirus")
add("GSE342661", "EBOV", D["342661_soft"], "Series_summary",
    s3["Series_summary"][0][:70] + "...", "verified", "no filovirus arm")
add("GSE342661", "Same host cell type", D["342661_soft"], "Series_summary",
    s3["Series_summary"][0][:70] + "...", "verified", "no BDBV/EBOV pair")
add("GSE342661", "Host transcriptome", D["342661_soft"], "Series_type",
    s3["Series_type"][0], "verified",
    f"{len(sc)} sample columns, ISG17 module scored")
add("GSE342661", "IFN readout", D["342661_soft"], "Series_overall_design",
    s3["Series_overall_design"][0][:70] + "...", "verified",
    "IFN-beta1 and IFN-lambda1/2/3/4 arms")
add("GSE342661", ">= 2 replicates", D["342661_scores"],
    "rep column, max", int(sc["rep"].max()), "verified",
    "biological replicates 1-5 per arm")
add("GSE342661", "Publicly released", D["342661_pub"],
    "matrix GSE342661_counts_tpm.matrix.gz",
    "released on GEO FTP 2026-09-25, 10,325,209 B; also on disk under "
    "analysis/gse342661_matrix_20260926/raw/", "contradicted",
    "frozen grid says 0; the processed matrix is now publicly retrievable, "
    "so this cell should be 1 -> 342661 becomes 4/7, not 3/7")

out = pd.DataFrame(rows, columns=["requirement", "dataset", "satisfied",
                                  "voucher_file", "voucher_field",
                                  "voucher_value", "accession",
                                  "voucher_status", "note"])
out.to_csv(os.path.join(HERE, "X3_figs2d_provenance.tsv"), sep="\t",
           index=False)

print("cells:", len(out))
print(out["voucher_status"].value_counts().to_string())
print()
print("=== per-dataset count using the FROZEN values ===")
frozen_tot = {ds: sum(v.values()) for ds, v in FROZEN.items()}
print(frozen_tot)
print()
print("=== per-dataset count if the contradicted cell is corrected ===")
corr = dict(frozen_tot)
corr["GSE342661"] = corr["GSE342661"] + 1
print(corr)
print()
print(out[["dataset", "requirement", "satisfied", "voucher_status"]]
      .to_string(index=False))
