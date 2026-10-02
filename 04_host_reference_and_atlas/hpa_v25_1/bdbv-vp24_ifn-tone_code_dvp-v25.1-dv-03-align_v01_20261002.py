# -*- coding: utf-8 -*-
"""dv_03_align -- DVP x Delta alignment (task `dvp_delta_alignment`).  FINAL pipeline.

FROZEN CRITERIA: PREREG.md (v1) + the revision recorded in §7 of PREREG (v2).

  Resolution note (established by dv_02_probe_genes.py, before any statistics):
    * dvp_cell_type.tsv      : 27 cell types, protein intensity only.
    * dvp_cell_type_group_data.tsv : 24 GROUPS, protein intensity AND `Matched nCPM`.
      The three collapsed pairs are (alveolar 1 + 2 -> alveolar cell types),
      (cd4 t + cd8 t -> t-cells), (secretory cells -> fallopian tube cells).
    * Column 0 is the Ensembl id; the gene SYMBOL is column 1.

  P1  three crosswalk variants (all token-set based, no arbitrary "first hit wins"):
        primary     = 1:1 for normalised-name equality
        sens_a      = primary + 1:1 for token-subset containment
        sensitivity = sens_a  + 1:many for DVP-tokens-subset-of-Delta-tokens families (<=3)
  P2  DVP `Matched nCPM` vs Delta-side HPA nCPM, six focal genes, with a SAME-SOURCE gate
      that runs BEFORE any correlation is interpreted.
  P3  recompute Delta' from DVP matched nCPM on the 24 groups; compare to Delta by RANK;
      10,000 within-pair label permutations.
  P4  protein testability: needs non-empty Intensity in >=5 matched groups.

Writes only inside analysis/dvp_delta_alignment_20260926/.
"""
import csv
import hashlib
import io
import json
import math
import os
import random
import sys
import zipfile

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "dvp_delta_alignment_20260926")
SCAN = os.path.join(ROOT, "analysis", "protein_layer_scan_20260926")
DVP = os.path.join(SCAN, "raw", "hpa_dvp")
LAND = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")

FOCAL = ["IFNAR1", "IFNAR2", "IFNLR1", "IL10RB", "ISG15", "MX1"]
EXTRA = ["KPNA5"]
ARMS = {
    "IFN_I_capacity": ["IFNAR1", "IFNAR2"],
    "IFN_III_capacity": ["IFNLR1", "IL10RB"],
    "ISG_priming": ["ISG15", "MX1"],
}
RNG_SEED = 20260927
N_PERM = 10000
LOG = []


def log(msg=""):
    print(msg)
    LOG.append(str(msg))


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
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    if da == 0 or db == 0:
        return float("nan")
    return num / (da * db)


def spearman(a, b):
    return pearson(rankdata(a), rankdata(b))


def perm_p(obs, null):
    ge = sum(1 for v in null if abs(v) >= abs(obs) - 1e-12)
    return (ge + 1) / (len(null) + 1)


def min_detectable_rho(n, power=0.80, alpha=0.05):
    if n < 4:
        return float("nan")
    return math.tanh((1.959963985 + 0.8416212336) / math.sqrt(n - 3))


def zscore(xs):
    n = len(xs)
    m = sum(xs) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))
    if sd == 0:
        return [0.0] * n
    return [(x - m) / sd for x in xs]


def read_tsv_zip(path, member):
    with zipfile.ZipFile(path) as zf:
        with zf.open(member) as fh:
            return list(csv.reader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t"))


def load_dvp_group():
    rows = read_tsv_zip(os.path.join(DVP, "dvp_cell_type_group_data.tsv.zip"),
                        "dvp_cell_type_group_data.tsv")
    assert rows[0][:5] == ["Gene", "Gene name", "Cell type name", "Intensity", "Matched nCPM"]
    data, types, ens = {}, [], {}
    for r in rows[1:]:
        ens_id, gene, ctype, inten, ncpm = r[0], r[1], r[2], r[3], r[4]
        ens[gene] = ens_id
        if ctype not in types:
            types.append(ctype)
        data.setdefault(gene, {})[ctype] = (float(inten) if inten.strip() else None,
                                            float(ncpm) if ncpm.strip() else None)
    return data, types, ens


def load_dvp_celltype():
    rows = read_tsv_zip(os.path.join(DVP, "dvp_cell_type.tsv.zip"), "dvp_cell_type.tsv")
    assert rows[0][:4] == ["Gene", "Gene name", "Cell type", "Intensity"]
    data, types = {}, []
    for r in rows[1:]:
        gene, ctype, inten = r[1], r[2], r[3]
        if ctype not in types:
            types.append(ctype)
        data.setdefault(gene, {})[ctype] = float(inten) if inten.strip() else None
    return data, types


def load_delta_wide():
    with open(os.path.join(LAND, "hpa_ifn_landscape_wide.tsv"), encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        genes = hdr[1:]
        out = {}
        for r in rdr:
            out[r[0]] = {g: (float(v) if v.strip() else None) for g, v in zip(genes, r[1:])}
    return out


def load_rg():
    with open(os.path.join(LAND, "restriction_gradient.tsv"), encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        out = {}
        for r in rdr:
            out[r[0]] = dict(zip(hdr[1:], r[1:]))
    return hdr, out


def load_existing():
    out = {}
    with open(os.path.join(SCAN, "dvp_vs_delta_match.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["dvp_cell_type"]] = (r["delta_match"], r["match_kind"])
    return out


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


def build_primary(dvp_types, delta_types):
    fwd, back = {}, {}
    for d in dvp_types:
        hits = [a for a in delta_types if normalise(a) == normalise(d)]
        if len(hits) == 1:
            fwd[d] = hits[0]
            back.setdefault(hits[0], []).append(d)
    return {d: a for d, a in fwd.items() if len(back[a]) == 1}


def build_subset(dvp_types, delta_types):
    fwd, back = {}, {}
    for d in dvp_types:
        td = tokset(d)
        hits = [a for a in delta_types if td and tokset(a) and
                (td <= tokset(a) or tokset(a) <= td)]
        if len(hits) == 1:
            fwd[d] = hits[0]
            back.setdefault(hits[0], []).append(d)
    return {d: a for d, a in fwd.items() if len(back[a]) == 1}


def build_family(dvp_types, delta_types, already):
    """1:many -- DVP tokens are a strict subset of the Delta name's tokens."""
    fam = {}
    for d in dvp_types:
        if d in already:
            continue
        td = tokset(d)
        hits = [a for a in delta_types if td and td < tokset(a)]
        if 1 <= len(hits) <= 3:
            fam[d] = sorted(hits)
    return fam


def sha256(path, n=16):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def pair_value(pairs, d, delta_wide, gene):
    """Delta-side value for a DVP group: mean over its Delta family (None if any missing)."""
    vals = [delta_wide[a].get(gene) for a in pairs[d]]
    if any(v is None for v in vals):
        return None
    return sum(vals) / len(vals)


def main():
    os.makedirs(OUT, exist_ok=True)
    grp, grp_types, ens = load_dvp_group()
    ct, ct_types = load_dvp_celltype()
    delta_wide = load_delta_wide()
    delta_types = list(delta_wide.keys())
    rg_hdr, rg = load_rg()
    dcol = next(c for c in rg_hdr[1:] if "delta" in c.lower())
    exist = load_existing()

    log("Delta cell types %d ; DVP groups %d ; DVP cell types %d" % (
        len(delta_types), len(grp_types), len(ct_types)))
    log("restriction_gradient Delta column: %r" % dcol)
    log("")

    prim = build_primary(grp_types, delta_types)
    sub = build_subset(grp_types, delta_types)
    sens_a = dict(prim)
    sens_a.update({d: a for d, a in sub.items() if d not in sens_a})
    fam = build_family(grp_types, delta_types, sens_a)
    sens = {d: [a] for d, a in sens_a.items()}
    for d, lst in fam.items():
        sens[d] = lst

    pairs_variants = {
        "primary": {d: [a] for d, a in prim.items()},
        "sens_a": {d: [a] for d, a in sens_a.items()},
        "sensitivity": sens,
    }
    for k, v in pairs_variants.items():
        log("P1 %-12s n=%d" % (k, len(v)))
    log("P1 existing table (usable rows) n=%d" % sum(
        1 for v in exist.values() if not v[1].startswith(("none", "loose"))))

    rows = ["dvp_cell_type\tdvp_protein_cell_types\tin_primary\tin_sens_a\tin_sensitivity\t"
            "delta_strict\tdelta_family\texisting_match\texisting_kind\tdisagreement"]
    for d in grp_types:
        fam_s = ";".join(sens.get(d, [])) if isinstance(sens.get(d), list) else sens.get(d, "")
        e, ek = exist.get(d, ("", ""))
        e_ok = e if not ek.startswith(("none", "loose")) else ""
        primary_delta = prim.get(d, "")
        disag = ""
        if primary_delta and e_ok and primary_delta != e_ok:
            disag = "primary=%s vs existing=%s" % (primary_delta, e_ok)
        elif not primary_delta and e_ok:
            disag = "existing-only (%s | %s)" % (e_ok, ek)
        elif primary_delta and not e_ok:
            disag = "primary-only (%s) existing=%s" % (primary_delta, ek)
        rows.append("\t".join([d,
                               "%d" % sum(1 for c in ct_types if tokset(c) & tokset(d)),
                               "yes" if d in prim else "no",
                               "yes" if d in sens_a else "no",
                               "yes" if d in sens else "no",
                               primary_delta, fam_s, e, ek, disag]))
    with open(os.path.join(OUT, "crosswalk.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")
    disag = [r for r in rows[1:] if r.split("\t")[9]]
    log("P1 rows disagreeing with the existing table: %d" % len(disag))
    for r in disag:
        log("    " + r.replace("\t", " | "))

    json.dump({
        "delta_cell_types": len(delta_types),
        "dvp_groups_24": grp_types,
        "dvp_cell_types_27": ct_types,
        "n_primary": len(prim), "n_sens_a": len(sens_a), "n_sensitivity": len(sens),
        "families": fam,
        "n_existing_usable": sum(1 for v in exist.values() if not v[1].startswith(("none", "loose"))),
        "n_disagreements": len(disag),
    }, open(os.path.join(OUT, "p1_crosswalk_summary.json"), "w", encoding="utf-8"),
        indent=2, ensure_ascii=False)

    # ------------------------------------------------------------------ P2
    p2 = ["panel\tgene\tn_groups\tn_both_present\tn_digit_identical\tsame_source_frac\t"
          "spearman\tsource_verdict\tconcordance_verdict"]
    p2_same = {}
    for panel, pairs in pairs_variants.items():
        for g in FOCAL + EXTRA:
            xs, ys, ident = [], [], 0
            for d in pairs:
                dv = grp.get(g, {}).get(d, (None, None))[1]
                av = pair_value(pairs, d, delta_wide, g)
                if dv is None or av is None:
                    continue
                xs.append(dv)
                ys.append(av)
                if abs(dv - av) < 1e-9:
                    ident += 1
            both = len(xs)
            frac = (ident / both) if both else float("nan")
            rho = spearman(xs, ys) if both >= 3 else float("nan")
            if both == 0:
                src, verdict = "NO-OVERLAP", "NO-OVERLAP"
            else:
                src = "SAME-SOURCE" if frac >= 0.95 else "DISTINCT-SOURCE"
                if src == "SAME-SOURCE":
                    verdict = "SAME-SOURCE (internal consistency only)"
                elif math.isnan(rho):
                    verdict = "NOT-TESTABLE (coverage)"
                elif rho >= 0.5:
                    verdict = "CONCORDANT"
                elif rho >= 0.2:
                    verdict = "WEAK"
                else:
                    verdict = "DISCORDANT"
            p2.append("\t".join([panel, g, str(len(pairs)), str(both), str(ident),
                                 ("%.4f" % frac if both else ""),
                                 ("%.4f" % rho if not math.isnan(rho) else ""), src, verdict]))
            p2_same["%s/%s" % (panel, g)] = {"n": both, "identical": ident, "frac": frac,
                                             "rho": rho, "source": src}
    with open(os.path.join(OUT, "p2_mrna_concordance.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(p2) + "\n")
    for line in p2:
        log("P2 " + line)
    # worked example for the same-source gate
    log("")
    log("P2 same-source worked example (sensitivity panel): macrophages / IFNLR1")
    a = sens.get("macrophages")
    if a:
        log("    DVP matched nCPM = %s ; Delta-side nCPM = %s" % (
            grp["IFNLR1"]["macrophages"][1], [delta_wide[x]["IFNLR1"] for x in (a if isinstance(a, list) else [a])]))

    # ------------------------------------------------------------------ P3
    complete = [d for d in grp_types
                if all(grp.get(g, {}).get(d, (None, None))[1] is not None for g in FOCAL)]
    log("")
    log("P3 DVP groups with all six focal Matched nCPM present: %d / %d" % (len(complete), len(grp_types)))
    zz = {g: dict(zip(complete, zscore([grp[g][d][1] for d in complete]))) for g in FOCAL}
    delta_p = {d: sum(sum(zz[g][d] for g in ARMS[arm]) / len(ARMS[arm]) for arm in ARMS)
               for d in complete}

    p3 = ["panel\tarm\tn_pairs\tspearman\tnull_mean\tnull_p95\tperm_p\tverdict"]
    p3_json = {}
    for panel, pairs in pairs_variants.items():
        keep = {}
        for d, lst in pairs.items():
            if d not in delta_p:
                continue
            vals = [rg.get(a, {}).get(dcol, "") for a in lst]
            if any(v == "" for v in vals):
                continue
            keep[d] = sum(float(v) for v in vals) / len(vals)
        n = len(keep)
        if n < 4:
            p3.append("\t".join([panel, "Delta", str(n), "", "", "", "", "NOT-TESTABLE (coverage)"]))
            continue
        obs = [delta_p[d] for d in keep]
        ref = [keep[d] for d in keep]
        rho = spearman(obs, ref)
        rnd = random.Random(RNG_SEED)
        null = []
        for _ in range(N_PERM):
            sh = obs[:]
            rnd.shuffle(sh)
            null.append(spearman(sh, ref))
        p = perm_p(rho, null)
        srt = sorted(null)
        sm = sum(null) / len(null)
        q95 = srt[max(0, int(0.95 * len(srt)) - 1)]
        if rho >= 0.5 and p < 0.05:
            v = "ALIGNED"
        elif rho >= 0.5:
            v = "ALIGNED-BUT-UNDERPOWERED"
        elif rho >= 0.2:
            v = "PARTIAL"
        else:
            v = "NOT-ALIGNED"
        p3.append("\t".join([panel, "Delta", str(n), "%.4f" % rho, "%.4f" % sm,
                             "%.4f" % q95, "%.5f" % p, v]))
        p3_json[panel] = {"n_pairs": n, "spearman": rho, "null_mean": sm, "null_p95": q95,
                          "perm_p": p, "verdict": v, "min_detectable_rho": min_detectable_rho(n),
                          "n_perm": N_PERM, "seed": RNG_SEED}
        for arm in list(ARMS):
            ov = [sum(zz[g][d] for g in ARMS[arm]) / len(ARMS[arm]) for d in keep]
            p3.append("\t".join([panel, arm, str(n), "%.4f" % spearman(ov, ref),
                                 "", "", "", "diag"]))
    with open(os.path.join(OUT, "p3_delta_recompute.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(p3) + "\n")
    json.dump(p3_json, open(os.path.join(OUT, "p3_delta_recompute.json"), "w", encoding="utf-8"),
              indent=2, ensure_ascii=False)
    for line in p3:
        log("P3 " + line)

    # ------------------------------------------------------------------ P4
    p4 = ["panel\tgene\tn_nonempty_matched_24\tn_nonempty_all_24\t"
          "n_nonempty_27celltype\tprotein_testable\trho_protein_vs_matched_ncpm\t"
          "rho_protein_vs_delta\tverdict"]
    for panel, pairs in pairs_variants.items():
        for g in FOCAL + EXTRA:
            det = [d for d in pairs if grp.get(g, {}).get(d, (None, None))[0] is not None]
            det_all = sum(1 for d in grp_types if grp.get(g, {}).get(d, (None, None))[0] is not None)
            det27 = sum(1 for d in ct_types if ct.get(g, {}).get(d) is not None)
            if len(det) < 5:
                p4.append("\t".join([panel, g, str(len(det)), str(det_all), str(det27), "no",
                                     "", "", "NOT-TESTABLE (coverage)"]))
                continue
            xs, ys, ds = [], [], []
            for d in det:
                xs.append(grp[g][d][0])
                ys.append(grp[g][d][1])
                vals = [rg.get(a, {}).get(dcol, "") for a in pairs[d]]
                if any(v == "" for v in vals):
                    continue
                ds.append(sum(float(v) for v in vals) / len(vals))
            r1 = spearman(xs, ys) if len(xs) >= 3 else float("nan")
            r2 = spearman(xs, ds) if len(ds) >= 3 else float("nan")
            p4.append("\t".join([panel, g, str(len(det)), str(det_all), str(det27), "yes",
                                 ("%.4f" % r1 if not math.isnan(r1) else ""),
                                 ("%.4f" % r2 if not math.isnan(r2) else ""), "TESTED"]))
    with open(os.path.join(OUT, "p4_protein_layer.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(p4) + "\n")
    for line in p4:
        log("P4 " + line)

    # ------------------------------------------------------------------ repro
    rp = ["# repro -- dvp_delta_alignment", "",
          "python %s" % sys.version.split()[0], "",
          "```", "python analysis/dvp_delta_alignment_20260926/dv_01_probe.py",
          "python analysis/dvp_delta_alignment_20260926/dv_02_probe_genes.py",
          "python analysis/dvp_delta_alignment_20260926/dv_03_align.py", "```", "",
          "| input (read-only) | sha256[:16] |", "| --- | --- |"]
    for rel in ["analysis/t4_ifn_landscape_20260925/hpa_ifn_landscape_wide.tsv",
                "analysis/t4_ifn_landscape_20260925/restriction_gradient.tsv",
                "analysis/protein_layer_scan_20260926/dvp_vs_delta_match.tsv",
                "analysis/protein_layer_scan_20260926/raw/hpa_dvp/dvp_cell_type.tsv.zip",
                "analysis/protein_layer_scan_20260926/raw/hpa_dvp/dvp_cell_type_group_data.tsv.zip"]:
        rp.append("| `%s` | `%s` |" % (rel, sha256(os.path.join(ROOT, rel.replace("/", os.sep)))))
    open(os.path.join(OUT, "repro.md"), "w", encoding="utf-8").write("\n".join(rp) + "\n")
    open(os.path.join(OUT, "dv_03_runlog.txt"), "w", encoding="utf-8").write("\n".join(LOG) + "\n")
    log("")
    log("done -> %s" % OUT)


if __name__ == "__main__":
    main()
