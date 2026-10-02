"""G7: add the PY-STAT1 nuclear-import (KPNA) axis to the frozen restriction gradient.

Pre-registration: report/G7_预注册_KPNA轴与乘数函数形式_20260926.md
Written before running; thresholds and quadrant rules are frozen there.

Inputs (all local, nothing downloaded):
    analysis/t4_ifn_landscape_20260925/hpa_ifn_landscape_wide.tsv   (154 cell types, nCPM)
    analysis/t4_ifn_landscape_20260925/restriction_gradient.tsv      (frozen Delta)

Outputs:
    analysis/g7_kpna_axis_20260926/kpna_axis.tsv
    analysis/g7_kpna_axis_20260926/g7_summary.json
    analysis/g7_kpna_axis_20260926/g7_report.txt
"""

import csv
import json
import math
import os
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAND = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")
OUT = os.path.join(ROOT, "analysis", "g7_kpna_axis_20260926")

# frozen axis definitions (pre-registration section 2.3)
AXES = OrderedDict([
    ("KPNA_axis", ["KPNA1", "KPNA5", "KPNA6"]),
    ("cargo_axis", ["STAT1", "STAT2", "IRF9"]),
])
MAX_AXIS_GENES = ["KPNA1", "KPNA5", "KPNA6"]

# Compartment labels for the two frozen lists (pre-registration section 2.6).
# REVISION 1 (2026-09-26, see report/G7_修订记录1_区室标签口径_20260926.md):
#   the pre-registered keyword rules (e.g. "intestin|colon") failed to capture the
#   actual HPA names for gut epithelium ("enterocytes", "colonocytes", "goblet cells",
#   ...).  Keyword matching was therefore replaced by EXPLICIT name lists so that the
#   membership test is auditable.  This touches the *labels* only: the axes, medians
#   and quadrant assignment are byte-identical to the pre-registered computation.
AMPLIFIED_RULES = OrderedDict([
    ("gut_epithelium", [
        "colonocytes", "enterocytes", "enteric stem cells",
        "enteric transient amplifying cells", "goblet cells", "paneth cells",
        "foveolar cells", "mucous neck cells", "gastric chief cells",
        "gastric progenitor cells", "parietal cells",
        "tuft cells",                      # ambiguous: gut and airway
    ]),
    ("urinary_renal_epithelium", [
        "urothelial cells", "podocytes", "proximal tubule cells",
        "distal convoluted tubule cells", "loop of henle epithelial cells",
        "renal collecting duct intercalated cells",
        "renal collecting duct principal cells", "renal connecting tubule cells",
        "papillary tip epithelial cells",
    ]),
    ("hepatic", [
        "hepatocytes", "cholangiocytes", "hepatic stellate cells", "kupffer cells",
    ]),
    ("airway_epithelium", [
        "alveolar cells type 1", "alveolar cells type 2", "transitional alveolar cells",
        "respiratory basal cells", "respiratory ciliated cells",
        "respiratory deuterosomal cells", "respiratory ionocytes",
        "respiratory secretory cells", "submucosal glandular cells",
    ]),
    ("myeloid", [
        "macrophages", "monocytes", "monocyte progenitors", "neutrophils",
        "neutrophil progenitors", "mast cells", "microglia", "kupffer cells",
        "hofbauer cells", "cdc", "pdcs",
    ]),
])
SANCTUARY_RULES = OrderedDict([
    ("retina_PPE_and_neurons", [
        "retinal pigment epithelial cells", "rod photoreceptor cells",
        "cone photoreceptor cells", "retinal amacrine cells",
        "retinal bipolar cells", "retinal ganglion cells",
        "retinal horizontal cells", "müller glia",
    ]),
    ("retina_PPE_only", ["retinal pigment epithelial cells"]),
    ("photoreceptors_only", ["rod photoreceptor cells", "cone photoreceptor cells"]),
    ("testis_somatic", ["sertoli cells", "leydig cells", "peritubular myoid cells"]),
    ("testis_germline", ["undifferentiated spermatogonia", "differentiating spermatogonia",
                         "early primary spermatocytes", "late primary spermatocytes",
                         "early spermatids", "late spermatids", "oocytes"]),
    ("choroid_plexus", ["choroid plexus epithelial cells"]),
    ("ependymal", ["ependymal cells"]),
])


def log10p1(raw):
    try:
        return math.log10(float(raw) + 1.0)
    except (TypeError, ValueError):
        return None


def zscores(values):
    vals = [v for v in values.values() if v is not None]
    if not vals:
        return {}, None, None
    mu = sum(vals) / len(vals)
    sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals))
    if sd == 0:
        return {k: 0.0 for k in values}, mu, 0.0
    return {k: (None if v is None else (v - mu) / sd) for k, v in values.items()}, mu, sd


def _ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(xs, ys):
    rx, ry = _ranks(xs), _ranks(ys)
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float("nan")


def quantile(sorted_vals, q):
    if not sorted_vals:
        return None
    idx = q * (len(sorted_vals) - 1)
    lo, hi = int(math.floor(idx)), int(math.ceil(idx))
    if lo == hi:
        return sorted_vals[lo]
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (idx - lo)


def matches(name, rules):
    """Explicit-membership test (revision 1): exact cell-type name in the frozen list."""
    return [label for label, names in rules.items() if name in names]


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def say(msg=""):
        lines.append(msg)
        print(msg)

    wide = list(csv.DictReader(open(os.path.join(LAND, "hpa_ifn_landscape_wide.tsv"),
                                   encoding="utf-8"), delimiter="\t"))
    grad = {r["cell_type"]: r for r in csv.DictReader(
        open(os.path.join(LAND, "restriction_gradient.tsv"), encoding="utf-8"), delimiter="\t")}
    cell_types = [r["cell_type"] for r in wide]
    missing_delta = [c for c in cell_types if c not in grad]
    say("=" * 78)
    say("G7  KPNA axis added to the frozen restriction gradient")
    say("pre-registration: report/G7_预注册_KPNA轴与乘数函数形式_20260926.md")
    say("=" * 78)
    say("cell types: %d   matched to Delta: %d   unmatched: %s"
        % (len(cell_types), len(cell_types) - len(missing_delta), missing_delta))

    gene_z, gene_raw, gene_stats = {}, {}, {}
    all_genes = sorted({g for gs in list(AXES.values()) + [MAX_AXIS_GENES] for g in gs})
    for g in all_genes:
        raw = {r["cell_type"]: log10p1(r.get(g, "")) for r in wide}
        zz, mu, sd = zscores(raw)
        gene_z[g], gene_raw[g] = zz, raw
        gene_stats[g] = {"mean_log10_nCPM": mu, "sd_log10_nCPM": sd,
                         "n_present": sum(1 for v in raw.values() if v is not None)}

    axis = {}
    for name, genes in AXES.items():
        axis[name] = {
            c: sum(gene_z[g][c] for g in genes if gene_z[g].get(c) is not None)
               / max(1, sum(1 for g in genes if gene_z[g].get(c) is not None))
            for c in cell_types
        }
    axis["kpna_max"] = {
        c: max((gene_z[g][c] for g in MAX_AXIS_GENES if gene_z[g].get(c) is not None),
               default=None)
        for c in cell_types
    }

    delta_z_raw, _, _ = zscores({c: (float(grad[c]["delta_restriction"]) if c in grad else None)
                                 for c in cell_types})

    # T1/T2/T3 pre-checks
    common = [c for c in cell_types if c in grad and axis["KPNA_axis"][c] is not None]
    d_common = [float(grad[c]["delta_restriction"]) for c in common]
    k_common = [axis["KPNA_axis"][c] for c in common]
    rho_delta = spearman(d_common, k_common)
    rho_pairs = {}
    for i, g1 in enumerate(MAX_AXIS_GENES):
        for g2 in MAX_AXIS_GENES[i + 1:]:
            pair = [c for c in cell_types if gene_z[g1][c] is not None and gene_z[g2][c] is not None]
            rho_pairs["%s~%s" % (g1, g2)] = spearman([gene_z[g1][c] for c in pair],
                                                     [gene_z[g2][c] for c in pair])
    kv = sorted(k_common)
    iqr = quantile(kv, 0.75) - quantile(kv, 0.25)

    say("")
    say("[T1] Spearman rho(Delta, KPNA_axis) = %+.3f -> %s"
        % (rho_delta, "PASS (|rho|<0.6: separable)" if abs(rho_delta) < 0.6
           else "FAIL (|rho|>=0.6: near-synonymous axes)"))
    say("[T2] pairwise rho within KPNA1/5/6:")
    for k, v in rho_pairs.items():
        say("      %s: %+.3f" % (k, v))
    say("[T3] KPNA_axis IQR over %d cell types = %.3f -> %s"
        % (len(kv), iqr, "PASS" if iqr >= 0.5 else "FAIL (near-constant)"))

    # additional reported correlations (layers B and C; no judgment attached)
    import_axis = {c: (axis["KPNA_axis"][c] + axis["cargo_axis"][c]) / 2.0 for c in cell_types}
    extra_rho = {
        "rho_delta_vs_cargo_axis": spearman(d_common, [axis["cargo_axis"][c] for c in common]),
        "rho_delta_vs_kpna_max": spearman(d_common, [axis["kpna_max"][c] for c in common]),
        "rho_delta_vs_import_axis": spearman(d_common, [import_axis[c] for c in common]),
    }
    for k, v in extra_rho.items():
        say("[extra] %s = %+.3f" % (k, v))

    # quadrants
    mD = quantile(sorted(d_common), 0.5)
    mK = quantile(kv, 0.5)
    say("")
    say("[quadrants] median Delta = %.3f   median KPNA_axis = %.3f" % (mD, mK))

    records = []
    for c in cell_types:
        if c not in grad:
            continue
        d = float(grad[c]["delta_restriction"])
        k = axis["KPNA_axis"][c]
        if k is None:
            q = "NA"
        elif d >= mD and k >= mK:
            q = "Q1_highD_highK"
        elif d >= mD and k < mK:
            q = "Q2_highD_lowK"
        elif d < mD and k >= mK:
            q = "Q3_lowD_highK"
        else:
            q = "Q4_lowD_lowK"
        records.append({
            "cell_type": c,
            "delta_restriction": d,
            "delta_rank": int(grad[c]["rank"]),
            "KPNA1_z": gene_z["KPNA1"][c],
            "KPNA5_z": gene_z["KPNA5"][c],
            "KPNA6_z": gene_z["KPNA6"][c],
            "KPNA1_nCPM": 10 ** gene_raw["KPNA1"][c] - 1 if gene_raw["KPNA1"][c] is not None else None,
            "KPNA5_nCPM": 10 ** gene_raw["KPNA5"][c] - 1 if gene_raw["KPNA5"][c] is not None else None,
            "KPNA6_nCPM": 10 ** gene_raw["KPNA6"][c] - 1 if gene_raw["KPNA6"][c] is not None else None,
            "KPNA_axis": k,
            "kpna_max": axis["kpna_max"][c],
            "cargo_axis": axis["cargo_axis"][c],
            "quadrant": q,
            "amplified_groups": ";".join(matches(c, AMPLIFIED_RULES)),
            "sanctuary_groups": ";".join(matches(c, SANCTUARY_RULES)),
        })

    counts = {}
    for r in records:
        counts[r["quadrant"]] = counts.get(r["quadrant"], 0) + 1
    say("[quadrant sizes] " + ", ".join("%s=%d" % (k, v) for k, v in sorted(counts.items())))

    q2 = sorted([r for r in records if r["quadrant"] == "Q2_highD_lowK"],
                key=lambda r: -r["delta_restriction"])
    say("")
    say("[Q2] high-Delta / low-KPNA cell types (n=%d) -- the only set whose predicted"
        " direction flips if the KPNA axis is included" % len(q2))
    say("     %-46s %7s %7s %7s  %s" % ("cell_type", "Delta", "KPNAax", "Drank", "overlap"))
    for r in q2:
        ov = ",".join(x for x in [r["amplified_groups"], r["sanctuary_groups"]] if x) or "-"
        say("     %-46s %7.2f %7.2f %7d  %s"
            % (r["cell_type"], r["delta_restriction"], r["KPNA_axis"], r["delta_rank"], ov))

    say("")
    say("[overlap with frozen compartment lists]")
    comp_hits = {}
    for label, rules, key in (("amplified", AMPLIFIED_RULES, "amplified_groups"),
                              ("sanctuary", SANCTUARY_RULES, "sanctuary_groups")):
        comp_hits[label] = {}
        for grp in rules:
            hits = [r for r in records if grp in r[key].split(";")]
            comp_hits[label][grp] = [r["cell_type"] for r in hits]
            by_q = {}
            for h in hits:
                by_q[h["quadrant"]] = by_q.get(h["quadrant"], 0) + 1
            say("   %-9s %-18s n=%2d  %s"
                % (label, grp, len(hits),
                   ", ".join("%s:%d" % (k, v) for k, v in sorted(by_q.items())) or "-"))
            say("        " + ("; ".join("%s(Δrank %d, %s)" % (h["cell_type"], h["delta_rank"],
                                                              h["quadrant"].split("_")[0])
                                        for h in sorted(hits, key=lambda r: r["delta_rank"]))
                              if hits else "-"))

    say("")
    say("[axis extremes] KPNA_axis Top-15 / Bottom-15")
    srt = sorted([r for r in records if r["KPNA_axis"] is not None],
                 key=lambda r: -r["KPNA_axis"])
    for r in srt[:15]:
        say("     %-46s KPNAax=%6.2f  KPNA1=%7.1f KPNA5=%7.1f KPNA6=%7.1f"
            % (r["cell_type"], r["KPNA_axis"], r["KPNA1_nCPM"], r["KPNA5_nCPM"], r["KPNA6_nCPM"]))
    say("     " + "-" * 72)
    for r in srt[-15:]:
        say("     %-46s KPNAax=%6.2f  KPNA1=%7.1f KPNA5=%7.1f KPNA6=%7.1f"
            % (r["cell_type"], r["KPNA_axis"], r["KPNA1_nCPM"], r["KPNA5_nCPM"], r["KPNA6_nCPM"]))

    # sensitivity: rank drift (secondary only)
    say("")
    say("[sensitivity] rank drift if the KPNA axis is folded in (NOT a replacement for Delta)")
    # Pre-registration 2.4 specifies z(Delta) + z(KPNA_axis): BOTH terms must be
    # standardised before combining.  The KPNA axis is a mean of three z-scores
    # and therefore has SD < 1, so using it raw silently down-weights it.
    # (Corrected 2026-09-26; see report/G7_修订记录2_敏感性口径_20260926.md.)
    kpna_z, _, _ = zscores({c: axis["KPNA_axis"][c] for c in cell_types})
    kmax_z, _, _ = zscores({c: axis["kpna_max"][c] for c in cell_types})
    adj = {
        "Delta_adj_equal": {c: delta_z_raw[c] + kpna_z[c] for c in common},
        "Delta_adj_min": {c: min(delta_z_raw[c], kpna_z[c]) for c in common},
        "Delta_adj_maxK": {c: delta_z_raw[c] + kmax_z[c] for c in common},
    }
    base_order = list(common)
    base_sorted = sorted(base_order, key=lambda c: -float(grad[c]["delta_restriction"]))
    base_rank = {c: i + 1 for i, c in enumerate(base_sorted)}
    base_top20 = set(base_sorted[:20])
    sens = {}
    for name, vals in adj.items():
        order = sorted(base_order, key=lambda c: -vals[c])
        new_rank = {c: i + 1 for i, c in enumerate(order)}
        rho = spearman([float(grad[c]["delta_restriction"]) for c in base_order],
                       [vals[c] for c in base_order])
        moved = sorted(base_order, key=lambda c: -abs(base_rank[c] - new_rank[c]))[:8]
        sens[name] = {"spearman_vs_delta": rho,
                      "top20_overlap": len(base_top20 & set(order[:20])),
                      "largest_moves": [{"cell_type": c, "delta_rank": base_rank[c],
                                         "new_rank": new_rank[c]} for c in moved]}
        say("   %-16s rho(Delta)=%+.3f  top20 overlap=%d/20"
            % (name, rho, len(base_top20 & set(order[:20]))))
        say("        largest moves: " + "; ".join(
            "%s %d->%d" % (m["cell_type"], m["delta_rank"], m["new_rank"])
            for m in sens[name]["largest_moves"]))

    # write outputs
    fields = list(records[0].keys())
    with open(os.path.join(OUT, "kpna_axis.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(records)

    summary = {
        "pre_registration": "report/G7_预注册_KPNA轴与乘数函数形式_20260926.md",
        "inputs": {
            "expression": "analysis/t4_ifn_landscape_20260925/hpa_ifn_landscape_wide.tsv",
            "delta": "analysis/t4_ifn_landscape_20260925/restriction_gradient.tsv",
        },
        "axes": dict(AXES),
        "gene_stats": gene_stats,
        "n_cell_types": len(records),
        "prechecks": {
            "T1_rho_delta_kpna": rho_delta,
            "T1_pass": bool(abs(rho_delta) < 0.6),
            "T2_pairwise_rho": rho_pairs,
            "T3_kpna_axis_iqr": iqr,
            "T3_pass": bool(iqr >= 0.5),
            "extra_rho": extra_rho,
        },
        "medians": {"delta": mD, "kpna_axis": mK},
        "quadrant_sizes": counts,
        "Q2_highD_lowK": q2,
        "compartment_hits": comp_hits,
        "sensitivity_rank_drift": sens,
    }
    with open(os.path.join(OUT, "g7_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(OUT, "g7_report.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    say("")
    say("wrote %s\\{kpna_axis.tsv, g7_summary.json, g7_report.txt}"
        % os.path.relpath(OUT, ROOT))


if __name__ == "__main__":
    main()
