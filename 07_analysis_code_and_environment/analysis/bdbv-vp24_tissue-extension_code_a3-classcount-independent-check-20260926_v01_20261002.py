"""Independent check of the A3 class-count correction and its 40-class sensitivity.

This does **not** re-read the sibling report's summary numbers.  Everything below is
recomputed from primary artefacts on disk:

* the two tissue tables (42 executed / 40 pre-registered) -- correlations are recomputed,
  including **exact** two-sided permutation p-values by enumerating all 8! orderings
  rather than trusting the seeded sampled values in the summaries;
* the immune fraction is rebuilt from ``tissue_composition.tsv`` by summing the
  ``IMMUNE_CLASSES`` columns parsed out of the pipeline module (AST), then Spearman-correlated
  with the independently measured immune-marker index;
* the anchor A1 top-1 calls are recomputed as per-tissue argmax of the composition table;
* the claim "exactly one leaf value of ``a3_summary.json`` differs from the pre-fix snapshot"
  is checked by a recursive leaf diff, and the other five files are checked byte-for-byte;
* ``CLASSES`` is parsed out of the pipeline module to confirm the 42 -> 40 removal.

Read-only w.r.t. every published artefact; writes only into
``analysis/a3_classcount40_independent_check_20260926/``.
"""

from __future__ import annotations

import ast
import csv
import json
import os
from itertools import permutations

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FORMAL = os.path.join(ROOT, "analysis", "a3_deconvolution_20260926")
SENS = os.path.join(ROOT, "analysis", "a3_classcount40_sensitivity_20260926")
SNAP = os.path.join(ROOT, "analysis", "_audit_tmp", "a3_pre_classcount_fix_20260926")
OUT = os.path.join(ROOT, "analysis", "a3_classcount40_independent_check_20260926")
PIPELINE = os.path.join(ROOT, "tools", "a3_deconvolution_20260926.py")
PREREG = os.path.join(ROOT, "report", "A3_预注册_组织构成去卷积_20260926.md")

TISSUES = ["liver", "kidney", "brain", "adrenal gland", "lymph node", "skin", "sex organ", "lung"]

LOG: list[str] = []


def say(msg: str = "") -> None:
    LOG.append(msg)
    print(msg, flush=True)


def read_table(path: str) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def module_constants(path: str, names: tuple[str, ...]) -> dict:
    """Pull literal constants out of a module without importing it."""
    tree = ast.parse(open(path, encoding="utf-8").read())
    found: dict = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id in names:
                    found[tgt.id] = ast.literal_eval(node.value)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            if node.target.id in names and node.value is not None:
                found[node.target.id] = ast.literal_eval(node.value)
    return found


def exact_spearman_perm_p(x: list[float], y: list[float], observed: float) -> float:
    """Two-sided exact permutation p over all n! re-orderings (n = 8 -> 40320)."""
    n = len(x)
    ry = stats.rankdata(y)
    rx = stats.rankdata(x)
    hits = 0
    total = 0
    for perm in permutations(range(n)):
        rho = np.corrcoef(rx, ry[list(perm)])[0, 1]
        if abs(rho) >= abs(observed) - 1e-12:
            hits += 1
        total += 1
    return hits / total, total


def leaf_diff(a, b, path="") -> list[tuple[str, object, object]]:
    """Recursive leaf-level difference list ('' when identical)."""
    diff: list[tuple[str, object, object]] = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                diff.append((path + "/" + str(k), a.get(k, "<absent>"), b.get(k, "<absent>")))
            else:
                diff.extend(leaf_diff(a[k], b[k], path + "/" + str(k)))
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            diff.append((path + "/<len>", len(a), len(b)))
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                diff.extend(leaf_diff(x, y, path + "/%d" % i))
    else:
        if a != b:
            diff.append((path, a, b))
    return diff


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    report: dict = {"role": "independent recomputation of the A3 class-count correction"}

    consts = module_constants(PIPELINE, ("CLASSES", "IMMUNE_CLASSES"))
    classes = list(consts["CLASSES"])
    immune = set(consts["IMMUNE_CLASSES"])
    say("=== 1. class lists ===")
    say("CLASSES            n=%d" % len(classes))
    say("IMMUNE_CLASSES     n=%d -> %s" % (len(immune), sorted(immune)))
    extra = [c for c in ("megakaryocytes_platelets", "erythrocytes") if c in classes]
    say("extra classes named by the correction: %s (both present: %s)" % (extra, len(extra) == 2))
    say("neither extra class is in IMMUNE_CLASSES: %s"
        % all(c not in immune for c in extra))
    report["classes"] = {"n_executed": len(classes), "n_immune": len(immune),
                         "extra": extra, "extra_in_immune": sorted(set(extra) & immune)}

    say("")
    say("=== 2. composition tables ===")
    comp_f = read_table(os.path.join(FORMAL, "tissue_composition.tsv"))
    comp_s = read_table(os.path.join(SENS, "tissue_composition.tsv"))
    cols_f = [c for c in comp_f[0] if c != "tissue"]
    cols_s = [c for c in comp_s[0] if c != "tissue"]
    say("42-run classes in table: %d ; 40-run: %d ; difference: %s"
        % (len(cols_f), len(cols_s), sorted(set(cols_f) - set(cols_s))))
    say("40-run set == executed set minus the two named classes: %s"
        % (sorted(cols_s) == sorted(set(cols_f) - set(extra))))
    report["composition_tables"] = {
        "n_classes_formal": len(cols_f), "n_classes_sens": len(cols_s),
        "dropped_in_table": sorted(set(cols_f) - set(cols_s)),
        "identical_to_named_removal": sorted(cols_s) == sorted(set(cols_f) - set(extra)),
    }

    def immune_from_table(rows: list[dict], cols: list[str]) -> dict:
        out = {}
        for r in rows:
            out[r["tissue"]] = sum(float(r[c]) for c in cols if c in immune)
        return out

    imm_f = immune_from_table(comp_f, cols_f)
    imm_s = immune_from_table(comp_s, cols_s)
    say("immune fraction recomputed by summing IMMUNE_CLASSES columns (42-run):")
    for t in TISSUES:
        say("   %-14s %.6f" % (t, imm_f[t]))

    say("")
    say("=== 3. A2 anchor recomputed from the tissue tables ===")
    anchors_f = json.load(open(os.path.join(FORMAL, "anchors.json"), encoding="utf-8"))
    anchors_s = json.load(open(os.path.join(SENS, "anchors.json"), encoding="utf-8"))
    idx_f = {r["tissue"]: r["immune_marker_index"] for r in anchors_f["A2"]}
    idx_s = {r["tissue"]: r["immune_marker_index"] for r in anchors_s["A2"]}
    say("immune-marker index (independent variable) identical across runs: %s"
        % all(abs(idx_f[t] - idx_s[t]) < 1e-12 for t in TISSUES))
    rho_f = stats.spearmanr([imm_f[t] for t in TISSUES], [idx_f[t] for t in TISSUES])
    rho_s = stats.spearmanr([imm_s[t] for t in TISSUES], [idx_s[t] for t in TISSUES])
    say("A2 rho recomputed: 42-run %.4f (reported %.3f) | 40-run %.4f (reported %.3f)"
        % (rho_f.statistic, anchors_f["anchors"]["A2_rho_immune"],
           rho_s.statistic, anchors_s["anchors"]["A2_rho_immune"]))
    say("   -> recomputed A2 drift |Delta| = %.4f (correction reports 0.143)"
        % abs(rho_f.statistic - rho_s.statistic))
    per_tissue_drift = {t: imm_s[t] - imm_f[t] for t in sorted(TISSUES, key=lambda t: abs(imm_s[t] - imm_f[t]), reverse=True)}
    say("largest per-tissue immune-fraction moves (40 minus 42): %s"
        % ", ".join("%s %+.5f" % (t, v) for t, v in list(per_tissue_drift.items())[:4]))
    report["A2_anchor"] = {
        "recomputed_rho_formal": rho_f.statistic, "reported_rho_formal": anchors_f["anchors"]["A2_rho_immune"],
        "recomputed_rho_sens": rho_s.statistic, "reported_rho_sens": anchors_s["anchors"]["A2_rho_immune"],
        "recomputed_drift": abs(rho_f.statistic - rho_s.statistic),
        "immune_fraction_formal": imm_f, "immune_fraction_sens": imm_s,
        "per_tissue_drift": per_tissue_drift,
        "both_below_threshold": bool(rho_f.statistic < 0.5 and rho_s.statistic < 0.5),
    }

    say("")
    say("=== 4. A1 anchor recomputed as per-tissue argmax ===")
    def top1(rows, cols):
        return {r["tissue"]: max(cols, key=lambda c: float(r[c])) for r in rows}
    t1_f, t1_s = top1(comp_f, cols_f), top1(comp_s, cols_s)
    same = all(t1_f[t] == t1_s[t] for t in TISSUES)
    say("per-tissue top-1 identical across runs: %s" % same)
    for t in TISSUES:
        say("   %-14s 42-run %-22s 40-run %-22s %s"
            % (t, t1_f[t], t1_s[t], "same" if t1_f[t] == t1_s[t] else "DIFFERS"))
    report["A1_top1"] = {"identical": same, "formal": t1_f, "sens": t1_s}

    say("")
    say("=== 5. correlations recomputed from the tissue tables + exact permutation p ===")
    tf = {r["tissue"]: r for r in read_table(os.path.join(FORMAL, "a3_tissue_table.tsv"))}
    ts = {r["tissue"]: r for r in read_table(os.path.join(SENS, "a3_tissue_table.tsv"))}
    checks = {}
    for label, key_x, key_y in (("primary delta_weighted vs E_adj", "delta_weighted", "E_adj"),
                                ("reproduction delta_median_g8 vs E_raw", "delta_median_g8", "E_raw")):
        row = {}
        for tag, tab, want in (("formal42", tf, -0.524), ("sens40", ts, -0.5)):
            x = [float(tab[t][key_x]) for t in TISSUES]
            y = [float(tab[t][key_y]) for t in TISSUES]
            rho = stats.spearmanr(x, y).statistic
            p_exact, n_perm = exact_spearman_perm_p(x, y, rho)
            row[tag] = {"rho": rho, "p_exact_perm": p_exact, "n_perm": n_perm}
            say("   %-38s %-9s rho = %+.4f   exact p = %.4f (%d permutations)"
                % (label, tag, rho, p_exact, n_perm))
        checks[label] = row
    say("   frozen G8 quantity delta_median_g8 identical across runs: %s"
        % all(abs(float(tf[t]["delta_median_g8"]) - float(ts[t]["delta_median_g8"])) < 1e-12 for t in TISSUES))
    report["correlations"] = checks
    report["frozen_delta_untouched"] = all(
        abs(float(tf[t]["delta_median_g8"]) - float(ts[t]["delta_median_g8"])) < 1e-12 for t in TISSUES)

    say("")
    say("=== 6. 'exactly one leaf differs' claim vs the pre-fix snapshot ===")
    cur = json.load(open(os.path.join(FORMAL, "a3_summary.json"), encoding="utf-8"))
    old = json.load(open(os.path.join(SNAP, "a3_summary.json"), encoding="utf-8"))
    diff = leaf_diff(old, cur)
    say("leaf differences: %d" % len(diff))
    for path, a, b in diff:
        say("   %s" % path)
        say("     old: %s" % str(a)[:130])
        say("     new: %s" % str(b)[:130])
    report["summary_json_leaf_diff"] = [
        {"path": p, "old": str(a)[:400], "new": str(b)[:400]} for p, a, b in diff]

    say("")
    say("=== 7. other five artefacts byte-identical to the snapshot ===")
    import hashlib
    byte_checks = {}
    for fn in ("anchors.json", "composition.tsv", "tissue_composition.tsv",
               "a3_tissue_table.tsv", "a3_correlations.tsv"):
        pf, ps = os.path.join(FORMAL, fn), os.path.join(SNAP, fn)
        if not (os.path.exists(pf) and os.path.exists(ps)):
            byte_checks[fn] = "absent"
            say("   %-26s ABSENT in one of the trees" % fn)
            continue
        hf = hashlib.sha256(open(pf, "rb").read()).hexdigest()
        hs = hashlib.sha256(open(ps, "rb").read()).hexdigest()
        byte_checks[fn] = (hf == hs)
        say("   %-26s identical: %-5s  sha256:%s" % (fn, hf == hs, hf[:16]))
    report["snapshot_byte_checks"] = {k: (v if isinstance(v, str) else bool(v)) for k, v in byte_checks.items()}

    say("")
    say("=== 8. pre-registration text: does it really say '40' and carry no class table? ===")
    txt = open(PREREG, encoding="utf-8").read()
    hits40 = txt.count("40 类") + txt.count("40 pre-registered")
    has_md_table = "\n| " in txt
    say("occurrences of '40 类' / '40 pre-registered': %d" % hits40)
    say("markdown tables present in the pre-registration: %s" % has_md_table)
    say("mentions of the two extra classes: %s"
        % {c: (c in txt) for c in extra})
    report["preregistration_text"] = {
        "count_40": hits40, "has_markdown_table": has_md_table,
        "mentions_extra_classes": {c: (c in txt) for c in extra},
        "count_hits_ge_6": txt.count("≥6") + txt.count(">=6"),
    }

    with open(os.path.join(OUT, "independent_check.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(OUT, "independent_check_log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    say("")
    say("wrote %s" % os.path.join(OUT, "independent_check.json"))


if __name__ == "__main__":
    main()
