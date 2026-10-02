#!/usr/bin/env python
"""A3 class-count sensitivity: does the 40-vs-42 class discrepancy change anything?

Background (registration gap, logged in report/手稿v20_协调者复核与A类收尾_20260926.md §3.1)
--------------------------------------------------------------------------------------
The A3 pre-registration (report/A3_预注册_组织构成去卷积_20260926.md §2.1) fixes the
cell-class table at **40 classes** and says "see the table below", but that table is not
present in the pre-registration file. The executed script
(tools/a3_deconvolution_20260926.py) defines **42** classes -- it adds
`megakaryocytes_platelets` and `erythrocytes`. The formal run therefore carries an
undeclared deviation from its own pre-registration. The coordinator did not re-run the
40-class version, so "does the two-class difference affect any gate or verdict?" was
recorded as **unverified**.

This script closes that gap. It imports the executed A3 module unchanged, removes exactly
the two extra classes, and re-runs the whole pipeline with a different output directory.
Nothing in the frozen Delta values, the G8 module, or the anchor definitions is touched.

Decision rule (fixed here, before the numbers are read)
------------------------------------------------------
The 40-class run is judged **judgement-preserving** iff all four of these hold against the
42-class formal run (analysis/a3_deconvolution_20260926/a3_summary.json):
  R1 pipeline reproduction still passes (rho(delta_median, E_raw) ~ -0.333, n = 8);
  R2 anchor A1 top-1 hits: same count;
  R3 anchor A2 rho(immune fraction, immune-marker index): same pass/fail, and |drho| <= 0.05;
  R4 primary rho(delta_weighted, E_adj): same sign, and |drho| <= 0.05;
  R5 the verdict string is identical.
If any fails, the discrepancy must be reported as *judgement-relevant* and the manuscript
must carry the 42-class number with an explicit deviation note.

Outputs (analysis/a3_classcount40_sensitivity_20260926/)
    comp40_*.json / *.tsv (full pipeline outputs, names unchanged inside the dir) and
    classcount40_comparison.json
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A3_TOOL = ROOT / "tools" / "a3_deconvolution_20260926.py"
FORMAL = ROOT / "analysis" / "a3_deconvolution_20260926" / "a3_summary.json"
OUT = ROOT / "analysis" / "a3_classcount40_sensitivity_20260926"
OUT.mkdir(parents=True, exist_ok=True)

# the two classes the executed script added on top of the pre-registered 40
DROPPED = ["megakaryocytes_platelets", "erythrocytes"]


def load_a3():
    spec = importlib.util.spec_from_file_location("a3_tool", A3_TOOL)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    formal = json.loads(FORMAL.read_text(encoding="utf-8"))

    a3 = load_a3()
    n_before = len(a3.CLASSES)
    reduced = {k: v for k, v in a3.CLASSES.items() if k not in DROPPED}
    assert n_before - len(reduced) == 2, (n_before, len(reduced))
    a3.CLASSES = reduced                     # pre-registered class set
    a3.OUT = OUT                             # keep the formal run untouched

    print(f"class-count sensitivity: {n_before} classes (executed) -> "
          f"{len(reduced)} classes (pre-registered), dropped = {DROPPED}")
    print(f"output dir: {OUT}")
    a3.main()

    new = json.loads((OUT / "a3_summary.json").read_text(encoding="utf-8"))

    def a2(tbl):
        return {r["tissue"]: r["manual_immune_fraction"] for r in tbl}

    rows = []
    for key, name, fmt in [
        ("n_classes", "n classes", "{}"),
        ("n_marker_genes", "n marker genes", "{}"),
        ("composition_R2_on_ISG", "R^2 (composition -> ISG_core)", "{:.3f}"),
    ]:
        rows.append({"quantity": name,
                     "formal_42": fmt.format(formal[key]),
                     "sensitivity_40": fmt.format(new[key])})
    for label, f, n in [
        ("anchor A1 top-1 hits", formal["anchors"]["A1_top1_hits"], new["anchors"]["A1_top1_hits"]),
        ("anchor A1 pass", formal["anchors"]["A1_pass"], new["anchors"]["A1_pass"]),
        ("anchor A2 rho (immune vs marker index)",
         formal["anchors"]["A2_rho_immune"], new["anchors"]["A2_rho_immune"]),
        ("anchor A2 pass", formal["anchors"]["A2_pass"], new["anchors"]["A2_pass"]),
        ("all anchors pass", formal["anchors"]["all_pass"], new["anchors"]["all_pass"]),
        ("primary rho (delta_weighted vs E_adj)",
         formal["correlations"][0]["rho"], new["correlations"][0]["rho"]),
        ("primary p_perm", formal["correlations"][0]["p_perm"], new["correlations"][0]["p_perm"]),
        ("pipeline reproduction rho",
         formal["pipeline_reproduction"]["rho"], new["pipeline_reproduction"]["rho"]),
        ("pipeline reproduction ok", formal["reproduction_ok"], new["reproduction_ok"]),
        ("verdict", formal["verdict"], new["verdict"]),
    ]:
        rows.append({"quantity": label, "formal_42": str(f), "sensitivity_40": str(n)})

    cmp_tbl = pd.DataFrame(rows)
    cmp_tbl.to_csv(OUT / "classcount40_comparison.tsv", sep="\t", index=False)
    print("\n== 42 -> 40 class comparison ==")
    print(cmp_tbl.to_string(index=False))

    rho_delta_a2 = abs(float(new["anchors"]["A2_rho_immune"]) - float(formal["anchors"]["A2_rho_immune"]))
    rho_delta_p = abs(float(new["correlations"][0]["rho"]) - float(formal["correlations"][0]["rho"]))
    checks = {
        "R1_reproduction_ok": bool(new["reproduction_ok"]),
        "R2_A1_hits_equal": new["anchors"]["A1_top1_hits"] == formal["anchors"]["A1_top1_hits"],
        "R3_A2_same_pass_and_small_drift": bool(
            new["anchors"]["A2_pass"] == formal["anchors"]["A2_pass"] and rho_delta_a2 <= 0.05),
        "R4_primary_same_sign_small_drift": bool(
            (new["correlations"][0]["rho"] is not None
             and float(new["correlations"][0]["rho"]) * float(formal["correlations"][0]["rho"]) > 0
             and rho_delta_p <= 0.05)),
        "R5_verdict_identical": new["verdict"] == formal["verdict"],
    }
    checks["judgement_preserving"] = all(checks.values())

    drift = {
        "A2_rho_abs_drift": round(rho_delta_a2, 4),
        "primary_rho_abs_drift": round(rho_delta_p, 4),
        "A1_top1_by_tissue": {
            "formal_42": {r["tissue"]: r["top1"] for r in formal["anchor_A1_table"]},
            "sensitivity_40": {r["tissue"]: r["top1"] for r in new["anchor_A1_table"]},
        },
        "A2_immune_fraction_by_tissue": {
            "formal_42": a2(formal["anchor_A2_table"]),
            "sensitivity_40": a2(new["anchor_A2_table"]),
        },
    }
    payload = {
        "role": "class-count sensitivity for the A3 registration gap (40 pre-registered vs 42 executed)",
        "formal_run": str(FORMAL.relative_to(ROOT)).replace("\\", "/"),
        "dropped_classes": DROPPED,
        "checks": checks,
        "drift": drift,
        "comparison_table": cmp_tbl.to_dict(orient="records"),
    }
    (OUT / "classcount40_comparison.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print("\n== decision rule ==")
    print(json.dumps(checks, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
