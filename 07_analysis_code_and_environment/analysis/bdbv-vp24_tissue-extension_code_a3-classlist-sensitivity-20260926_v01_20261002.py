"""How fragile is the A3 gate to the *choice* of the class list?

`report/A3_修订记录1_细胞类数_20260926.md` closed the "40 pre-registered vs 42 executed"
gap by deleting the two extra classes and re-running.  That answers the letter of the gap,
but not the sharper question: the pre-registration never contained the class table, so the
40-class set is a *reconstruction*, and the gate it feeds is a threshold comparison
(A2 requires Spearman >= 0.5, executed value 0.405).  A gate that sits 0.095 below its
threshold can be flipped by the class list alone.

This script therefore re-runs the A3 composition step under many class lists:

* ``--mode baseline``  the executed 42 classes (sanity: must reproduce the published values)
* ``--mode loo``       every 41-class list (drop one class) -- 42 fits
* ``--mode pairs``     every 40-class list (drop any two classes) -- C(42,2) = 861 fits

The deconvolution follows the published pipeline exactly: HPA single-cell nCPM column
profiles -> per-class share -> per-class top-200 markers (nCPM >= 10) -> per-sample NNLS on
per-sample-normalised bulk fractions -> proportions renormalised to sum 1.  Reported per
subset: the A2 Spearman (immune fraction vs the bulk immune-marker index, which does **not**
depend on the class list), the A1 top-1 calls, and whether the gate would pass.

Read-only w.r.t. every published artefact; writes only into
``analysis/a3_classlist_sensitivity_20260926/``.
"""

from __future__ import annotations

import argparse
import ast
import itertools
import json
import os
import time

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "a3_classlist_sensitivity_20260926")
PIPELINE = os.path.join(ROOT, "tools", "a3_deconvolution_20260926.py")
G8_TOOL = os.path.join(ROOT, "tools", "t4_g8_tissue_ifn_env_20260926.py")
HPA_SC = os.path.join(ROOT, "data", "hpa", "rna_single_cell_type.tsv")
ANCHORS = os.path.join(ROOT, "analysis", "a3_deconvolution_20260926", "anchors.json")
TISSUES = ["liver", "kidney", "brain", "adrenal gland", "lymph node", "skin", "sex organ", "lung"]

LOG: list[str] = []


def say(msg: str = "") -> None:
    LOG.append(msg)
    print(msg, flush=True)


def module_constants(path: str, names: tuple[str, ...]) -> dict:
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="baseline", choices=["baseline", "loo", "pairs"])
    ap.add_argument("--max-subsets", type=int, default=0, help="0 = no cap")
    ap.add_argument("--seed", type=int, default=20260926)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    consts = module_constants(PIPELINE, ("CLASSES", "IMMUNE_CLASSES", "IMMUNE_MARKERS",
                                         "ANCHOR_TOP1", "ANCHOR_MIN_HITS",
                                         "ANCHOR_MIN_RHO_IMMUNE", "TOP_MARKERS"))
    classes_all = list(consts["CLASSES"])
    immune = list(consts["IMMUNE_CLASSES"])
    anchor_top1 = {t: set(v) for t, v in consts["ANCHOR_TOP1"].items()}
    min_hits = int(consts["ANCHOR_MIN_HITS"])
    min_rho = float(consts["ANCHOR_MIN_RHO_IMMUNE"])
    top_markers = int(consts["TOP_MARKERS"])

    # immune-marker index per tissue: class-list independent, so taken from the published run
    a = json.load(open(ANCHORS, encoding="utf-8"))
    marker_index = {r["tissue"]: r["immune_marker_index"] for r in a["A2"]}
    say("immune-marker index reused (class-list independent): %s"
        % {t: round(marker_index[t], 4) for t in TISSUES})

    # ---- reference (same convention as the published pipeline)
    say("loading HPA single-cell reference ...")
    hpa = pd.read_csv(HPA_SC, sep="\t", usecols=["Gene name", "Cell type", "nCPM"])
    hpa = hpa.dropna(subset=["Gene name", "Cell type"])
    hpa["nCPM"] = pd.to_numeric(hpa["nCPM"], errors="coerce")
    wide = hpa.pivot_table(index="Gene name", columns="Cell type", values="nCPM",
                           aggfunc="mean", observed=True).fillna(0.0)
    say("HPA reference: %d genes x %d cell types" % wide.shape)

    def resolve(name: str):
        if name in wide.columns:
            return name
        if name == "müller glia":
            hits = [c for c in wide.columns if isinstance(c, str) and c.endswith("ller glia")]
            return hits[0] if hits else None
        return None

    def profiles_for(subset: list[str]) -> pd.DataFrame:
        cols = {}
        for cls in subset:
            resolved = [r for r in (resolve(m) for m in classes_all_members[cls]) if r]
            if resolved:
                cols[cls] = wide[resolved].mean(axis=1)
        return pd.DataFrame(cols)

    def markers_for(prof: pd.DataFrame) -> list[str]:
        share = prof.div(prof.sum(axis=1).replace(0, np.nan), axis=0)
        keep: set[str] = set()
        for cls in prof.columns:
            ok = prof[cls] >= 10.0
            s = share.loc[ok, cls].dropna().sort_values(ascending=False)
            keep.update(s.head(top_markers).index)
        return sorted(keep)

    classes_all_members = dict(consts["CLASSES"])

    # ---- bulk, via the published G8 loader (same gene mapping / filtering)
    say("loading bulk expression via the published G8 helper ...")
    import importlib.util
    spec = importlib.util.spec_from_file_location("g8_tool", G8_TOOL)
    g8 = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(g8)
    g8.log = lambda *a, **k: None
    meta = g8.load_sample_meta()
    z, depth, logcpm = g8.load_expression()
    meta = meta[meta["column"].isin(z.columns)].copy()
    cpm = (2.0 ** logcpm) - 1.0
    all8 = meta[meta["tissue"].isin(TISSUES)].copy()
    say("samples in the 8 primary tissues: %d" % len(all8))

    prof_full = profiles_for(classes_all)
    common = prof_full.index.intersection(cpm.index)
    say("genes shared with the macaque matrix: %d" % len(common))
    bulk_frac = cpm.loc[common].div(cpm.loc[common].sum(axis=0), axis=1)
    sample_cols = list(all8["column"])
    tissue_of = dict(zip(all8["column"], all8["tissue"]))

    def evaluate(subset: list[str]) -> dict:
        prof_all = profiles_for(subset)
        # marker selection happens on the FULL gene space (as published); only then are the
        # markers filtered down to the genes shared with the bulk matrix
        mk = [g for g in markers_for(prof_all) if g in common]
        prof = prof_all.loc[common]
        prof_frac = prof.div(prof.sum(axis=0), axis=1)
        A = prof_frac.loc[mk].to_numpy()
        per_sample = {}
        for col in sample_cols:
            b = bulk_frac.loc[mk, col].to_numpy()
            p, _ = nnls(A, b)
            tot = p.sum()
            per_sample[col] = p / tot if tot > 0 else p
        comp = pd.DataFrame([{**{c: float(v) for c, v in zip(prof_frac.columns, per_sample[c])},
                              "column": c, "tissue": tissue_of[c]} for c in sample_cols])
        tissue_comp = comp.groupby("tissue")[list(prof_frac.columns)].mean()
        imm = {t: float(tissue_comp.loc[t, [c for c in immune if c in tissue_comp.columns]].sum())
               for t in TISSUES}
        rho_a2 = float(spearmanr([imm[t] for t in TISSUES],
                                 [marker_index[t] for t in TISSUES]).statistic)
        top1 = {t: tissue_comp.loc[t].idxmax() for t in TISSUES}
        hits = sum(1 for t in TISSUES if top1[t] in anchor_top1[t])
        return {"n_classes": len(prof_frac.columns), "n_markers": len(mk),
                "A1_hits": hits, "A1_pass": bool(hits >= min_hits),
                "A2_rho": rho_a2, "A2_pass": bool(rho_a2 >= min_rho),
                "gate_pass": bool(hits >= min_hits and rho_a2 >= min_rho),
                "immune_fraction": imm, "top1": top1}

    # ---- subset plan
    if args.mode == "baseline":
        plan = [("baseline_42", classes_all)]
    elif args.mode == "loo":
        plan = [("drop:" + c, [x for x in classes_all if x != c]) for c in classes_all]
    else:
        pair = ["megakaryocytes_platelets", "erythrocytes"]
        plan = [("drop:" + ",".join(sorted(p)), [x for x in classes_all if x not in p])
                for p in itertools.combinations(classes_all, 2)]
        plan.sort(key=lambda kv: (kv[0] != "drop:" + ",".join(sorted(pair)), kv[0]))
    if args.max_subsets:
        plan = plan[:args.max_subsets]
    say("subsets to evaluate: %d (mode=%s)" % (len(plan), args.mode))

    rows = []
    t0 = time.time()
    for i, (label, subset) in enumerate(plan, 1):
        res = evaluate(subset)
        rows.append({"label": label, "dropped": ";".join(sorted(set(classes_all) - set(subset))),
                     **{k: v for k, v in res.items() if k not in ("immune_fraction", "top1")}})
        if i % 10 == 0 or i == len(plan):
            el = time.time() - t0
            say("  [%d/%d] %-46s A2=%+.3f A1=%d/8  elapsed %.1fs (%.2fs/subset)"
                % (i, len(plan), label[:46], res["A2_rho"], res["A1_hits"], el, el / i))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "subsets_%s.tsv" % args.mode), sep="\t", index=False)

    gate_pass = int(df["gate_pass"].sum())
    say("")
    say("=== summary (mode=%s, %d class lists) ===" % (args.mode, len(df)))
    say("A2 rho: min %+.3f  max %+.3f  median %+.3f"
        % (df["A2_rho"].min(), df["A2_rho"].max(), df["A2_rho"].median()))
    say("A2 alone would pass (>= %.2f) in %d/%d lists" % (min_rho, int(df["A2_pass"].sum()), len(df)))
    say("A1 hits: min %d max %d ; A1 pass in %d/%d"
        % (df["A1_hits"].min(), df["A1_hits"].max(), int(df["A1_pass"].sum()), len(df)))
    say("FULL GATE (A1 and A2) would pass in %d/%d class lists" % (gate_pass, len(df)))
    say("worst-case A2: %s" % df.loc[df["A2_rho"].idxmin(), ["label", "dropped", "A2_rho"]].to_dict())
    say("best-case A2:  %s" % df.loc[df["A2_rho"].idxmax(), ["label", "dropped", "A2_rho"]].to_dict())
    if args.mode == "pairs":
        pre = df[df["dropped"] == "erythrocytes;megakaryocytes_platelets"]
        if len(pre):
            say("pre-registered reading (drop both named classes): A2 = %+.3f, gate = %s"
                % (float(pre.iloc[0]["A2_rho"]), bool(pre.iloc[0]["gate_pass"])))
    say("wrote %s" % os.path.join(OUT, "subsets_%s.tsv" % args.mode))
    with open(os.path.join(OUT, "classlist_sensitivity_log.txt"), "a", encoding="utf-8") as fh:
        fh.write("\n".join(["# mode=%s %s" % (args.mode, time.strftime("%Y-%m-%dT%H:%M:%S"))] + LOG) + "\n")


if __name__ == "__main__":
    main()
