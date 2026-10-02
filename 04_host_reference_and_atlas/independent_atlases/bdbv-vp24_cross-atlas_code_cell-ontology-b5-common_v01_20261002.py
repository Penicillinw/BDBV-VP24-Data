"""Shared helpers for the B5 harmonised-atlas workstream.

Everything in here is deterministic and side-effect free except for reads.

Ontology direction convention (important):
    ancestors(t)  = t plus all CL terms above t
    covering_atlas_terms(candidate) = { atlas term t : candidate in ancestors(t) }
i.e. atlas terms that are the candidate itself or a MORE SPECIFIC term under it.
The previous probe script had this inverted, which inflated apparent coverage;
the bug is corrected here and is logged in FINDINGS.md.
"""

from __future__ import annotations

import csv
import math
import os
import re
from collections import defaultdict

import h5py
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
RAW = os.path.join(OUT, "raw")
SC = os.path.join(ROOT, "data", "scRNAseq_raw")
WIDE = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925", "hpa_ifn_landscape_wide.tsv")

MODULE = {
    "IFN_I_capacity": ["IFNAR1", "IFNAR2"],
    "IFN_III_capacity": ["IFNLR1", "IL10RB"],
    "ISG_priming": ["ISG15", "MX1"],
}
SIX = [g for v in MODULE.values() for g in v]
QC_FLOOR = 50  # cells per (atlas, CL term) to enter the pseudobulk

# HPA's own shorthand for a handful of cell types (frozen, from the HPA table
# column values; used only to build the ontology query string).
HPA_ALIASES = {
    "b-cells": "B cell",
    "t-cells": "T cell",
    "nk-cells": "natural killer cell",
    "pdcs": "plasmacytoid dendritic cell, human",
    "cdc": "conventional dendritic cell",
    "microglia": "microglial cell",
    "erythrocytes": "erythrocyte",
    "platelets": "platelet",
    "hepatocytes": "hepatocyte",
    "colonocytes": "colonocyte",
    "podocytes": "podocyte",
    "melanocytes": "melanocyte",
    "cardiomyocytes": "cardiac muscle cell",
    "somatotrophs": "somatotroph",
    "lactotrophs": "mammotroph",
    "corticotrophs": "corticotroph",
    "thyrotrophs": "thyrotroph",
    "gonadotrophs": "gonadtroph",
    "mucous neck cells": "mucous neck cell of gastric gland",
}

# atlas file -> (cell type label column, library size column)
ATLAS_COLS = {
    "gut_APC_lamina_propria.h5ad": ("cell_type", None),
    "HBCA_microglia.h5ad": ("cell_type", "total_UMIs"),
    "HBCA_vascular.h5ad": ("cell_type", "total_UMIs"),
    "HLiCA_cholangiocyte.h5ad": ("cell_type", "nCount_RNA"),
    "HLiCA_endothelial.h5ad": ("cell_type", "nCount_RNA"),
    "HLiCA_hepatocyte.h5ad": ("cell_type", "nCount_RNA"),
    "HLiCA_myeloid.h5ad": ("cell_type", "nCount_RNA"),
    "Krasnow_lung_10X.h5ad": ("cell_type", "nUMI"),
    "synovium_RA.h5ad": ("cell_type", "total_counts"),
    "TS_Vasculature_cellxgene_source.h5ad": ("cell_type", "total_counts"),
    "TS_bone_marrow.h5ad": ("cell_type", "total_counts"),
    "TS_kidney.h5ad": ("cell_type", "total_counts"),
    "TS_large_intestine.h5ad": ("cell_type", "total_counts"),
    "TS_muscle.h5ad": ("cell_type", "total_counts"),
    "TS_prostate.h5ad": ("cell_type", "total_counts"),
    "TS_skin.h5ad": ("cell_type", "total_counts"),
    "TS_small_intestine.h5ad": ("cell_type", "total_counts"),
    "TS_testis.h5ad": ("cell_type", "total_counts"),
    "TS_uterus.h5ad": ("cell_type", "total_counts"),
}


def parse_obo(path: str) -> dict[str, dict[str, object]]:
    terms: dict[str, dict[str, object]] = {}
    cur: dict[str, object] | None = None
    in_term = False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                if cur and cur.get("id"):
                    terms[str(cur["id"])] = cur
                cur = {"name": "", "is_a": [], "synonyms": [], "obsolete": False, "parents_alt": []}
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                if cur and cur.get("id"):
                    terms[str(cur["id"])] = cur
                cur = None
                in_term = False
                continue
            if not in_term or cur is None or not line:
                continue
            if line.startswith("id: "):
                cur["id"] = line[4:].strip()
            elif line.startswith("name: "):
                cur["name"] = line[6:].strip()
            elif line.startswith("is_a: "):
                cur["is_a"].append(line[6:].split("!")[0].strip())
            elif line.startswith("relationship: part_of "):
                cur["parents_alt"].append(line.split()[2].strip())
            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "(.*?)" (\w+)', line)
                if m:
                    cur["synonyms"].append(m.group(1))
            elif line.startswith("is_obsolete: true"):
                cur["obsolete"] = True
    if cur and cur.get("id"):
        terms[str(cur["id"])] = cur
    return terms


class Ontology:
    def __init__(self, obo_path: str):
        self.terms = parse_obo(obo_path)
        parent: dict[str, set[str]] = defaultdict(set)
        for tid, rec in self.terms.items():
            for p in list(rec["is_a"]) + list(rec["parents_alt"]):  # type: ignore[arg-type]
                parent[tid].add(p)
        self.parent = parent
        self._anc: dict[str, set[str]] = {}

    def ancestors(self, tid: str) -> set[str]:
        """tid plus every CL term above it (reflexive transitive closure)."""
        if tid in self._anc:
            return self._anc[tid]
        seen = {tid}
        stack = [tid]
        while stack:
            cur = stack.pop()
            for p in self.parent.get(cur, ()):  # type: ignore[arg-type]
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        self._anc[tid] = seen
        return seen

    def name(self, tid: str) -> str:
        rec = self.terms.get(tid)
        return str(rec["name"]) if rec else ""


def norm(text: str) -> str:
    text = text.lower().strip()
    text = text.replace("alpha-beta", "alphabeta").replace("-", " ")
    text = re.sub(r"[^a-z0-9 ]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    words = []
    for w in text.split(" "):
        if w.endswith("s") and not w.endswith("ss") and len(w) > 3:
            w = w[:-1]
        words.append(w)
    return " ".join(sorted(words))


def decode_h5_node(node) -> np.ndarray:
    """Decode an h5ad obs/var dataset, or categorical group, into an array of str."""
    if isinstance(node, h5py.Group):
        cats = node["categories"][()]
        codes = node["codes"][()]
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in cats]
        return np.array([cats[c] if 0 <= c < len(cats) else "" for c in codes], dtype=object)
    data = node[()]
    if isinstance(data, bytes):
        data = [data]
    return np.array([d.decode() if isinstance(d, bytes) else str(d) for d in data], dtype=object)


def load_hpa() -> dict[str, dict[str, float]]:
    rows = {}
    for r in csv.DictReader(open(WIDE, encoding="utf-8"), delimiter="\t"):
        rows[r["cell_type"]] = {k: float(v) if v not in ("", None) else 0.0 for k, v in r.items() if k != "cell_type"}
    return rows


def hpa_delta(row: dict[str, float]) -> float:
    comps = []
    for genes in MODULE.values():
        comps.append(sum(math.log10(row.get(g, 0.0) + 1.0) for g in genes) / len(genes))
    return sum(comps) / len(comps)


def atlas_inventory() -> list[dict[str, str]]:
    return list(
        csv.DictReader(open(os.path.join(RAW, "atlas_celltype_inventory.tsv"), encoding="utf-8"), delimiter="\t")
    )


def spearman(a: list[float], b: list[float]) -> float:
    def rank(x: list[float]) -> list[float]:
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    da = math.sqrt(sum((v - ma) ** 2 for v in ra))
    db = math.sqrt(sum((v - mb) ** 2 for v in rb))
    return num / (da * db) if da and db else float("nan")


def spearman_p_value(rho: float, n: int) -> float:
    """Two-sided t-approximation (documented; permutation p is the headline)."""
    if n < 3 or not math.isfinite(rho) or abs(rho) == 1.0:
        return 0.0 if abs(rho) == 1.0 and n >= 3 else float("nan")
    t = rho * math.sqrt((n - 2) / (1 - rho * rho))
    # two-sided p from Student t via the incomplete beta function
    from scipy import stats

    return float(2 * stats.t.sf(abs(t), df=n - 2))


def zscore(vals: dict[str, float]) -> dict[str, float]:
    xs = list(vals.values())
    mu = sum(xs) / len(xs)
    sd = math.sqrt(sum((v - mu) ** 2 for v in xs) / len(xs)) or 1.0
    return {k: (v - mu) / sd for k, v in vals.items()}


def min_detectable_rho(n: int, power: float = 0.8, alpha: float = 0.05) -> float:
    """Smallest |Spearman rho| detectable with the given power (Fisher z approximation)."""
    if n < 4:
        return float("nan")
    from scipy import stats

    z_alpha = stats.norm.isf(alpha / 2)
    z_power = stats.norm.isf(1 - power)
    z_r = (z_alpha + z_power) / math.sqrt(n - 3)
    return float(math.tanh(z_r))


def achieved_power(rho: float, n: int, alpha: float = 0.05) -> float:
    if n < 4 or not math.isfinite(rho) or abs(rho) >= 1:
        return float("nan")
    from scipy import stats

    z_alpha = stats.norm.isf(alpha / 2)
    z_r = math.atanh(abs(rho)) * math.sqrt(n - 3)
    return float(stats.norm.sf(z_alpha - z_r) + stats.norm.cdf(-z_alpha - z_r))
