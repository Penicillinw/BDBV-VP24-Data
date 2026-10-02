"""LRS-01  Reference-side coverage scan (focus: the low-end compartments).

Question this answers
---------------------
The cross-atlas reproduction of the restriction gradient connects a reference
cell type to an atlas label through Cell Ontology identifiers.  A reference
class can therefore fail to appear in the reproduction for two very different
reasons: because the frozen reference-side mapping never gave it an ontology
term (a *reference-side* gap, in principle repairable), or because no atlas
label denotes it (an *atlas-side* gap, not repairable without new data).
H2C found this for the retina by hand: `rod photoreceptor cells` and
`cone photoreceptor cells` had no frozen mapping at all, and supplying one
connected both to real atlas labels.  This script generalises that diagnosis
to every reference class, and splits the failures further into a third,
strictly mechanical category: pairs whose within-atlas z is identically zero
because the atlas has no contrast at all inside the reference class.

Frozen criteria (fixed before any value in this analysis was read; any change
requires a separate revision record, per the project freeze file section 5)
-----------------------------------------------------------------------------
R1  Ontology.  Cell Ontology read from
    analysis/b5_harmonised_atlas_20260926/raw/cl.obo.  For every [Term] block
    keep id, name, synonyms, is_obsolete, and parents = is_a UNION
    relationship: part_of.  Qualifier suffixes are stripped by taking the
    identifier token only, which is the repair of the b5_common.parse_obo
    defect logged in the freeze file section 4 (1,320 is_a edges were being
    dropped silently).  ancestors(t) is the reflexive transitive closure.
R2  Covering rule (identical to B5 / H2B).  An atlas label with ontology term
    a covers the reference class with term c iff c is in ancestors(a), i.e.
    the atlas label is c itself or a more specific term below it.
R3  Floor.  An atlas label enters a pair only if n_cells >= 50.
R4  Connection.  A reference class is CONNECTED in a scope iff at least one
    floor-passing atlas label of that scope covers it.
R5  Mapping route.  A reference class present in the frozen B5 mapping keeps
    its recorded route (exact_label / synonym / curated).
R6  Witness search (mechanical; no hand-picked ontology ids except the frozen
    alias table below).  For a reference class absent from the frozen mapping,
    build the witness set = the ontology term of every floor-passing atlas
    label in the scope PLUS every ancestor of those terms.  A witness term is
    a match if it satisfies one of these frozen tiers, evaluated in order:
        T1 exact        norm(name) == norm(reference label)
        T2 synonym      norm(synonym) == norm(reference label)
        T2.5 alias      norm(name) == norm(frozen ALIASES pointer for the label)
        T3 containment  every token of the shorter normalised name is a token
                        of the longer one, and the shorter name carries at
                        least one token outside the frozen STOP_TOKENS set.
                        The stop-token condition is what keeps the ontology
                        root ("cell") from matching every reference label that
                        happens to contain the word cell; without it T3 is a
                        false-positive generator and the band is worthless.
    norm() = lower-case, NFKD ASCII-fold, non-alphanumerics to space, collapse,
    singularise tokens of length > 3 ending in "s" but not "ss", sort tokens.
    Obsolete terms are never witnesses.  The witness with the lowest tier,
    then the largest n_cells, then the smallest term id, is reported.
R7  Classification of a (reference class, scope) cell.  T3 is reported
    SEPARATELY from T1/T2/T2.5 because a containment match is a name-similarity
    proxy for the ancestry rule R2, not a guarantee of it; only T1/T2/T2.5
    assert that the atlas label really sits below a term naming the reference
    class.  The split is deliberate: it bounds recall honestly, because H2C's
    hand-curated rod/cone repair (retinal rod cell / retinal cone cell) is
    NOT reachable by any mechanical tier and shows that the T3 band is not an
    upper bound on what curation can add.
        MAPPED_CONNECTED              frozen mapping, >= 1 floor-passing covering label
        MAPPED_UNCONNECTED            frozen mapping, no floor-passing covering label
        UNMAPPED_NAME_IDENTICAL       no frozen mapping, T1/T2/T2.5 witness
        UNMAPPED_CONTAINMENT_CANDIDATE no frozen mapping, T3 witness only
        UNMAPPED_NO_WITNESS           no frozen mapping, no witness in this scope
R8  Mechanically degenerate pair.  A (atlas, reference class) pair is flagged
    MECHANICAL_ZERO iff its reported z is forced to 0 by construction, i.e.
    the covering labels are the entire floor-passing label set of that atlas,
    or the atlas has exactly one floor-passing label.  Such a pair carries no
    information about the reference value, whatever that value is.
R9  Scopes.  Every reference class is reported under BOTH scopes:
        S19 = the 19 harmonised atlases of B5,
        S26 = S19 plus the 7 low-end atlases downloaded in H2.
R10 LOW_END set (frozen).  The union of
        (a) the bottom quartile of the frozen Delta ranking, rank >= 116, and
        (b) the six compartments the H2 hunt went looking for:
            urothelial cells, retinal pigment epithelial cells,
            rod photoreceptor cells, cone photoreceptor cells,
            Sertoli cells, choroid plexus epithelial cells.

Nothing in the frozen Delta ranking, the frozen mapping file or any atlas
pseudobulk is modified; this script only reads them and writes into
analysis/lowend_reference_scan_20260926/.
"""

from __future__ import annotations

import csv
import json
import os
import re
import unicodedata
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
RAW_B5 = os.path.join(B5, "raw")
H2 = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926")
T4 = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")

OBO = os.path.join(RAW_B5, "cl.obo")
REF = os.path.join(T4, "restriction_gradient.tsv")
MAP = os.path.join(B5, "label_mapping.tsv")
INV19 = os.path.join(RAW_B5, "atlas_celltype_inventory.tsv")
LAB7 = os.path.join(H2, "raw", "new_atlas_labels.tsv")

FLOOR = 50
BOTTOM_QUARTILE_RANK = 116  # rank >= 116 -> bottom 39 of 154
LOWEND_HUNT = [
    "urothelial cells",
    "retinal pigment epithelial cells",
    "rod photoreceptor cells",
    "cone photoreceptor cells",
    "sertoli cells",
    "choroid plexus epithelial cells",
]

# R6 / T2.5: frozen alias table.  Only transliteration- and spelling-level
# aliases are allowed here; every entry must state its evidence.  Anything
# that would need a biological judgement call belongs in a separate curate step
# and must not be smuggled in through this table.
ALIASES = {
    # HPA spells the retinal glia with a diaeresis, CL with the ASCII
    # transliteration; after NFKD folding the two still differ (uller/mueller).
    "müller glia": ("Mueller cell", "ASCII transliteration of the same eponym"),
    "muller glia": ("Mueller cell", "ASCII transliteration of the same eponym"),
}

_ALIAS_BY_NORM: dict[str, str] = {}

# Tokens that on their own carry no cell-type information.  A containment match
# whose shorter side consists only of these is discarded (R6 / T3).
STOP_TOKENS = {"cell", "human", "type"}


def alias_target_norm(ref_label: str) -> str:
    """Lazy lookup so that this can be defined above norm()."""
    if not _ALIAS_BY_NORM:
        _ALIAS_BY_NORM.update({norm(k): v[0] for k, v in ALIASES.items()})
    return _ALIAS_BY_NORM.get(norm(ref_label), "")


# --------------------------------------------------------------------------
# ontology
# --------------------------------------------------------------------------
def parse_obo(path: str) -> dict[str, dict[str, object]]:
    """R1: qualifier-safe parser (repair of the b5_common defect)."""
    terms: dict[str, dict[str, object]] = {}
    cur: dict[str, object] | None = None
    in_term = False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                if cur and cur.get("id"):
                    terms[str(cur["id"])] = cur
                cur = {"name": "", "synonyms": [], "parents": [], "obsolete": False}
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
                cur["parents"].append(line.split()[1])  # type: ignore[union-attr]
            elif line.startswith("relationship: part_of "):
                cur["parents"].append(line.split()[2])  # type: ignore[union-attr]
            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "(.*?)" (\w+)', line)
                if m:
                    cur["synonyms"].append(m.group(1))  # type: ignore[union-attr]
            elif line.startswith("is_obsolete: true"):
                cur["obsolete"] = True
    if cur and cur.get("id"):
        terms[str(cur["id"])] = cur
    return terms


def build_ancestors(terms: dict[str, dict[str, object]]) -> dict[str, set[str]]:
    parent = {t: set(v["parents"]) for t, v in terms.items()}  # type: ignore[arg-type]
    memo: dict[str, set[str]] = {}

    def closure(t: str) -> set[str]:
        if t in memo:
            return memo[t]
        seen, stack = {t}, [t]
        while stack:
            cur = stack.pop()
            for p in parent.get(cur, ()):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        memo[t] = seen
        return seen

    for t in list(parent):
        closure(t)
    return memo


# --------------------------------------------------------------------------
# normalisation (R6)
# --------------------------------------------------------------------------
def norm(text: str) -> str:
    text = unicodedata.normalize("NFKD", str(text))
    text = "".join(c for c in text if not unicodedata.combining(c))
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


def tokens(text: str) -> set[str]:
    return set(norm(text).split(" ")) - {""}


def tier_match(ref_label: str, term_name: str, term_syns: list[str]) -> int | None:
    """Return the frozen tier (1, 2, 2.5 or 3) or None."""
    if norm(term_name) == norm(ref_label):
        return 1
    if any(norm(s) == norm(ref_label) for s in term_syns):
        return 2
    alias_target = alias_target_norm(ref_label)
    if alias_target and norm(term_name) == norm(alias_target):
        return 2.5
    a, b = tokens(ref_label), tokens(term_name)
    if not a or not b:
        return None
    short, long_ = (a, b) if len(a) <= len(b) else (b, a)
    if short <= long_ and (short - STOP_TOKENS):
        return 3
    return None


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------
def load_reference() -> dict[str, dict[str, float | int]]:
    out = {}
    for r in csv.DictReader(open(REF, encoding="utf-8"), delimiter="\t"):
        out[r["cell_type"]] = {
            "delta": float(r["delta_restriction"]),
            "rank": int(r["rank"]),
        }
    return out


def load_mapping() -> dict[str, dict[str, str]]:
    out = {}
    for r in csv.DictReader(open(MAP, encoding="utf-8"), delimiter="\t"):
        out[r["hpa_cell_type"]] = {
            "cl_id": r["cl_id"],
            "cl_label": r["cl_label"],
            "route": r["mapping_route"],
            "testable": r["testable"],
        }
    return out


def load_labels() -> list[dict[str, object]]:
    labels: list[dict[str, object]] = []
    for r in csv.DictReader(open(INV19, encoding="utf-8"), delimiter="\t"):
        labels.append({
            "origin": "S19",
            "atlas": r["atlas"],
            "cl_id": r["cl_id"],
            "cl_label": r["cl_label"],
            "n_cells": int(float(r["n_cells"])),
        })
    for r in csv.DictReader(open(LAB7, encoding="utf-8"), delimiter="\t"):
        if not r["cl_id"].startswith("CL:"):
            continue  # labels the atlas could not be mapped to an ontology term
        labels.append({
            "origin": "new7",
            "atlas": r["atlas"],
            "cl_id": r["cl_id"],
            "cl_label": r["cl_label"],
            "n_cells": int(float(r["n_cells"])),
        })
    return labels


def in_scope(label: dict[str, object], scope: str) -> bool:
    """S19 = the 19 harmonised atlases; S26 = S19 plus the 7 H2 low-end atlases."""
    return True if scope == "S26" else label["origin"] == "S19"


# --------------------------------------------------------------------------
# main scan
# --------------------------------------------------------------------------
def main() -> None:
    terms = parse_obo(OBO)
    anc = build_ancestors(terms)
    ref = load_reference()
    mapping = load_mapping()
    labels = load_labels()

    print(f"ontology terms {len(terms)}; reference classes {len(ref)}; "
          f"frozen mapping rows {len(mapping)}; atlas labels {len(labels)}")

    # R2/R3: for each scope, which reference classes each atlas label covers.
    # reference class -> cl term (frozen mapping) -> covering atlas labels
    rows: list[dict[str, object]] = []
    for hpa, info in sorted(ref.items()):
        for scope in ("S19", "S26"):
            scope_labels = [l for l in labels if in_scope(l, scope)]
            mapped = hpa in mapping
            cl_id = mapping[hpa]["cl_id"] if mapped else ""
            route = mapping[hpa]["route"] if mapped else ""
            covering = []
            if mapped:
                for l in scope_labels:
                    if int(l["n_cells"]) < FLOOR:
                        continue
                    if cl_id in anc.get(str(l["cl_id"]), {str(l["cl_id"])}):
                        covering.append(l)
            if mapped:
                cls = "MAPPED_CONNECTED" if covering else "MAPPED_UNCONNECTED"
                witness = ""
                wit_tier = ""
                wit_cells = ""
                wit_labels: list[str] = []
            else:
                # R6: search the witness set = terms of floor-passing labels + ancestors
                best: tuple[tuple[float, int, str], float, str, str, dict] | None = None
                for l in scope_labels:
                    if int(l["n_cells"]) < FLOOR:
                        continue
                    chain = anc.get(str(l["cl_id"]), {str(l["cl_id"])})
                    for t in chain:
                        rec = terms.get(t)
                        if not rec or rec.get("obsolete"):
                            continue
                        tier = tier_match(hpa, str(rec["name"]), list(rec["synonyms"]))  # type: ignore[arg-type]
                        if tier is None:
                            continue
                        # Prefer the most specific witness term (deepest in the
                        # ontology) before the largest label: a coarse ancestor
                        # such as "retinal cell" is technically a witness for
                        # several reference classes and carries less information
                        # than the label's own term would.
                        depth = len(anc.get(t, {t}))
                        key = (tier, -depth, -int(l["n_cells"]), t)
                        if best is None or key < best[0]:
                            best = (key, tier, t, str(rec["name"]), l)
                if best is None:
                    cls = "UNMAPPED_NO_WITNESS"
                    witness = ""
                    wit_tier = ""
                    wit_cells = ""
                else:
                    _, tier, t, tname, l = best
                    cls = ("UNMAPPED_NAME_IDENTICAL" if tier <= 2.5
                           else "UNMAPPED_CONTAINMENT_CANDIDATE")
                    witness = f"{t} {tname} @ {l['atlas']}:{l['cl_label']}"
                    wit_tier = tier
                    wit_cells = int(l["n_cells"])
                    # A tier-3 witness term is often coarser than the reference
                    # label, so several atlas labels can sit underneath it (e.g.
                    # both rod and cone sit under "photoreceptor cell").  Report
                    # every matched label, not just the single best one.
                    wit_labels = []
                    for l2 in scope_labels:
                        if int(l2["n_cells"]) < FLOOR:
                            continue
                        if t in anc.get(str(l2["cl_id"]), {str(l2["cl_id"])}):
                            wit_labels.append(
                                f"{l2['atlas'][:18]}:{l2['cl_label']}({l2['n_cells']})")
            # Independent of the ontology route: is there an atlas label that is
            # merely NAMED like the reference class?  A hit that is not also a
            # covering label means the cells are present but the label was never
            # annotated with a term that connects, i.e. an annotation gap rather
            # than an absence of the cell type.
            covering_ids = {(str(l["atlas"]), str(l["cl_id"])) for l in covering}
            name_hits = []
            coarser_hits = []
            for l in scope_labels:
                if int(l["n_cells"]) < FLOOR:
                    continue
                if (str(l["atlas"]), str(l["cl_id"])) in covering_ids:
                    continue
                if tier_match(hpa, str(l["cl_label"]), []) is not None:
                    name_hits.append(f"{l['atlas'][:18]}:{l['cl_label']}({l['n_cells']})")
                # Granularity mismatch: the atlas label is COARSER than the
                # reference class (the atlas term is a strict ancestor of the
                # reference term), so the frozen descendant direction cannot
                # connect it even though the cell type is plainly present.
                if mapped and str(l["cl_id"]) in (anc.get(cl_id, {cl_id}) - {cl_id}):
                    coarser_hits.append(f"{l['atlas'][:18]}:{l['cl_label']}({l['n_cells']})")
            rows.append({
                "reference_cell_type": hpa,
                "scope": scope,
                "reference_delta": round(float(ref[hpa]["delta"]), 4),
                "reference_rank": ref[hpa]["rank"],
                "low_end": "yes" if (int(ref[hpa]["rank"]) >= BOTTOM_QUARTILE_RANK
                                     or hpa in LOWEND_HUNT) else "no",
                "frozen_mapping": "yes" if mapped else "no",
                "frozen_route": route,
                "frozen_cl_id": cl_id,
                "frozen_testable": mapping[hpa]["testable"] if mapped else "",
                "classification": cls,
                "n_covering_labels": len(covering),
                "n_covering_atlases": len({str(l["atlas"]) for l in covering}),
                "covering_atlases": ";".join(sorted({str(l["atlas"]) for l in covering})),
                "covering_cells": sum(int(l["n_cells"]) for l in covering),
                "covering_terms": ";".join(
                    f"{l['atlas'][:18]}:{l['cl_id']}({l['n_cells']})" for l in covering[:8]),
                "witness_tier": wit_tier,
                "witness": witness,
                "witness_cells": wit_cells,
                "witness_labels": ";".join(wit_labels[:8]),
                "n_name_level_hits": len(name_hits),
                "name_level_hits": ";".join(name_hits[:6]),
                "n_coarser_label_hits": len(coarser_hits),
                "coarser_label_hits": ";".join(coarser_hits[:6]),
            })

    with open(os.path.join(HERE, "reference_coverage_scan.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- summary
    summary: dict[str, object] = {"floor_cells": FLOOR,
                                  "bottom_quartile_rank": BOTTOM_QUARTILE_RANK,
                                  "n_reference_classes": len(ref),
                                  "n_frozen_mapped": len(mapping),
                                  "n_atlas_labels": len(labels)}
    for scope in ("S19", "S26"):
        sub = [r for r in rows if r["scope"] == scope]
        counts: dict[str, int] = defaultdict(int)
        for r in sub:
            counts[str(r["classification"])] += 1
        low = [r for r in sub if r["low_end"] == "yes"]
        low_counts: dict[str, int] = defaultdict(int)
        for r in low:
            low_counts[str(r["classification"])] += 1
        summary[scope] = {
            "classification": dict(sorted(counts.items())),
            "n_low_end": len(low),
            "low_end_classification": dict(sorted(low_counts.items())),
        }

    # ---- degeneracy diagnostic (R8)
    pb_files = [(os.path.join(RAW_B5, "atlas_pseudobulk.tsv"), "S19"),
                (os.path.join(H2, "raw", "new_atlas_pseudobulk.tsv"), "new7")]
    diag_rows: list[dict[str, object]] = []
    for path, origin in pb_files:
        recs = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
        by_atlas: dict[str, list[dict[str, str]]] = defaultdict(list)
        for r in recs:
            by_atlas[r["atlas"]].append(r)
        for atlas, rs in sorted(by_atlas.items()):
            floor_rs = [r for r in rs if int(float(r["n_cells"])) >= FLOOR]
            diag_rows.append({
                "origin": origin, "atlas": atlas, "n_labels_total": len(rs),
                "n_labels_ge_floor": len(floor_rs),
                "atlas_has_contrast": "no" if len(floor_rs) <= 1 else "yes",
            })
    with open(os.path.join(HERE, "atlas_contrast_diagnostic.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(diag_rows[0].keys()))
        w.writeheader()
        w.writerows(diag_rows)
    summary["atlases_without_contrast"] = [
        r["atlas"] for r in diag_rows if r["atlas_has_contrast"] == "no"]

    # ---- what the 7 low-end atlases actually closed (low-end subset only)
    by_class: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for r in rows:
        by_class[str(r["reference_cell_type"])][str(r["scope"])] = r
    flipped, still_open, low_rows = [], [], []
    for hpa, sc in sorted(by_class.items()):
        if sc["S19"]["low_end"] != "yes":
            continue
        low_rows.append(hpa)
        a, b = str(sc["S19"]["classification"]), str(sc["S26"]["classification"])
        if a == "MAPPED_UNCONNECTED" and b == "MAPPED_CONNECTED":
            flipped.append({"reference_cell_type": hpa,
                            "reference_rank": sc["S26"]["reference_rank"],
                            "reference_delta": sc["S26"]["reference_delta"],
                            "atlases": sc["S26"]["covering_atlases"]})
        if b != "MAPPED_CONNECTED":
            still_open.append({"reference_cell_type": hpa,
                               "reference_rank": sc["S26"]["reference_rank"],
                               "classification_S26": b,
                               "frozen_route": sc["S26"]["frozen_route"] or "(no mapping)",
                               "name_level_hits": sc["S26"]["n_name_level_hits"]})
    summary["low_end_closure"] = {
        "n_low_end": len(low_rows),
        "n_connected_in_S19": sum(1 for h in low_rows
                                  if by_class[h]["S19"]["classification"] == "MAPPED_CONNECTED"),
        "n_connected_in_S26": sum(1 for h in low_rows
                                  if by_class[h]["S26"]["classification"] == "MAPPED_CONNECTED"),
        "n_closed_by_the_7_new_atlases": len(flipped),
        "closed": flipped,
        "n_still_not_connected_in_S26": len(still_open),
        "still_not_connected": still_open,
    }

    # ---- granularity diagnostic: classes whose atlas label exists but is
    # coarser than the reference class (frozen direction cannot reach it)
    gran = [r for r in rows if r["scope"] == "S26" and int(r["n_coarser_label_hits"]) > 0]
    summary["granularity_mismatch_S26"] = {
        "n_classes": len(gran),
        "n_low_end": sum(1 for r in gran if r["low_end"] == "yes"),
        "classes": sorted(str(r["reference_cell_type"]) for r in gran),
    }

    # ---- reference-side gap list
    gaps = [r for r in rows if str(r["classification"]).startswith("UNMAPPED")]
    with open(os.path.join(HERE, "reference_side_gaps.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(gaps)

    with open(os.path.join(HERE, "scan_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    # ---- self-check A: reproduce the frozen pipeline's pair count.  The frozen
    # agreement script (b5_05_agreement.py) restricts to testable == yes and the
    # hierarchy route; the post-OBO-fix replay (b5b_obo_fix_20260926) reports
    # 162 pairs over 18 atlases and 39 distinct reference classes.  This script
    # uses the same rule, so it must land on exactly those numbers; if it does
    # not, nothing downstream of the scan is trustworthy.
    selfcheck: dict[str, object] = {}
    for scope, expect in (("S19", {"pairs": 162, "atlases": 18, "classes": 39}),):
        sub = [r for r in rows if r["scope"] == scope and r["frozen_testable"] == "yes"]
        classes_with_pairs = [r for r in sub if int(r["n_covering_atlases"]) > 0]
        n_pairs = sum(int(r["n_covering_atlases"]) for r in sub)
        atlases = {a for r in sub for a in str(r["covering_atlases"]).split(";") if a}
        selfcheck[scope] = {
            "recomputed_pairs": n_pairs,
            "recomputed_atlases_with_pairs": len(atlases),
            "recomputed_classes_with_pairs": len(classes_with_pairs),
            "frozen_b5b_reference": expect,
            "reproduces_frozen_pipeline": (
                n_pairs == expect["pairs"]
                and len(atlases) == expect["atlases"]
                and len(classes_with_pairs) == expect["classes"]),
        }
    summary["selfcheck_A_vs_frozen_pipeline"] = selfcheck

    # ---- self-check B against the frozen mapping's own coverage columns.  Any
    # difference is the residual effect of the b5_common OBO parse defect (the
    # frozen columns were computed with the buggy parser).
    frozen_cols = {}
    for r in csv.DictReader(open(MAP, encoding="utf-8"), delimiter="\t"):
        frozen_cols[r["hpa_cell_type"]] = {
            "n_ge_floor": int(r["n_covering_ge_floor"] or 0),
            "testable": r["testable"],
        }
    diffs = []
    for r in rows:
        if r["scope"] != "S19" or r["frozen_mapping"] != "yes":
            continue
        mine = int(r["n_covering_labels"])
        theirs = frozen_cols[str(r["reference_cell_type"])]["n_ge_floor"]
        if mine != theirs:
            diffs.append({
                "reference_cell_type": r["reference_cell_type"],
                "frozen_covering_ge_floor": theirs,
                "recomputed_with_fixed_parser": mine,
            })
    summary["selfcheck_vs_frozen_mapping_ge_floor"] = {
        "n_compared": sum(1 for r in rows if r["scope"] == "S19"
                          and r["frozen_mapping"] == "yes"),
        "n_identical": sum(1 for r in rows if r["scope"] == "S19"
                           and r["frozen_mapping"] == "yes") - len(diffs),
        "n_differing": len(diffs),
        "differences": diffs,
    }
    with open(os.path.join(HERE, "scan_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print("self-check vs frozen mapping:",
          summary["selfcheck_vs_frozen_mapping_ge_floor"]["n_identical"], "identical,",
          summary["selfcheck_vs_frozen_mapping_ge_floor"]["n_differing"], "differing")
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
