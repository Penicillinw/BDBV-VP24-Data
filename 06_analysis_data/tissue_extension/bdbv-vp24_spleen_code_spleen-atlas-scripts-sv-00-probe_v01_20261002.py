"""SRB-VERIFY probe: inspect the inputs needed for an independent recompute.

Read-only.  Purpose: fix the exact column names, h5ad layout and label-map
semantics BEFORE writing the verification scripts (judgement criteria first).

Writes: analysis/srb_verify_20260926/raw/probe_structures.json
"""

import json
import os

import h5py

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUT = os.path.join(ROOT, "analysis", "srb_verify_20260926")
RAW = os.path.join(OUT, "raw")
os.makedirs(RAW, exist_ok=True)

T4 = os.path.join(ROOT, "analysis", "t4_ifn_landscape_20260925")
B5 = os.path.join(ROOT, "analysis", "b5_harmonised_atlas_20260926")
B6 = os.path.join(ROOT, "analysis", "b6_coverage_20260926")
SRB = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926")
H2 = os.path.join(ROOT, "analysis", "h2_lowend_atlas_hunt_20260926")
HCL = os.path.join(H2, "raw", "atlas_h5ad",
                   "2adb1f8a_Construction_of_a_human_cell_landscape_at_single.h5ad")


def head(path, n=1):
    with open(path, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = [fh.readline().rstrip("\n").split("\t") for _ in range(n)]
    return hdr, rows


info = {}

for name in ["hpa_ifn_landscape_wide.tsv", "restriction_gradient.tsv"]:
    p = os.path.join(T4, name)
    hdr, rows = head(p)
    info[name] = {"n_cols": len(hdr), "header": hdr, "first_row": rows[0]}

p = os.path.join(T4, "restriction_gradient_summary.json")
with open(p, encoding="utf-8") as fh:
    summ = json.load(fh)
info["restriction_gradient_summary.keys"] = list(summ.keys())
info["restriction_gradient_summary.gene_stats"] = summ.get("gene_stats")

for name in ["label_mapping.tsv", "compartment_coverage.tsv"]:
    p = os.path.join(B5, name)
    if os.path.exists(p):
        hdr, rows = head(p)
        info["b5/" + name] = {"n_cols": len(hdr), "header": hdr, "first_row": rows[0]}
    else:
        info["b5/" + name] = "MISSING"

if os.path.isdir(B6):
    info["b6_files"] = sorted(os.listdir(B6))
else:
    info["b6_files"] = "MISSING"

for name in ["hcl_spleen_pseudobulk.tsv", "srb_anchor_check.tsv",
             "srb_spleen_members.tsv", "srb_null.json",
             "srb_independent_test.json", "ts_spleen_pseudobulk.tsv",
             "srb_platform_offset.tsv", "srb_placement.json"]:
    p = os.path.join(SRB, "raw", name)
    if p.endswith(".json"):
        with open(p, encoding="utf-8") as fh:
            info["srb/" + name] = json.load(fh)
    else:
        hdr, rows = head(p, 3)
        info["srb/" + name] = {"header": hdr, "rows": rows}

info["hcl_path"] = HCL
info["hcl_exists"] = os.path.exists(HCL)
info["hcl_bytes"] = os.path.getsize(HCL) if os.path.exists(HCL) else None

with h5py.File(HCL, "r") as f:
    keys = []
    f.visit(lambda k: keys.append(k))
    info["hcl_top"] = sorted({k.split("/")[0] for k in keys})
    info["hcl_raw_children"] = sorted(f["raw"].keys()) if "raw" in f else None
    info["hcl_X_type"] = type(f["raw/X"]).__name__ if "raw" in f else type(f["X"]).__name__
    node = f["raw/X"] if ("raw" in f and "X" in f["raw"]) else f["X"]
    if hasattr(node, "shape"):
        info["hcl_X_shape"] = list(node.shape)
        info["hcl_X_kind"] = "dense_dataset"
    else:
        info["hcl_X_children"] = sorted(node.keys())
        info["hcl_X_shape"] = list(node.attrs.get("shape", []))
        info["hcl_X_kind"] = "sparse_group"
        for ch in ("data", "indices", "indptr"):
            if ch in node:
                info["hcl_X_" + ch + "_dtype"] = str(node[ch].dtype)
                info["hcl_X_" + ch + "_len"] = int(node[ch].shape[0])
        if "data" in node and node["data"].shape[0] and len(node["data"].shape) > 1:
            info["hcl_X_encoding_type"] = str(node.attrs.get("encoding-type", ""))
    vg = f["raw/var"] if "raw" in f else f["var"]
    info["hcl_var_children"] = sorted(vg.keys()) if hasattr(vg, "keys") else None
    info["hcl_var_attrs"] = {k: str(v) for k, v in vg.attrs.items()}
    ob = f["obs"]
    info["hcl_obs_children"] = sorted(ob.keys())[:60]
    info["hcl_obs_n"] = int(ob.attrs.get("_index", "_index") is not None and len(ob[ob.attrs["_index"]]))
    info["hcl_obs_index_name"] = ob.attrs.get("_index")
    for col in ["tissue", "cell_type_ontology_term_id", "author_cell_type",
                "cell_type", "donor_id"]:
        if col in ob:
            node2 = ob[col]
            if isinstance(node2, h5py.Group):
                codes = node2["codes"][:]
                cats = [(c.decode() if isinstance(c, bytes) else str(c))
                        for c in node2["categories"][:]]
                uniq = sorted(set(cats))
                n = len(codes)
                kind = "categorical"
            else:
                arr = node2[:]
                vals = [(a.decode() if isinstance(a, bytes) else str(a)) for a in arr]
                uniq = sorted(set(vals))
                n = len(arr)
                kind = "array"
            info["hcl_obs_" + col] = {
                "kind": kind,
                "n": n,
                "n_unique": len(uniq),
                "first10": uniq[:10],
                "spleen_terms": [u for u in uniq if "spleen" in u.lower()][:20],
            }

with open(os.path.join(RAW, "probe_structures.json"), "w", encoding="utf-8") as fh:
    json.dump(info, fh, ensure_ascii=False, indent=1, default=str)
print("wrote raw/probe_structures.json")
print(json.dumps({k: v for k, v in info.items()
                  if k.startswith("hcl_") and k != "hcl_path"}, ensure_ascii=False, indent=1)[:4000])
