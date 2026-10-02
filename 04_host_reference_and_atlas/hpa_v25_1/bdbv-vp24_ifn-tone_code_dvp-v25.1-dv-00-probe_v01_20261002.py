# -*- coding: utf-8 -*-
"""dv_00_probe: read-only reconnaissance of the HPA DVP Deep Visual Proteomics matrices.

Owner: subagent `/root/dvp_delta_alignment`.
Writes nothing outside analysis/dvp_delta_alignment_20260926/.
"""
import io
import os
import sys
import zipfile

BASE = r"G:\本迪布焦研究\analysis\protein_layer_scan_20260926\raw\hpa_dvp"
OUT = r"G:\本迪布焦研究\analysis\dvp_delta_alignment_20260926"

FILES = [
    "dvp_cell_type_group_data.tsv.zip",
    "dvp_cell_type.tsv.zip",
    "dvp_sample_data.tsv.zip",
]


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = []
    lines.append("python %s" % sys.version.replace("\n", " "))
    for z in FILES:
        path = os.path.join(BASE, z)
        lines.append("")
        lines.append("=== %s (exists=%s, bytes=%s)" % (z, os.path.exists(path),
                                                     os.path.getsize(path) if os.path.exists(path) else "-"))
        with zipfile.ZipFile(path) as zf:
            lines.append("    members: %s" % zf.namelist())
            for name in zf.namelist():
                if name.endswith("/"):
                    continue
                with zf.open(name) as fh:
                    txt = io.TextIOWrapper(fh, encoding="utf-8")
                    for i, line in enumerate(txt):
                        if i >= 4:
                            break
                        lines.append("    [%s] %s" % (i, line.rstrip()[:400]))
                break
    out = "\n".join(lines)
    print(out)
    with open(os.path.join(OUT, "dv_00_probe.txt"), "w", encoding="utf-8") as fh:
        fh.write(out + "\n")


if __name__ == "__main__":
    main()
