"""Step 1: inspect the GSE342661 TPM matrix (header, columns, row symbols).

Frozen criterion (written before seeing the values):
  * the matrix is usable only if every column maps back to a treatment condition
    parsed from the sample title, and the row identifier column is a gene symbol
    column we can spot-check by name.

Read-only. Writes only g342_matrix_shape.json in this directory.
"""

import gzip
import json
import pathlib
import re

RAW = pathlib.Path(__file__).resolve().parent / "raw" / "GSE342661_counts_tpm.matrix.gz"
OUT = pathlib.Path(__file__).resolve().parent / "g342_matrix_shape.json"


def main() -> None:
    with gzip.open(RAW, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        first = fh.readline()
        n_rows = 1
        probes = []
        for line in fh:
            n_rows += 1
            if len(probes) < 5:
                probes.append(line.split("\t", 1)[1] if False else line.split("\t")[0])

    cols = header[1:]
    n_rep_time = re.compile(r"^(?P<time>\d+h)\s*,\s*(?P<trt>[^,]+?)\s*,\s*biol rep\s*(?P<rep>\d+)\s*,\s*(?P<route>[ab]+)\s*$")
    parsed = {}
    unmapped = []
    for c in cols:
        m = n_rep_time.match(c.strip())
        if m:
            key = (m.group("time"), m.group("trt").strip(), m.group("route"))
            parsed[key] = parsed.get(key, 0) + 1
        else:
            unmapped.append(c)

    shape = {
        "file": str(RAW),
        "n_columns_including_id": len(header),
        "first_header_cell": header[0],
        "n_sample_columns": len(cols),
        "n_data_rows": n_rows,
        "first_row_id": first.split("\t")[0],
        "first_row_id_column_name": header[0],
        "first_sample_labels": cols[:6],
        "last_sample_labels": cols[-6:],
        "n_titles_matched": sum(parsed.values()),
        "n_titles_unmatched": len(unmapped),
        "unmatched_titles": unmapped[:20],
        "design_grid": {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in sorted(parsed.items())},
        "distinct_times": sorted({k[0] for k in parsed}),
        "distinct_treatments": sorted({k[1] for k in parsed}),
        "distinct_routes": sorted({k[2] for k in parsed}),
    }
    OUT.write_text(json.dumps(shape, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({k: v for k, v in shape.items() if k != "design_grid"}, indent=2, ensure_ascii=False))
    print("grid cells:", len(shape["design_grid"]))
    for k, v in shape["design_grid"].items():
        print(" ", k, v)


if __name__ == "__main__":
    main()
