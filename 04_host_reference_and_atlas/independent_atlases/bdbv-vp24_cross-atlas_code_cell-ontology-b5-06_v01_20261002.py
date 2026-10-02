"""B5 step 6 - pointer file for the 19 cached atlases.

Records where each h5ad came from (all 19 are CZ CELLxGENE Discover datasets,
schema 7.1.0) so the harmonised-label analysis is reproducible from the source
collections without shipping the multi-GB matrices inside the deliverable.
"""

from __future__ import annotations

import csv
import glob
import os

import h5py

from b5_common import RAW, SC


def main() -> None:
    rows = []
    for path in sorted(glob.glob(os.path.join(SC, "*.h5ad"))):
        with h5py.File(path, "r") as h:
            uns = h["uns"]
            rec = {"atlas_file": os.path.basename(path), "bytes": os.path.getsize(path)}
            for key in ("schema_version", "title", "citation"):
                if key in uns:
                    v = uns[key][()]
                    rec[key] = v.decode() if isinstance(v, bytes) else str(v)
                else:
                    rec[key] = ""
            rows.append(rec)
    out = os.path.join(RAW, "atlas_provenance.tsv")
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out} ({len(rows)} atlases)")


if __name__ == "__main__":
    main()
