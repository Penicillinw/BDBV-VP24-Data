"""SRB step 0 - fetch the independent adult spleen atlas (Tabula Sapiens - Spleen).

Frozen purpose (see FINDINGS.md R7): the rebuilt spleen reference is estimated
from the local Human Cell Landscape (HCL) file; the *test* needs a spleen atlas
from a different study.  The chosen test resource is

    CELLxGENE dataset cee11228-9f0b-4e57-afe2-cfe15ee56312
    title "Tabula Sapiens - Spleen", 70,448 cells, donors TSP2/TSP7/TSP14/TSP25...
    h5ad 2,531.7 MB
    https://datasets.cellxgene.cziscience.com/5498457f-6a18-4766-b9a3-a2c8c4662703.h5ad

Independence: different study, different platform mix (10x 3' v3 / 10x 5' v2 /
Smart-seq2) and different donors from HCL (microwell-seq).  It shares only the
public data portal (CELLxGENE) with the atlas set used earlier, which is not a
source of dependence between reference and test.

Writes only into this task directory.  Safe to re-run: resumes with HTTP Range.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "analysis", "spleen_reference_rebuild_20260926")
RAW = os.path.join(OUT, "raw")

URL = "https://datasets.cellxgene.cziscience.com/5498457f-6a18-4766-b9a3-a2c8c4662703.h5ad"
DEST = os.path.join(RAW, "TS_spleen.h5ad")
EXPECT = 2531.73 * 1024 * 1024  # nominal, from the CELLxGENE asset listing
UA = {"User-Agent": "Mozilla/5.0 (research; BDBV-VP24 spleen-reference-rebuild)"}


def head_size() -> int | None:
    req = urllib.request.Request(URL, headers=UA, method="HEAD")
    with urllib.request.urlopen(req, timeout=120) as r:
        n = r.headers.get("Content-Length")
    return int(n) if n else None


def download() -> None:
    os.makedirs(RAW, exist_ok=True)
    total = head_size()
    have = os.path.getsize(DEST) if os.path.exists(DEST) else 0
    if total is not None and have == total:
        print(f"already complete: {have} bytes", flush=True)
        return
    headers = dict(UA)
    if have:
        headers["Range"] = f"bytes={have}-"
    req = urllib.request.Request(URL, headers=headers)
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=600) as r, open(DEST, "ab") as fh:
        got = have
        while True:
            chunk = r.read(1 << 22)
            if not chunk:
                break
            fh.write(chunk)
            got += len(chunk)
            if total and got % (1 << 28) < (1 << 22):
                pct = 100.0 * got / total
                rate = got / max(1e-9, time.time() - t0) / 1e6
                print(f"  {got/1e9:6.3f} GB / {total/1e9:.3f} GB  ({pct:5.1f}%)  {rate:6.1f} MB/s", flush=True)
    print(f"downloaded {os.path.getsize(DEST)} bytes to {DEST}", flush=True)


def sha16(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()[:16]


def main() -> None:
    download()
    size = os.path.getsize(DEST)
    info = {"url": URL, "dest": os.path.relpath(DEST, ROOT), "bytes": size,
            "sha256_first16": sha16(DEST),
            "cellxgene_dataset_id": "cee11228-9f0b-4e57-afe2-cfe15ee56312",
            "title": "Tabula Sapiens - Spleen", "cell_count": 70448,
            "retrieved": time.strftime("%Y-%m-%d")}
    with open(os.path.join(RAW, "ts_spleen_download.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=1)
    print(json.dumps(info, indent=1), flush=True)


if __name__ == "__main__":
    main()
