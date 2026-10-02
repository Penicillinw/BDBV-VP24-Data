"""H2 line: resolve CellxGene asset URLs for the selected low-end atlases and
download the h5ad files into ./raw/atlas_h5ad/.

Selection (chosen for compartment coverage and modest size):
  3e524d45  Human adult ureter                         urothelium (adult, primary)
  8e5978bd  Adult bladder urothelial subset            urothelium (adult, snRNA)
  e6dad530  All cell types of human eye                RPE + rod + cone (adult)
  7b75b2c4  Photoreceptor cells of fovea/periphery     rod + cone (adult)
  cfa755c1  Supercluster: Choroid plexus               choroid plexus epithelium
  43dc52cb  Human embryonic male gonadal tissue         Sertoli (fetal caveat)
  2adb1f8a  Construction of a human cell landscape      Sertoli (adult, coarse labels)

Usage:  python h2_download_targets.py [--only ID ...] [--list]
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
DEST = RAW / "atlas_h5ad"
UA = {"User-Agent": "Mozilla/5.0 (research; BDBV-VP24 lowend-atlas-download)"}

TARGETS = {
    "3e524d45": "Human adult ureter",
    "8e5978bd": "Adult bladder urothelial subset",
    "e6dad530": "All cell types of human eye",
    "7b75b2c4": "Photoreceptor cells of the human fovea and peripheral retina",
    "cfa755c1": "Supercluster: Choroid plexus",
    "43dc52cb": "Human embryonic male gonadal tissue",
    "2adb1f8a": "Construction of a human cell landscape at single-cell level",
}


def safe_name(text: str) -> str:
    """Windows-safe file name component."""
    text = re.sub(r'[<>:"/\\|?*]', "_", text)
    return re.sub(r"\s+", "_", text).strip("_")


def resolve() -> dict[str, dict]:
    """Map short dataset-id prefix -> {title, url, size, dataset_id}."""
    out: dict[str, dict] = {}
    detail_dir = RAW / "cellxgene_collection_details"
    for fp in sorted(detail_dir.glob("*.json")):
        coll = json.loads(fp.read_text(encoding="utf-8"))
        for ds in coll.get("datasets") or []:
            did = ds.get("dataset_id") or ""
            key = did.split("-")[0]
            if key not in TARGETS:
                continue
            h5ad = [a for a in (ds.get("assets") or [])
                    if (a.get("filetype") or "").upper() == "H5AD"]
            if not h5ad:
                continue
            out[key] = {
                "dataset_id": did,
                "title": ds.get("title"),
                "collection": coll.get("name"),
                "doi": coll.get("doi"),
                "url": h5ad[0]["url"],
                "size": h5ad[0].get("filesize"),
                "cell_count": ds.get("cell_count"),
            }
    return out


def download(url: str, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
        total = int(r.headers.get("Content-Length") or 0)
        got = 0
        last = 0
        with tmp.open("wb") as fh:
            while True:
                chunk = r.read(1 << 22)
                if not chunk:
                    break
                fh.write(chunk)
                got += len(chunk)
                if total and got - last > (200 << 20):
                    last = got
                    print(f"    {dest.name}: {got/1e9:.2f}/{total/1e9:.2f} GB", flush=True)
    shutil.move(str(tmp), str(dest))
    print(f"    {dest.name}: done, {dest.stat().st_size/1e9:.2f} GB", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    DEST.mkdir(parents=True, exist_ok=True)
    info = resolve()
    (RAW / "download_manifest.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    if args.list:
        for k, v in info.items():
            print(f"{k}  {v['url'][:110]}  {v['size']}")
        return

    todo = args.only or list(TARGETS)
    for key in todo:
        if key not in info:
            print(f"!! no asset resolved for {key}", file=sys.stderr)
            continue
        rec = info[key]
        name = f"{key}_{safe_name(TARGETS[key])[:48]}.h5ad"
        dest = DEST / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"[skip] {name} already present ({dest.stat().st_size/1e9:.2f} GB)")
            continue
        print(f"[get ] {name}  size={rec['size']}  cells={rec['cell_count']}")
        download(rec["url"], dest)


if __name__ == "__main__":
    main()
