"""G6c: extract the GSE306664 RAW archive (131 per-sample 10x fixed-RNA h5 files)
and report the matrix layout actually present, before any biology is computed.
"""

from __future__ import annotations

import os
import sys
import tarfile

ROOT = r"G:\本迪布焦研究"
BASE = os.path.join(ROOT, "analysis", "g6c_primary_immune_20260926")
RAW = os.path.join(BASE, "raw")
TAR = os.path.join(RAW, "GSE306664_RAW.tar")
DEST = os.path.join(RAW, "gse306664_h5")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(DEST, exist_ok=True)
    with tarfile.open(TAR) as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(".h5")]
        print(f"{len(members)} h5 members")
        have = set(os.listdir(DEST))
        for i, member in enumerate(members, 1):
            name = os.path.basename(member.name)
            if name in have:
                continue
            tar.extract(member, path=DEST)
            src = os.path.join(DEST, member.name)
            if os.path.abspath(src) != os.path.abspath(os.path.join(DEST, name)):
                os.replace(src, os.path.join(DEST, name))
            if i % 20 == 0:
                print(f"  extracted {i}/{len(members)}", flush=True)
    print("extracted:", len(os.listdir(DEST)))

    import h5py  # imported late so extraction failures surface first

    sample = sorted(os.listdir(DEST))[0]
    with h5py.File(os.path.join(DEST, sample), "r") as fh:
        def walk(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"   {name}: shape={obj.shape} dtype={obj.dtype}")
        print(f"\nstructure of {sample}:")
        fh.visititems(walk)
        try:
            print("feature names:", [x.decode() for x in fh["features/name"][:5]])
        except Exception as exc:  # noqa: BLE001
            print("no features/name:", exc)


if __name__ == "__main__":
    main()
