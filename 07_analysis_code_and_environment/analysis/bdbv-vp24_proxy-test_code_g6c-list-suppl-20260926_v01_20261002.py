"""G6c helper: list GEO supplementary files (name + size) for the accessions
given on the command line. Usage:

    python tools/g6c_list_suppl_20260926.py GSE327707 GSE72502
"""

from __future__ import annotations

import re
import sys
from urllib.request import Request, urlopen


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    for acc in sys.argv[1:]:
        prefix = acc[:-3] + "nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/suppl/"
        try:
            html = urlopen(  # noqa: S310
                Request(url, headers={"User-Agent": "CodexResearch/1.0"}), timeout=120
            ).read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            print(f"== {acc}: ERROR {exc}")
            continue
        print(f"== {acc}")
        for line in html.splitlines():
            if 'href="' not in line or "Parent Directory" in line:
                continue
            name = line.split('href="', 1)[1].split('"', 1)[0]
            if name.startswith("/") or name.startswith("http"):
                continue
            size = re.findall(r"([0-9.]+[KMG]?)\s*$", line.strip())
            print(f"   {name}  {size[0] if size else ''}")
        listing = html
        for m in re.finditer(r"(?m)^<a href=\"([^\"?/][^\"]*)\">([^<]+)</a>", listing):
            pass


if __name__ == "__main__":
    main()
