"""Final availability check for GSE342661: per-GSM download links, FTP listing, plausible names."""

from __future__ import annotations

import re
from urllib.request import Request, urlopen

UA = {"User-Agent": "CodexResearch/1.0"}


def get(url: str) -> str:
    try:
        return urlopen(Request(url, headers=UA), timeout=120).read().decode("utf-8", "replace")
    except Exception as exc:  # noqa: BLE001
        return f"ERROR: {exc}"


def head(url: str) -> str:
    try:
        r = urlopen(Request(url, headers=UA, method="HEAD"), timeout=60)
        return f"OK {r.headers.get('Content-Length')} bytes"
    except Exception as exc:  # noqa: BLE001
        return f"none ({type(exc).__name__})"


def main() -> None:
    print("== FTP series dir")
    for path in ["", "suppl/", "matrix/", "soft/"]:
        t = get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE342nnn/GSE342661/" + path)
        files = [] if t.startswith("ERROR") else re.findall(r'href="([^"?/][^"]*)"', t)
        print(f"  {path or './'}: {files}")

    print("\n== plausible supplementary names")
    for name in ["GSE342661_counts.txt.gz", "GSE342661_TPM.txt.gz",
                 "GSE342661_processed_data.txt.gz", "GSE342661_RAW.tar",
                 "GSE342661_series_matrix.txt.gz"]:
        print("  suppl/" + name + ": "
              + head("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE342nnn/GSE342661/suppl/" + name))
        print("  matrix/" + name + ": "
              + head("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE342nnn/GSE342661/matrix/" + name))

    print("\n== GSM page (HTML) download links")
    html = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9935959")
    if html.startswith("ERROR"):
        print("  ", html[:160])
    else:
        links = re.findall(r'href="([^"]*(?:ftp|suppl|\.gz|\.tar)[^"]*)"', html)
        print("  candidate links:", sorted(set(links))[:12])


if __name__ == "__main__":
    main()
