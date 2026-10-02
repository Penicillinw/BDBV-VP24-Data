"""G6c helper: resumable download of large GEO supplementary files via urllib
(the local curl/schannel stack fails TLS against ftp.ncbi.nlm.nih.gov).

Usage: python tools/g6c_download_20260926.py <url> <dest> [expected_bytes]
"""

from __future__ import annotations

import os
import sys
import time
from urllib.request import Request, urlopen


def main() -> None:
    url, dest = sys.argv[1], sys.argv[2]
    expected = int(sys.argv[3]) if len(sys.argv) > 3 else None
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    for attempt in range(1, 9):
        have = os.path.getsize(dest) if os.path.exists(dest) else 0
        if expected is not None and have >= expected:
            print(f"complete: {have} bytes")
            return
        headers = {"User-Agent": "CodexResearch/1.0"}
        if have:
            headers["Range"] = f"bytes={have}-"
        try:
            req = Request(url, headers=headers)
            with urlopen(req, timeout=180) as resp:  # noqa: S310
                mode = "ab" if have and resp.status == 206 else "wb"
                if mode == "wb":
                    have = 0
                total = resp.headers.get("Content-Length")
                total = int(total) + have if total else expected
                start = time.time()
                last = start
                with open(dest, mode) as fh:
                    while True:
                        chunk = resp.read(1 << 20)
                        if not chunk:
                            break
                        fh.write(chunk)
                        have += len(chunk)
                        now = time.time()
                        if now - last > 15:
                            rate = have / max(now - start, 1e-6) / 1e6
                            pct = f"{100 * have / total:.1f}%" if total else "?"
                            print(f"{have:,} bytes ({pct}) {rate:.1f} MB/s", flush=True)
                            last = now
                print(f"done pass: {have:,} bytes", flush=True)
                if expected is None or have >= expected:
                    return
        except Exception as exc:  # noqa: BLE001 - transient network failures
            print(f"attempt {attempt} failed: {exc}", flush=True)
            time.sleep(5 * attempt)
    raise SystemExit("download incomplete after retries")


if __name__ == "__main__":
    main()
