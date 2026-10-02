"""GSE342661: are the five biological replicates five donors? Read BioSample attributes,
and get SRA run sizes to judge feasibility."""

from __future__ import annotations

import re
from urllib.request import Request, urlopen


def get(url: str, tries: int = 3) -> str:
    last = None
    for a in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            with urlopen(req, timeout=180) as r:  # noqa: S310
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            import time

            time.sleep(2 * (a + 1))
    return f"ERROR: {last}"


def main() -> None:
    samples = ["SAMN62241143", "SAMN62241142", "SAMN62241141", "SAMN62241140"]
    for s in samples:
        t = get("https://www.ncbi.nlm.nih.gov/biosample/" + s + "?report=full&format=text")
        print(f"\n===== {s}")
        if t.startswith("ERROR"):
            print("  ", t[:120])
            continue
        for line in t.splitlines():
            if re.match(r"^/", line) or line.startswith(("Identifiers", "Organism", "Title",
                                                          "Description", "Attributes")):
                print("  ", line.strip()[:190])

    print("\n===== SRA run info")
    t = get("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1508721")
    if t.startswith("ERROR"):
        print("  ", t[:160])
    else:
        lines = [l for l in t.splitlines() if l.strip()]
        print("  rows:", len(lines))
        if lines:
            print("  header:", lines[0][:200])
            for l in lines[1:4]:
                print("  row:", l[:200])


if __name__ == "__main__":
    main()
