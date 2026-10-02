"""GSE342661: read the series matrix header and tabulate the full design.

Frozen criterion: the dataset counts as column-mappable only if every sample
carries explicit treatment / time / route fields in the series matrix header
(so no hole-name-to-condition guesswork is needed).

Read-only apart from raw/.
"""

from __future__ import annotations

import collections
import gzip
import io
import os
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
os.makedirs(RAW, exist_ok=True)
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE342nnn/GSE342661/matrix/"
       "GSE342661_series_matrix.txt.gz")

req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
blob = urllib.request.urlopen(req, timeout=180).read()
with open(os.path.join(RAW, "GSE342661_series_matrix.txt.gz"), "wb") as fh:
    fh.write(blob)
print("series matrix bytes:", len(blob))

txt = gzip.decompress(blob).decode("utf-8", "replace")
lines = txt.splitlines()
print("matrix header lines:", sum(1 for l in lines if l.startswith("!Sample")))
print("data rows:", sum(1 for l in lines if l and not l.startswith("!") and "\t" in l))

meta = collections.defaultdict(list)
for ln in lines:
    if not ln.startswith("!Sample"):
        continue
    parts = [p.strip('"') for p in ln.split("\t")]
    meta[parts[0].lstrip("!")] = parts[1:]

titles = meta.get("Sample_title", [])
chars = meta.get("Sample_characteristics_ch1", [])
print("samples:", len(titles))

# characteristics rows are repeated per distinct key; collect them
keys = collections.Counter()
per_key = collections.defaultdict(list)
for row in chars:
    key = row.split(":")[0].strip()
    keys[key] += 1
    per_key[key].append(row.split(":", 1)[1].strip() if ":" in row else "")
print("characteristics keys:", dict(keys))

for key in ("treatment", "time", "route"):
    vals = per_key.get(key)
    if vals is None:
        print(f"  {key}: ABSENT")
        continue
    print(f"  {key}: {dict(sorted(collections.Counter(vals).items()))}")

# full design grid
combos = collections.Counter()
t, tm, rt = per_key.get("treatment"), per_key.get("time"), per_key.get("route")
if t and tm and rt and len(t) == len(tm) == len(rt) == len(titles):
    for a, b, c in zip(t, tm, rt):
        combos[(a, b, c)] += 1
    print("\ndesign grid (treatment / time / route -> n):")
    for k in sorted(combos):
        print("   %-10s %-5s %-14s %d" % (k[0], k[1], k[2], combos[k]))
    print("cells with n<3:", [k for k, v in combos.items() if v < 3])
else:
    print("\n!! characteristics do not carry treatment/time; the sample TITLES do.")
    import re
    pat = re.compile(r"^(\d+h),\s*(IFN[A-Za-z0-9]+),\s*biol rep (\d+),\s*(a|ab)$")
    parsed = [(m.group(1), m.group(2), int(m.group(3)), m.group(4))
              for m in (pat.match(x) for x in titles) if m]
    print("titles parsed :", len(parsed), "of", len(titles))
    unparsed = [x for x in titles if not pat.match(x)]
    if unparsed:
        print("unparsed:", unparsed[:5])
    route_of = {"a": "apical-only", "ab": "apical-basal"}
    grid = collections.Counter((tm_, tr_, route_of[rt_]) for tm_, tr_, _r, rt_ in parsed)
    print("\ndesign grid (time / treatment / route -> n):")
    for k in sorted(grid, key=lambda x: (x[0], x[1], x[2])):
        print("   %-5s %-8s %-14s %d" % (k[0], k[1], k[2], grid[k]))
    print("cells:", len(grid), " cells with n<3:", [k for k, v in grid.items() if v < 3])
    reps = collections.Counter((tm_, tr_, rt_) for tm_, tr_, _r, rt_ in parsed)
    print("max replicate index per cell:", max(r for *_x, r in
                                               [(a, b, c, d) for a, b, c, d in parsed]),
          "(replicate numbers 1..N within each cell)")
    print("unique replicates per cell:",
          sorted({r for *_x, r in parsed}))
