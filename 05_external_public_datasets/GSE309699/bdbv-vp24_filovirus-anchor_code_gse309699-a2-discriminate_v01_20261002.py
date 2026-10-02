"""Independent discriminator for the GSE309699 column map.

The digit/time rule failed (every group has a monotone digit trend), so the map is
settled with a signature that is orthogonal to overall ISG level: the IFN-gamma
signature (CXCL9/10/11, GBP1/2/5, IDO1, IRF1). Whatever group carries that signature is
IFN-gamma regardless of how high its general ISG panel is.
"""
from __future__ import annotations

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
GAMMA = ["CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP2", "GBP5", "IDO1", "IRF1"]
TYPEI = ["ISG15", "MX1", "RSAD2", "IFIT1", "IFIT3", "IFI6", "USP18", "OAS1"]


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "module_levels.tsv"),
                                    encoding="utf-8"), delimiter="\t"))
    sel = [r for r in rows if r["block"] == "all"]
    gs = GAMMA + TYPEI + ["STAT1", "IFNLR1", "IFNAR1", "IFNAR2"]
    print("group  n  " + " ".join(f"{g[:7]:>7s}" for g in gs))
    for r in sel:
        print(f"{r['group']:5s} {r['n']:>2s}  "
              + " ".join(f"{float(r[g]):7.2f}" for g in gs))

    print("\ngamma-signature mean minus type-I mean (higher = more gamma-like):")
    out = []
    for r in sel:
        gm = sum(float(r[g]) for g in GAMMA) / len(GAMMA)
        tm = sum(float(r[g]) for g in TYPEI) / len(TYPEI)
        out.append((r["group"], r["n"], gm, tm, gm - tm))
    for g, n, gm, tm, d in sorted(out, key=lambda x: -x[4]):
        print(f"  {g} (n={n}): gamma={gm:.2f} typeI={tm:.2f} diff={d:+.2f}")

    print("\nper time block (digit<=4 assumed 6 h, digit>=5 assumed 24 h):")
    by = {(r["group"], r["block"]): r for r in rows}
    print(f"{'group':6s}{'block':>10s}{'n':>3s}{'panel':>8s}{'ISG15':>8s}"
          f"{'MX1':>8s}{'RSAD2':>8s}{'IFIT1':>8s}{'CXCL10':>8s}{'GBP1':>8s}")
    for g in ("A", "B", "M", "G", "K", "N"):
        for blk in ("digit<=4", "digit>=5"):
            if (g, blk) not in by:
                continue
            r = by[(g, blk)]
            print(f"{g:6s}{blk:>10s}{r['n']:>3s}"
                  + "".join(f"{float(r[k]):8.2f}" for k in
                            ("panel_mean", "ISG15", "MX1", "RSAD2", "IFIT1", "CXCL10", "GBP1")))


if __name__ == "__main__":
    main()
