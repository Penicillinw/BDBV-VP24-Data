"""B1 positive control - score the FoldX mutation cycle against PMID 39894818.

Pre-registered sign predictions (ddG_int = mutant - paired wild type,
positive = weaker binding), taken verbatim from the Europe PMC abstract of
PMID 39894818 (Commun Biol 2025, DOI 10.1038/s42003-025-07613-y):
SUDV binds hNPC1 more strongly than EBOV; Ile79 / Ala141 / Pro148 enhance
binding, Gln142 reduces it.

  SUDV -> EBOV : I79V > 0, A141V > 0, P148A > 0, Q142S < 0
  EBOV -> SUDV : V79I < 0, V141A < 0, A148P < 0, S142Q > 0

Reads analysis/b1_positive_control_20260924/foldx_mutation_results.tsv and
writes .../b1_sign_test.tsv. Reader-only: no FoldX call, no file mutation.
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

BASE = Path(r"G:\本迪布焦研究") / "analysis" / "b1_positive_control_20260924"

PREDICTED = {
    ("sudv", "I79V"): +1,
    ("sudv", "A141V"): +1,
    ("sudv", "P148A"): +1,
    ("sudv", "Q142S"): -1,
    ("sudv_copy2", "I79V"): +1,
    ("sudv_copy2", "A141V"): +1,
    ("sudv_copy2", "P148A"): +1,
    ("sudv_copy2", "Q142S"): -1,
    ("ebov", "V79I"): -1,
    ("ebov", "V141A"): -1,
    ("ebov", "A148P"): -1,
    ("ebov", "S142Q"): +1,
}


def main() -> None:
    src = BASE / "foldx_mutation_results.tsv"
    with src.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    agg: dict[tuple[str, str], list[float]] = {}
    for r in rows:
        agg.setdefault((r["complex"], r["mutation"]), []).append(float(r["ddG_int"]))

    out_rows = []
    for key, vals in sorted(agg.items()):
        mean = statistics.mean(vals)
        sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
        exp = PREDICTED[key]
        ok = (mean > 0) == (exp > 0)
        out_rows.append(
            {
                "complex": key[0],
                "mutation": key[1],
                "n_runs": len(vals),
                "mean_ddG_int": round(mean, 3),
                "sd_ddG_int": round(sd, 3),
                "expected_sign": "+" if exp > 0 else "-",
                "observed_sign": "+" if mean > 0 else "-",
                "sign_match": ok,
                "magnitude_over_noise": round(abs(mean) / sd, 2) if sd > 0 else None,
            }
        )

    out = BASE / "b1_sign_test.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    hdr = "{:<11} {:<7} {:>8} {:>7} {:>8} {:>8} {:>6}".format(
        "complex", "mut", "mean", "sd", "expect", "obs", "match"
    )
    print(hdr)
    for r in out_rows:
        print(
            "{:<11} {:<7} {:>+8.3f} {:>7.3f} {:>8} {:>8} {:>6}".format(
                r["complex"],
                r["mutation"],
                r["mean_ddG_int"],
                r["sd_ddG_int"],
                r["expected_sign"],
                r["observed_sign"],
                str(r["sign_match"]),
            )
        )

    n_ok = sum(1 for r in out_rows if r["sign_match"])
    print("sign matches: {}/{}".format(n_ok, len(out_rows)))

    per_position = {}
    for r in out_rows:
        pos = r["mutation"][1:-1]
        per_position.setdefault(pos, []).append(r["sign_match"])
    for pos, oks in sorted(per_position.items(), key=lambda x: int(x[0])):
        print("position {}: {}/{} directions match".format(pos, sum(oks), len(oks)))
    print("wrote", out)


if __name__ == "__main__":
    main()
