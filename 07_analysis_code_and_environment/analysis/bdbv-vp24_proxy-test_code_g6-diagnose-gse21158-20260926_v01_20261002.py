"""G6 diagnostic: is the baseline-vs-induction correlation driven by a bimodal
responder / non-responder split rather than by a graded physiological relationship?

A bimodal split in cancer cell lines is the signature expected if some lines carry
lesions in the interferon pathway (JAK/STAT/IFNAR) rather than if baseline tone
physiologically tunes responsiveness. This decides whether GSE21158 is usable
material for the claim at all.
"""

from __future__ import annotations

import json
import math
import os

ROOT = r"G:\本迪布焦研究"
OUT = os.path.join(ROOT, "analysis", "g6_ifn_responsiveness_20260926")


def main() -> None:
    rel = json.load(open(os.path.join(OUT, "g6_reliability.json"), encoding="utf-8"))
    lines = rel["cell_lines"]
    base = rel["baseline_values"]
    ind = rel["induction24_values"]
    rows = sorted(zip(lines, base, ind), key=lambda t: t[2])

    responders = [t for t in rows if t[2] >= 0.3]
    non = [t for t in rows if t[2] < 0.3]

    def mean(xs):
        return sum(xs) / len(xs) if xs else float("nan")

    out = {
        "n_lines": len(rows),
        "induction_sorted": [[l, round(b, 3), round(i, 3)] for l, b, i in rows],
        "responders_threshold_0.3": {
            "n": len(responders),
            "lines": [t[0] for t in responders],
            "mean_baseline": round(mean([t[1] for t in responders]), 3),
            "mean_induction": round(mean([t[2] for t in responders]), 3),
        },
        "non_responders_threshold_0.3": {
            "n": len(non),
            "lines": [t[0] for t in non],
            "mean_baseline": round(mean([t[1] for t in non]), 3),
            "mean_induction": round(mean([t[2] for t in non]), 3),
        },
        "gap_in_induction": round(rows[len(non)][2] - rows[len(non) - 1][2], 3) if non else None,
        "phenotype_is_binary": bool(non and rows[len(non)][2] > 5 * rows[len(non) - 1][2]),
    }

    def spearman(a, b):
        def rank(x):
            order = sorted(range(len(x)), key=lambda i: x[i])
            r = [0.0] * len(x)
            i = 0
            while i < len(order):
                j = i
                while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                    j += 1
                avg = (i + j) / 2 + 1
                for k in range(i, j + 1):
                    r[order[k]] = avg
                i = j + 1
            return r
        ra, rb = rank(a), rank(b)
        n = len(a)
        ma, mb = sum(ra) / n, sum(rb) / n
        num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
        da = math.sqrt(sum((v - ma) ** 2 for v in ra))
        db = math.sqrt(sum((v - mb) ** 2 for v in rb))
        return num / (da * db) if da and db else float("nan")

    if responders and non:
        r_base = [t[1] for t in responders]
        r_ind = [t[2] for t in responders]
        n_base = [t[1] for t in non]
        n_ind = [t[2] for t in non]
        out["sign_flip_test"] = {
            "spearman_all_10": round(spearman([t[1] for t in rows], [t[2] for t in rows]), 3),
            "spearman_responders_only_n5": round(spearman(r_base, r_ind), 3),
            "spearman_non_responders_only_n5": round(spearman(n_base, n_ind), 3),
            "interpretation": "if the sign differs between strata, the pooled correlation "
                              "reflects the bimodal split rather than a graded relationship",
        }
    with open(os.path.join(OUT, "g6_bimodality_diagnostic.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
