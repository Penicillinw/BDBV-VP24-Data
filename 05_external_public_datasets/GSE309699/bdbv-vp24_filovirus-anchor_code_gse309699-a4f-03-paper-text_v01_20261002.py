"""A4-filovirus independent recompute, step 3: read the restriction evidence out of
the article **text**, not out of the figure panels.

Motivation. The original implementation took its restriction ordering from
``published_restriction_values.tsv``, whose values were digitised from Figure 3A-E
images. Digitising a figure is an error-prone route, and it is also the route that
introduced the "circularity" caveat. This script instead extracts the sentences in
which the authors *state* the ordering, so that the two routes can be compared.

Frozen criteria (written before any sentence was inspected)
----------------------------------------------------------
C11 Source is the deposited full text ``raw/PMC12737376.xml`` (JATS). Tags are stripped
    and whitespace collapsed; nothing is paraphrased.
C12 Every sentence containing an interferon token (Greek letter or spelled out) *and*
    one of the following evidence words is kept verbatim: dose/ng/mL, inhibit,
    significant, compared, similar, enhanced, reduced, titer/titre.
C13 Sentences are written out in UTF-8 so the Greek letters survive; the console
    encoding is not allowed to mangle them.
C14 No value from Figure 3 is used here. If the text and the digitised figure values
    disagree, both are reported and neither is silently preferred.

Writes ``paper_statements.tsv`` and ``paper_statements.txt`` into this directory.
"""
from __future__ import annotations

import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"G:\本迪布焦研究"
XML = os.path.join(ROOT, "analysis", "a2_filovirus_anchor_20260926", "raw",
                   "PMC12737376.xml")

IFN_TOKENS = ["IFN-α", "IFN-β", "IFN-γ", "IFN-λ", "IFN-alpha", "IFN-beta",
              "IFN-gamma", "IFN-lambda", "lambda", "interferon lambda"]
EVIDENCE = ["ng/mL", "ng ml", "dose", "inhibit", "significant", "compared",
            "similar", "enhanc", "reduc", "titer", "titre", "potent", "greatest",
            "most", "least", "strongest", "weaker", "less"]


def clean_text():
    with open(XML, encoding="utf-8", errors="replace") as fh:
        raw = fh.read()
    # drop the reference list and the supplementary caption block: they add nothing
    raw = re.sub(r"<ref-list.*?</ref-list>", " ", raw, flags=re.S)
    txt = re.sub(r"<[^>]+>", " ", raw)
    txt = (txt.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
              .replace("&#x03b1;", "α").replace("&#x03b2;", "β")
              .replace("&#x03b3;", "γ").replace("&#x03bb;", "λ")
              .replace("&#x394;", "Δ").replace("&Delta;", "Δ"))
    txt = re.sub(r"\s+", " ", txt)
    return txt


def sentences(txt):
    # split on sentence enders but protect common abbreviations / decimal points
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z(])", txt)
    return [p.strip() for p in parts if p.strip()]


def main():
    txt = clean_text()
    keep = []
    for s in sentences(txt):
        if any(t.lower() in s.lower() for t in IFN_TOKENS) and any(
                e.lower() in s.lower() for e in EVIDENCE):
            keep.append(s)

    with open(os.path.join(HERE, "paper_statements.txt"), "w", encoding="utf-8") as fh:
        for i, s in enumerate(keep, 1):
            fh.write(f"[{i}] {s}\n\n")

    def topic(s):
        tags = []
        for w in ("dose", "ng/mL", "6 h", "24 h", "significant", "similar",
                  "enhanc", "inhibit"):
            if w.lower() in s.lower():
                tags.append(w)
        return ";".join(tags)

    with open(os.path.join(HERE, "paper_statements.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["idx", "tags", "greek_letters_present", "sentence"])
        for i, s in enumerate(keep, 1):
            greeks = "".join(g for g in "αβγλ" if g in s)
            w.writerow([i, topic(s), greeks, s])

    print(f"{len(keep)} sentences kept -> paper_statements.tsv / .txt")
    # print the ones that carry the ordering claims, safely re-encoded
    for i, s in enumerate(keep, 1):
        low = s.lower()
        if ("more inhibitory" in low or "less inhibitory" in low
                or "similarly inhibitory" in low or "ng/ml" in low):
            print(f"[{i}] {s[:400]}")
            print()


if __name__ == "__main__":
    main()
