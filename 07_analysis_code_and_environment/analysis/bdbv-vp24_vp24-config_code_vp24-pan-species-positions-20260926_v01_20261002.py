"""Six-species Orthoebolavirus VP24 residue census at the four BDBV-defining positions.

For every ICTV-recognised Orthoebolavirus species, retrieves up to `PER_SPECIES_CAP`
VP24 protein records from NCBI, aligns each to the EBOV VP24 reference (NC_002549.1
CDS translation) with a Needleman-Wunsch global aligner, and reports the residues
occupying EBOV-equivalent positions 83, 135, 140, 141 (plus 184 and 217).

Purpose: test whether BDBV is the only species whose VP24 carries S/Q/H/A at those
four positions, and quantify how many sequences back that statement.

Read-only w.r.t. project data; writes only into
analysis/vp24_pan_species_20260925/coordinator_crosscheck/.
"""

from __future__ import annotations

import csv
import json
import os
import re
import time
from collections import Counter
from collections import defaultdict
from urllib.parse import quote
from urllib.error import URLError
from urllib.request import Request, urlopen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "vp24_pan_species_20260925", "coordinator_crosscheck")
EBOV_REF_ACC = "NC_002549.1"
POSITIONS = (83, 135, 140, 141, 184, 217)
FOUR = (83, 135, 140, 141)
PER_SPECIES_CAP = 120
SPECIES_CAP = {"EBOV": 10000, "SUDV": 10000, "BDBV": 10000}
BATCH = 500

SPECIES = {
    "BDBV": "Bundibugyo ebolavirus[Organism]",
    "EBOV": "Zaire ebolavirus[Organism]",
    "SUDV": "Sudan ebolavirus[Organism]",
    "TAFV": "Tai Forest ebolavirus[Organism]",
    "RESTV": "Reston ebolavirus[Organism]",
    "BOMV": "Bombali ebolavirus[Organism]",
}

# The NCBI `... AND VP24` query also returns chains from PDB entries whose
# *entry-level* annotation mentions VP24 while the chain itself is another
# protein (e.g. 8Y9J chains A/B are Nucleoprotein). Keep only VP24 records.
VP24_TITLE = re.compile(
    r"vp24|virus protein 24|viral protein 24|membrane-associated protein|"
    r"membrane associated protein|membrane-associated",
    re.IGNORECASE,
)
NOT_VP24_TITLE = re.compile(r"nucleoprotein|nucleocapsid|glycoprotein|polymerase", re.IGNORECASE)


def get(url: str, tries: int = 4, data: bytes | None = None) -> str:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "CodexResearch/1.0"})
            if data is not None:
                req.data = data
            with urlopen(req, timeout=300) as resp:  # noqa: S310 (fixed https host)
                return resp.read().decode("utf-8", "replace")
        except (URLError, OSError, ConnectionError) as exc:  # transient NCBI hiccups
            last = exc
            time.sleep(3 * (attempt + 1))
    raise SystemExit(f"NCBI request failed after {tries} attempts: {url}\n{last}")


def ebov_vp24() -> str:
    text = get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={EBOV_REF_ACC}&rettype=fasta_cds_aa&retmode=text"
    )
    for record in text.split(">")[1:]:
        lines = record.splitlines()
        if "[gene=VP24]" in lines[0]:
            return "".join(lines[1:]).strip().replace("*", "")
    raise SystemExit("EBOV VP24 CDS not found")


def species_vp24(code: str, term: str) -> tuple[str, list[tuple[str, str, str]], list[str]]:
    """Return (total_hits, [(accession, sequence, header), ...], excluded_accessions)."""
    cap = SPECIES_CAP.get(code, PER_SPECIES_CAP)
    query = quote(f"{term} AND VP24")
    search = json.loads(
        get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=protein&retmax={cap}&retmode=json&term={query}"
        )
    )["esearchresult"]
    total, ids = search["count"], search["idlist"]
    if not ids:
        return total, [], []

    out: list[tuple[str, str, str]] = []
    excluded: list[str] = []
    for start in range(0, len(ids), BATCH):
        chunk = ids[start : start + BATCH]
        body = ("db=protein&rettype=fasta&retmode=text&id=" + ",".join(chunk)).encode()
        text = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", data=body)
        for record in text.split(">")[1:]:
            lines = record.splitlines()
            header = lines[0]
            acc = header.split()[0]
            seq = "".join(lines[1:]).strip().replace("*", "")
            if not seq:
                continue
            if NOT_VP24_TITLE.search(header) or not VP24_TITLE.search(header):
                excluded.append(acc)
                continue
            out.append((acc, seq, header))
        print(f"    ... {code}: fetched {min(start + BATCH, len(ids))}/{len(ids)}")
    return total, out, excluded


def needle(a: str, b: str) -> list[tuple[str, str]]:
    """Global alignment (a = reference, b = query) with linear gap costs."""
    n, m = len(a), len(b)
    gap = -6

    def sub(x: str, y: str) -> int:
        return 2 if x == y else -1

    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * gap
    for j in range(1, m + 1):
        score[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score[i][j] = max(
                score[i - 1][j - 1] + sub(a[i - 1], b[j - 1]),
                score[i - 1][j] + gap,
                score[i][j - 1] + gap,
            )
    i, j = n, m
    pairs: list[tuple[str, str]] = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i - 1][j - 1] + sub(a[i - 1], b[j - 1]):
            pairs.append((a[i - 1], b[j - 1]))
            i, j = i - 1, j - 1
        elif i > 0 and score[i][j] == score[i - 1][j] + gap:
            pairs.append((a[i - 1], "-"))
            i -= 1
        else:
            pairs.append(("-", b[j - 1]))
            j -= 1
    return pairs[::-1]


def positions_of(ref: str, seq: str) -> dict[int, str]:
    if len(seq) == len(ref):
        # 251-aa orthologues align 1:1 with the EBOV reference; verify by identity
        ident = sum(1 for a, b in zip(ref, seq) if a == b) / len(ref)
        if ident >= 0.5:
            return {p: seq[p - 1] for p in POSITIONS}
    pairs = needle(ref, seq)
    out: dict[int, str] = {}
    ref_idx = 0
    for x, y in pairs:
        if x == "-":
            continue
        ref_idx += 1
        if ref_idx in POSITIONS:
            out[ref_idx] = y
    return out


_CACHE: dict[str, dict[int, str]] = {}


def cached_positions(ref: str, seq: str) -> dict[int, str]:
    """Alignment is the expensive step; identical sequences are aligned once."""
    if seq not in _CACHE:
        _CACHE[seq] = positions_of(ref, seq)
    return _CACHE[seq]


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    ref = ebov_vp24()
    print(f"EBOV VP24 reference ({EBOV_REF_ACC}): {len(ref)} aa")

    rows: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "ebov_reference": EBOV_REF_ACC,
        "ebov_vp24_length": len(ref),
        "per_species_cap": PER_SPECIES_CAP,
        "per_species_cap_overrides": SPECIES_CAP,
        "positions": list(POSITIONS),
        "four_position_motif_definition": "83/135/140/141 in EBOV VP24 numbering",
        "species": {},
    }

    for code, term in SPECIES.items():
        total, seqs, excluded = species_vp24(code, term)
        motifs: Counter[str] = Counter()
        motif_examples: dict[str, list[str]] = defaultdict(list)
        per_position: dict[int, Counter[str]] = {p: Counter() for p in POSITIONS}
        for acc, seq, header in seqs:
            if len(seq) < 200:
                # partial fragments cannot support a four-position call
                rows.append(
                    {
                        "species": code,
                        "accession": acc,
                        "length_aa": len(seq),
                        **{f"pos{p}": "" for p in POSITIONS},
                        "four_position_motif": "NA_partial",
                        "header": header[:120],
                    }
                )
                continue
            pos = cached_positions(ref, seq)
            motif = "".join(pos.get(p, "?") for p in FOUR)
            motifs[motif] += 1
            motif_examples[motif].append(acc)
            for p in POSITIONS:
                per_position[p][pos.get(p, "?")] += 1
            rows.append(
                {
                    "species": code,
                    "accession": acc,
                    "length_aa": len(seq),
                    **{f"pos{p}": pos.get(p, "?") for p in POSITIONS},
                    "four_position_motif": motif,
                    "header": header[:120],
                }
            )
        summary["species"][code] = {
            "esearch_hits": int(total),
            "records_returned": len(seqs),
            "sequences_assessed": sum(1 for _, s, _ in seqs if len(s) >= 200),
            "n_full_length_251": sum(1 for _, s, _ in seqs if len(s) == 251),
            "n_partial_skipped": sum(1 for _, s, _ in seqs if len(s) < 200),
            "n_excluded_non_vp24_chain": len(excluded),
            "excluded_non_vp24_examples": excluded[:8],
            "four_position_motif_counts": dict(motifs.most_common()),
            "non_canonical_motif_examples": {
                m: a[:6] for m, a in motif_examples.items() if m != "SQHA" and len(motifs) > 1
            },
            "per_position_residue_counts": {
                str(p): dict(per_position[p].most_common()) for p in POSITIONS
            },
            "carries_SQHA": "SQHA" in motifs,
        }
        print(f"{code}: hits={total} assessed={len(seqs)} motifs={dict(motifs.most_common(4))}")
        with open(
            os.path.join(OUT, "vp24_six_species_positions.tsv"), "w", encoding="utf-8", newline=""
        ) as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    summary["species_with_SQHA"] = [
        c for c, d in summary["species"].items() if d["carries_SQHA"]
    ]
    summary["claim_bdbv_exclusive_supported"] = summary["species_with_SQHA"] == ["BDBV"]

    with open(os.path.join(OUT, "vp24_six_species_positions.tsv"), "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with open(os.path.join(OUT, "vp24_six_species_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    print()
    print(json.dumps(summary["species_with_SQHA"], ensure_ascii=False))
    print("exclusive claim supported:", summary["claim_bdbv_exclusive_supported"])


if __name__ == "__main__":
    main()
