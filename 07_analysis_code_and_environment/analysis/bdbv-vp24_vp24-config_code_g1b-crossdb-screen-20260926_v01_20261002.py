"""G1b: repair the two G1 retrieval routes that produced no sequence-level evidence.

`tools/g1_crossdb_census_20260926.py` (run 2026-09-26T01:18-01:21) left two holes,
both of which are *failures*, not negative results:

* **ENA / INSDC** -- the portal rejected ``fields=accession,description,sequence_length``
  with HTTP 400 for ``result=sequence``, so every ENA call returned an empty list and the
  residue-level screen never ran (record *counts* were still obtained through a fallback).
* **BV-BRC** -- the G1 payload carries only ``aa_sequence_md5``, so G1 could only ask
  "does BV-BRC hold a protein identical to <this species' own UniProt backbone with
  S/Q/H/A engineered in>".  That is a narrow proxy; the actual positional screen of
  BV-BRC proteins was never done.

This script fixes both and re-runs them.

* ENA: enumerate records with the accepted field list, prefer records whose description
  says "complete genome", fetch the EMBL flat file, pull the VP24 CDS translation, align
  it to the EBOV reference and call positions 83/135/140/141 by alignment.
  **Sampling is explicit and reported** (large archives are capped, small ones are complete).
* BV-BRC: screen *every* distinct VP24 ``aa_sequence_md5`` position by position.  BV-BRC's
  ``genome_feature`` row carries the md5 but **not** the sequence, and the
  ``protein``/``protein_sequence`` data types answered HTTP 404 on 2026-09-26 (API v1.9.3).
  BV-BRC ``protein_id`` values are GenBank/RefSeq accessions, so the sequence text is
  recovered from NCBI E-utilities ``efetch(db=protein)`` and the md5 is re-verified.

Read-only w.r.t. published project data; writes only into
``analysis/g1b_crossdb_screen_20260926/`` (``raw/`` keeps every payload + sha256).
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "g1b_crossdb_screen_20260926")
RAW = os.path.join(OUT, "raw")
G1_RAW = os.path.join(ROOT, "analysis", "g1_crossdb_census_20260926", "raw")
UA = {"User-Agent": "bdbv-vp24-census/1.1 (research script)"}

# EBOV VP24 reference (Zaire ebolavirus Mayinga-76), UniProt P0C779.
EBOV_VP24 = (
    "MAKATGRYNLISPKKDLEKGVVLSDLCNFLVSQTIQGWKVYWAGIEFDVTHKGMALLHRL"
    "KTNDFAPAWSMTRNLFPHLFQNPNSTIESPLWALRVILAAGIQDQLIDQSLIEPLAGALG"
    "LISDWLLTTNTNHFNMRTQRVKEQLSLKMLSLIRSNILKFINKLDALHVVNYNGLLSSIE"
    "IGTQNHTIIITRTNMGFLVELQEPDKSAMNRMKPGPAKFSLLHESTLKAFTQGSSTRMQS"
    "LILEFNSSLAI"
)
PRIMARY = (83, 135, 140, 141)
BDBV_TYPE = "SQHA"

SPECIES = {
    "BDBV": 565995,
    "EBOV": 186538,
    "SUDV": 186540,
    "RESTV": 186539,
    "TAFV": 186541,
    "BOMV": 2010960,
}
# Records kept per species on the ENA route.  BDBV/TAFV/BOMV are complete; the large
# archives (EBOV 3,680 / SUDV 164 / RESTV 61) are capped and the cap is reported.
ENA_CAP = {"BDBV": 200, "EBOV": 60, "SUDV": 60, "RESTV": 60, "TAFV": 200, "BOMV": 200}

LOG: list[str] = []
# `--reuse-raw` reads payloads already archived under ``analysis/g1b_crossdb_screen_20260926/raw``
# instead of re-downloading them.  Used to re-score an archive after a parser fix without
# hammering the archives again; every reused payload keeps the sha256 of the original download.
REUSE_RAW = False


def say(msg: str = "") -> None:
    LOG.append(msg)
    print(msg, flush=True)


def fetch(url: str, tag: str, tries: int = 3, sleep_after: float = 0.0) -> bytes:
    """GET with retries; the raw payload and its sha256 are archived."""
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, tag)
    if REUSE_RAW and os.path.exists(path):
        with open(path, "rb") as fh:
            body = fh.read()
        say("  [cache] %-46s %9d bytes" % (tag, len(body)))
        return body
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as r:
                body = r.read()
            with open(path, "wb") as fh:
                fh.write(body)
            say("  [fetch] %-46s %9d bytes  sha256:%s"
                % (tag, len(body), hashlib.sha256(body).hexdigest()[:16]))
            if sleep_after:
                time.sleep(sleep_after)
            return body
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1.5 + 2.5 * attempt)
    raise RuntimeError("fetch failed for %s: %s" % (url, last))


def md5_of(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()


_ALIGNER = None


def align_positions(seq: str, positions=PRIMARY, reference: str = EBOV_VP24) -> dict:
    """Map reference numbering onto `seq` by global alignment (BLOSUM62, affine gaps)."""
    global _ALIGNER
    if _ALIGNER is None:
        a = PairwiseAligner()
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -11
        a.extend_gap_score = -1
        _ALIGNER = a
    try:
        aln = _ALIGNER.align(reference, seq)[0]
    except Exception:  # noqa: BLE001
        return {p: None for p in positions}
    a_str, b_str = str(aln[0]), str(aln[1])
    ref_i = 0
    out: dict[int, str | None] = {}
    for ca, cb in zip(a_str, b_str):
        if ca != "-":
            ref_i += 1
        if ca != "-" and ref_i in positions:
            out[ref_i] = cb if cb != "-" else None
    for p in positions:
        out.setdefault(p, None)
    return out


def calls_key(calls: dict) -> str:
    return "".join(calls[p] or "-" for p in PRIMARY)


def is_bdbv_type(calls: dict) -> bool:
    return all(calls.get(p) == c for p, c in zip(PRIMARY, BDBV_TYPE))


# --------------------------------------------------------------------- ENA / INSDC
def ena_list(taxid: int, tag: str, limit: int = 200000) -> list[dict]:
    """Record table for a taxon.  `sequence_length` is rejected (HTTP 400) for
    ``result=sequence``; the accepted field list is accession + description."""
    params = {
        "result": "sequence",
        "query": "tax_eq(%d)" % taxid,
        "fields": "accession,description",
        "format": "tsv",
        "limit": str(limit),
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(params)
    body = fetch(url, tag)
    lines = body.decode("utf-8", "replace").strip().splitlines()
    if len(lines) <= 1:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln.strip()]


def embl_feature_blocks(text: str) -> list[tuple[str, str]]:
    """[(feature_key, block_text)] parsed from an EMBL flat file.

    EMBL puts the feature key in columns 6-20 and qualifiers (and the wrapped
    continuation lines of a long qualifier value such as ``/translation``) from
    column 22 on.  A line whose columns 6-20 are blank is therefore a **qualifier
    continuation**, not a new feature: treating it as a key truncates every
    translation that wraps onto more than one line, which is what made the first
    G1b run report "VP24 CDS called = 0" for all six species.
    """
    blocks: list[tuple[str, str]] = []
    key: str | None = None
    buf: list[str] = []
    for line in text.splitlines():
        if not line.startswith("FT"):
            continue
        body = line[2:]
        if not body.strip():
            continue
        # Key field = flat-file columns 6-20 -> body indices 3..18; qualifiers and
        # their wrapped continuations start at flat-file column 22 -> body index 19.
        keyfield = body[3:19].strip()
        if not keyfield:
            if key is not None:
                buf.append(body)
            continue
        if key is not None:
            blocks.append((key, "\n".join(buf)))
        key = keyfield.split()[0]
        buf = [body]
    if key is not None:
        blocks.append((key, "\n".join(buf)))
    return blocks


def ena_vp24_translations(accession: str, tag: str) -> list[tuple[str, str, str]]:
    """[(gene, product, translation)] for VP24 CDS features of an EMBL record."""
    url = "https://www.ebi.ac.uk/ena/browser/api/embl/%s?lineLimit=0" % accession
    body = fetch(url, tag, sleep_after=0.25)
    text = body.decode("utf-8", "replace")
    blocks = embl_feature_blocks(text)

    out: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    for fkey, chunk in blocks:
        if fkey != "CDS":
            continue
        gene = re.search(r'/gene="([^"]+)"', chunk)
        product = re.search(r'/product="([^"]+)"', chunk)
        trans = re.search(r'/translation="([^"]+)"', chunk, re.S)
        gname = gene.group(1) if gene else ""
        pname = product.group(1) if product else ""
        if not trans:
            continue
        if "VP24" not in gname.upper() and "VP24" not in pname.upper():
            continue
        seq = re.sub(r"\s+", "", trans.group(1))
        if seq in seen:
            continue
        seen.add(seq)
        out.append((gname, pname, seq))
    return out


def route_ena() -> dict:
    say("=" * 78)
    say("EN  ENA / INSDC -- VP24 CDS translations parsed from EMBL flat files")
    say("    (field list repaired: the G1 run sent `sequence_length`, rejected HTTP 400)")
    result: dict[str, dict] = {}
    for sp, taxid in SPECIES.items():
        recs = ena_list(taxid, "g1b_ena_%s_list.tsv" % sp)
        total = len(recs)
        # Complete genomes first: they are the records that certainly carry the VP24 CDS.
        def rank(r: dict) -> tuple:
            d = str(r.get("description", "")).lower()
            return (0 if "complete genome" in d or "complete sequence" in d else 1,
                    str(r.get("accession", "")))

        recs = sorted(recs, key=rank)
        cap = ENA_CAP[sp]
        sample = recs[:cap]
        calls: list[dict] = []
        failures: list[str] = []
        for rec in sample:
            acc = rec["accession"]
            try:
                trans = ena_vp24_translations(acc, "g1b_ena_%s_%s.embl" % (sp, acc))
            except Exception as exc:  # noqa: BLE001
                failures.append("%s (%s)" % (acc, exc))
                continue
            for gname, pname, seq in trans:
                pos = align_positions(seq)
                calls.append({
                    "accession": acc,
                    "description": rec.get("description", ""),
                    "gene": gname,
                    "product": pname,
                    "length": len(seq),
                    "positions": pos,
                    "bdbv_type_at_all_four": is_bdbv_type(pos),
                })
        bdbv_like = [c for c in calls if c["bdbv_type_at_all_four"]]
        distinct = {c["length"]: None for c in calls}
        result[sp] = {
            "n_ena_records_total": total,
            "n_sampled": len(sample),
            "n_vp24_cds_called": len(calls),
            "n_with_all_four_bdbv_type": len(bdbv_like),
            "complete_genomes_first": True,
            "embm_fetch_failures": failures,
            "calls": calls,
        }
        say("  %-6s records=%5d  sampled=%4d  VP24 CDS called=%4d  carrying %s: %d"
            % (sp, total, len(sample), len(calls), BDBV_TYPE, len(bdbv_like)))
        say("        EMBL fetch failures: %d" % len(failures))
        for key, grp in sorted(_group_calls(calls).items(), key=lambda kv: -kv[1]["n"]):
            say("        %s x%-4d len=%s  %s"
                % (key, grp["n"], grp["lengths"], grp["example"]))
    return result


def _group_calls(calls: list[dict]) -> dict:
    groups: dict[str, dict] = {}
    for c in calls:
        key = calls_key(c["positions"])
        g = groups.setdefault(key, {"n": 0, "lengths": sorted({c["length"]}), "example": ""})
        g["n"] += 1
        if not g["example"]:
            g["example"] = "%s %s" % (c["accession"], c["product"][:26])
    return groups


# ----------------------------------------------------------------------- BV-BRC
# BV-BRC keeps VP24 as a ``genome_feature`` row: it carries ``aa_sequence_md5`` and
# ``aa_length`` but no sequence text, and on 2026-09-26 the ``protein`` /
# ``protein_sequence`` data types answered HTTP 404 (BV-BRC API v1.9.3).  A BV-BRC
# ``protein_id`` is a GenBank/RefSeq protein accession, so the sequence text is
# recovered from NCBI E-utilities and the md5 is re-checked before a call is used.
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _add_fasta(out: dict[str, str], chunk: list[str]) -> None:
    """Parse one FASTA record held as [header, seqlines...] into `out`."""
    acc = chunk[0].split()[0]
    if "|" in acc:
        parts = [p for p in acc.split("|") if p]
        if len(parts) >= 2:
            acc = parts[1]
    seq = "".join(chunk[1:]).replace("*", "").replace(" ", "")
    if acc and seq:
        out[acc] = seq


def ncbi_protein_sequences(accessions: list[str], tag_base: str) -> dict[str, str]:
    """{'protein_accession': aa_sequence} from NCBI efetch(db=protein)."""
    out: dict[str, str] = {}
    batch = 40
    for i in range(0, len(accessions), batch):
        params = urllib.parse.urlencode({
            "db": "protein",
            "id": ",".join(accessions[i:i + batch]),
            "rettype": "fasta",
            "retmode": "text",
        })
        body = fetch("%s?%s" % (NCBI_EFETCH, params),
                     "%s_%03d.fasta" % (tag_base, i // batch), sleep_after=0.4)
        chunk: list[str] = []
        for line in body.decode("utf-8", "replace").splitlines():
            if line.startswith(">"):
                if chunk:
                    _add_fasta(out, chunk)
                chunk = [line[1:].strip()]
            elif chunk:
                chunk.append(line.strip())
        if chunk:
            _add_fasta(out, chunk)
    return out


def route_bvbrc() -> dict:
    say("=" * 78)
    say("BR  BV-BRC -- every distinct VP24 protein sequence, screened position by position")
    result: dict[str, dict] = {}
    for sp, taxid in SPECIES.items():
        path = os.path.join(G1_RAW, "bvbrc_%s_features.json" % sp)
        if not os.path.exists(path):
            say("  %-6s SKIPPED (no cached feature table)" % sp)
            result[sp] = {"error": "no cached feature table"}
            continue
        with open(path, encoding="utf-8") as fh:
            rows = json.load(fh)
        vp24 = [r for r in rows
                if str(r.get("gene", "")).upper().endswith("VP24")
                or "VP24" in str(r.get("product", "")).upper()]
        by_md5: dict[str, list[dict]] = defaultdict(list)
        for r in vp24:
            if r.get("aa_sequence_md5"):
                by_md5[r["aa_sequence_md5"]].append(r)
        representative = {m: max(rs, key=lambda r: str(r.get("date_modified", "")))
                          for m, rs in by_md5.items()}
        pids = sorted({r["protein_id"] for r in representative.values()
                       if r.get("protein_id")})
        seqs = ncbi_protein_sequences(pids, "g1b_bvbrc_%s_protein" % sp)
        calls: list[dict] = []
        missing: list[str] = []
        md5_mismatch: list[str] = []
        for m, rec in representative.items():
            pid = rec.get("protein_id")
            seq = seqs.get(pid, "")
            if not seq:
                missing.append(pid or m)
                continue
            pos = align_positions(seq)
            consistent = md5_of(seq) == m
            if not consistent:
                md5_mismatch.append("%s (BV-BRC md5 %s, NCBI md5 %s)"
                                    % (pid, m, md5_of(seq)))
            calls.append({
                "aa_sequence_md5": m,
                "protein_id": pid,
                "accession": rec.get("accession"),
                "length": len(seq),
                "md5_consistent": consistent,
                "positions": pos,
                "bdbv_type_at_all_four": is_bdbv_type(pos),
            })
        bdbv_like = [c for c in calls if c["bdbv_type_at_all_four"]]
        # residue distribution at each of the four positions across distinct proteins
        dist = {str(p): dict(sorted(Counter(
            c["positions"].get(p) or "-" for c in calls).items(), key=lambda kv: -kv[1]))
            for p in PRIMARY}
        result[sp] = {
            "n_vp24_features": len(vp24),
            "n_distinct_aa_md5": len(by_md5),
            "n_sequences_screened": len(calls),
            "n_sequences_unavailable": len(missing),
            "n_with_all_four_bdbv_type": len(bdbv_like),
            "sequence_source": "NCBI efetch(db=protein) by BV-BRC protein_id",
            "md5_mismatches": md5_mismatch,
            "residue_distribution": dist,
            "calls": calls,
        }
        say("  %-6s VP24 features=%5d distinct proteins=%3d screened=%3d (unavailable %d)"
            % (sp, len(vp24), len(by_md5), len(calls), len(missing)))
        say("        md5 re-verified against BV-BRC: %d/%d identical, %d mismatch"
            % (len(calls) - len(md5_mismatch), len(calls), len(md5_mismatch)))
        say("        carrying %s at all four: %d" % (BDBV_TYPE, len(bdbv_like)))
        for p in PRIMARY:
            say("        pos %3d: %s" % (p, result[sp]["residue_distribution"][str(p)]))
    return result


# ------------------------------------------------------------------------ UniProt
def route_uniprot_local() -> dict:
    """Re-read the UniProt FASTAs already archived by G1 (no new network calls)."""
    say("=" * 78)
    say("UP  UniProt -- re-scored from the G1 payloads (no new downloads)")
    result: dict[str, dict] = {}
    for sp in SPECIES:
        path = os.path.join(G1_RAW, "uniprot_%s_vp24.fasta" % sp)
        if not os.path.exists(path):
            result[sp] = {"error": "no cached fasta"}
            continue
        calls = []
        for rec in SeqIO.parse(path, "fasta"):
            seq = str(rec.seq)
            pos = align_positions(seq)
            calls.append({"id": rec.id, "length": len(seq), "positions": pos,
                          "bdbv_type_at_all_four": is_bdbv_type(pos)})
        dist = {str(p): dict(sorted(Counter(c["positions"].get(p) or "-" for c in calls).items(),
                                    key=lambda kv: -kv[1])) for p in PRIMARY}
        result[sp] = {"n_entries": len(calls),
                      "n_with_all_four_bdbv_type": sum(c["bdbv_type_at_all_four"] for c in calls),
                      "residue_distribution": dist, "calls": calls}
        say("  %-6s entries=%3d  carrying %s: %d   pos83=%s pos135=%s pos140=%s pos141=%s"
            % (sp, len(calls), BDBV_TYPE, result[sp]["n_with_all_four_bdbv_type"],
               dist["83"], dist["135"], dist["140"], dist["141"]))
    return result


def probe() -> None:
    """Cheap reachability check for the two repaired routes."""
    say("PROBE -- ENA field list and NCBI protein route for BV-BRC protein_ids")
    recs = ena_list(565995, "g1b_probe_ena_bdbv.tsv", limit=20)
    say("  ENA tax_eq(565995) rows=%d  first=%s" % (len(recs), recs[0] if recs else None))
    # TAFV VP24 YP_003815430.1 is a BV-BRC protein_id whose md5 is on record.
    expect = {"YP_003815430.1": "06d2b6d630dc9c0b45a387a4f4cab31f"}
    seqs = ncbi_protein_sequences(sorted(expect), "g1b_probe_ncbi")
    for pid, want in expect.items():
        s = seqs.get(pid, "")
        say("  NCBI protein %s len=%d md5=%s  BV-BRC md5=%s  match=%s"
            % (pid, len(s), md5_of(s)[:16] if s else "-", want[:16],
               md5_of(s) == want if s else False))


def main() -> None:
    global REUSE_RAW
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true", help="only check route reachability")
    ap.add_argument("--skip-ena", action="store_true")
    ap.add_argument("--reuse-raw", action="store_true",
                    help="re-score the archived payloads instead of re-downloading them")
    args = ap.parse_args()
    REUSE_RAW = args.reuse_raw

    os.makedirs(OUT, exist_ok=True)
    if args.probe:
        probe()
        with open(os.path.join(OUT, "g1b_probe_log.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(LOG) + "\n")
        return

    started = time.strftime("%Y-%m-%dT%H:%M:%S")
    say("G1b cross-database screen -- started %s" % started)
    say("reference: EBOV VP24 P0C779 (%d aa), positions %s, BDBV-type %s"
        % (len(EBOV_VP24), PRIMARY, BDBV_TYPE))
    summary: dict = {
        "started": started,
        "reference": {"id": "P0C779", "length": len(EBOV_VP24), "md5": md5_of(EBOV_VP24)},
        "positions": list(PRIMARY),
        "bdbv_type": BDBV_TYPE,
        "upstream": "tools/g1_crossdb_census_20260926.py (2026-09-26T01:18:54-01:21:07)",
        "ena_cap_per_species": ENA_CAP,
        "reused_archived_payloads": bool(args.reuse_raw),
        "uniprot": route_uniprot_local(),
        "bvbrc": route_bvbrc(),
    }
    if not args.skip_ena:
        summary["ena"] = route_ena()
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")

    with open(os.path.join(OUT, "g1b_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(OUT, "g1b_log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    say("")
    say("wrote %s" % os.path.join(OUT, "g1b_summary.json"))


if __name__ == "__main__":
    main()
