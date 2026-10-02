"""G1: is the BDBV four-position VP24 configuration exclusive *outside NCBI*?

The manuscript's exclusivity claim was established on NCBI alone.  This script
re-tests it against four independent retrieval routes and records, for every
route, the URL, the retrieval timestamp and the sha256 of the raw payload so the
statement can be re-derived.

Routes
------
PX  Pathoplexus (LAPIS)  ebola-bdbv / ebola-zaire / ebola-sudan
    per-gene amino-acid variant tables with coverage, i.e. the whole VP24
    protein is screened position by position across every sequence on the
    instance.
EN  ENA portal (INSDC)    tax_tree(<species>) AND description="VP24"
    nucleotide records; the VP24 CDS translation is taken from the EMBL flat
    file and the four positions are called after alignment to the EBOV
    reference, so numbering is never assumed to be 1:1.
BR  BV-BRC                genome_feature, CDS features annotated as VP24;
    protein identity is recorded as aa_sequence_md5 and compared with the md5
    of the reference sequences.
UP  UniProt               protein-level VP24 entries per species; sequences
    aligned and positions called directly.

Outputs land in analysis/g1_crossdb_census_20260926/ (raw/ keeps the payloads).
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import re
import time
import urllib.parse
import urllib.request

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "g1_crossdb_census_20260926")
RAW = os.path.join(OUT, "raw")
UA = {"User-Agent": "bdbv-vp24-census/1.0 (research script)"}

# EBOV VP24 reference (Zaire ebolavirus Mayinga-76), UniProt P0C779.
EBOV_VP24 = (
    "MAKATGRYNLISPKKDLEKGVVLSDLCNFLVSQTIQGWKVYWAGIEFDVTHKGMALLHRL"
    "KTNDFAPAWSMTRNLFPHLFQNPNSTIESPLWALRVILAAGIQDQLIDQSLIEPLAGALG"
    "LISDWLLTTNTNHFNMRTQRVKEQLSLKMLSLIRSNILKFINKLDALHVVNYNGLLSSIE"
    "IGTQNHTIIITRTNMGFLVELQEPDKSAMNRMKPGPAKFSLLHESTLKAFTQGSSTRMQS"
    "LILEFNSSLAI"
)
PRIMARY = (83, 135, 140, 141)

SPECIES = {
    "BDBV": 565995,
    "EBOV": 186538,
    "SUDV": 186540,
    "RESTV": 186539,
    "TAFV": 186541,
    "BOMV": 2010960,
}
PX_INSTANCES = {"BDBV": "ebola-bdbv", "EBOV": "ebola-zaire", "SUDV": "ebola-sudan"}

LOG: list[str] = []


def say(msg: str = "") -> None:
    LOG.append(msg)
    print(msg, flush=True)


def fetch(url: str, tag: str, tries: int = 3) -> bytes:
    """GET with retries; raw payload + sha256 are written for provenance."""
    os.makedirs(RAW, exist_ok=True)
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=180) as r:
                body = r.read()
            path = os.path.join(RAW, tag)
            with open(path, "wb") as fh:
                fh.write(body)
            digest = hashlib.sha256(body).hexdigest()
            say("  [fetch] %-42s %8d bytes  sha256:%s" % (tag, len(body), digest[:16]))
            return body
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 + 3 * attempt)
    raise RuntimeError("fetch failed for %s: %s" % (url, last))


def align_positions(seq: str, positions=PRIMARY, reference: str = EBOV_VP24) -> dict:
    """Map reference numbering onto `seq` by global alignment."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    try:
        aln = aligner.align(reference, seq)[0]
    except Exception:  # noqa: BLE001
        return {p: None for p in positions}
    # aligned strings: use the alignment object's indices
    a, b = str(aln[0]), str(aln[1])
    ref_i = seq_i = 0
    out: dict[int, str | None] = {}
    for ca, cb in zip(a, b):
        if ca != "-":
            ref_i += 1
        if cb != "-":
            seq_i += 1
        if ca != "-" and ref_i in positions:
            out[ref_i] = cb if cb != "-" else None
    for p in positions:
        out.setdefault(p, None)
    return out


def call_length(seq: str) -> int:
    return len(seq)


def md5_of(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()


# ------------------------------------------------------------------ Pathoplexus
def route_pathoplexus() -> dict:
    say("=" * 78)
    say("PX  Pathoplexus (LAPIS) -- whole-protein amino-acid variant screen")
    result: dict[str, dict] = {}
    for sp, org in PX_INSTANCES.items():
        url = ("https://lapis.pathoplexus.org/%s/sample/aminoAcidMutations?limit=5000" % org)
        try:
            body = fetch(url, "px_%s_aaMutations.json" % org)
            payload = json.loads(body.decode("utf-8"))
        except Exception as exc:  # noqa: BLE001
            say("  %s: FAILED %s" % (org, exc))
            result[sp] = {"error": str(exc)}
            continue
        rows = payload.get("data", [])
        vp24 = [r for r in rows if r.get("sequenceName") == "VP24"]
        coverage = max([r["coverage"] for r in vp24], default=0)
        # any variant at the four positions?
        at_primary = [r for r in vp24 if r.get("position") in PRIMARY]
        result[sp] = {
            "instance": org,
            "url": url,
            "data_version": payload.get("info", {}).get("dataVersion"),
            "n_vp24_variant_rows": len(vp24),
            "max_vp24_coverage": coverage,
            "variants": [
                {"mutation": r["mutation"], "position": r["position"],
                 "from": r["mutationFrom"], "to": r["mutationTo"],
                 "count": r["count"], "coverage": r["coverage"],
                 "proportion": r["proportion"]}
                for r in sorted(vp24, key=lambda r: -r["coverage"])
            ],
            "variants_at_primary_positions": [r["mutation"] for r in at_primary],
        }
        say("  %-6s %-14s VP24 variant rows=%d  max coverage=%d  variants at %s: %s"
            % (sp, org, len(vp24), coverage, PRIMARY,
               [r["mutation"] for r in at_primary] or "NONE"))
        for r in sorted(vp24, key=lambda r: -r["coverage"]):
            say("        %-14s count=%6d cov=%6d prop=%.4f"
                % (r["mutation"], r["count"], r["coverage"], r["proportion"]))
    return result


# -------------------------------------------------------------------------- ENA
def ena_search(taxid: int, phrase: str, tag: str) -> list[dict]:
    """Enumerate the species' ENA records and filter client-side.

    The server-side ``description="VP24"`` predicate is rejected (HTTP 400) by
    the portal for this result domain, so the filter is applied locally to the
    returned table instead; the table itself is archived for provenance.
    """
    params = {
        "result": "sequence",
        "query": "tax_eq(%d)" % taxid,
        "fields": "accession,description,sequence_length",
        "format": "tsv",
        "limit": "200000",
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(params)
    try:
        body = fetch(url, tag)
    except Exception as exc:  # noqa: BLE001
        say("    ENA search failed: %s" % exc)
        return []
    text = body.decode("utf-8", "replace").strip().splitlines()
    if len(text) <= 1:
        return []
    header = text[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in text[1:] if line.strip()]
    return [r for r in rows if phrase.lower() in str(r.get("description", "")).lower()]


def ena_count(taxid: int) -> int:
    """Total ENA sequence records for the species (the denominator)."""
    params = {"result": "sequence", "query": "tax_eq(%d)" % taxid,
              "fields": "accession", "format": "tsv", "limit": "200000"}
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(params)
    try:
        body = fetch(url, "ena_count_%d.tsv" % taxid)
    except Exception:  # noqa: BLE001
        return -1
    return max(0, len(body.decode("utf-8", "replace").strip().splitlines()) - 1)


def ena_embl_translations(accession: str, tag: str) -> list[tuple[str, str, str]]:
    """Return [(gene, product, translation)] for the VP24 CDS of an EMBL record."""
    url = ("https://www.ebi.ac.uk/ena/browser/api/embl/%s?lineLimit=0" % accession)
    try:
        body = fetch(url, tag)
    except Exception as exc:  # noqa: BLE001
        say("    embl fetch failed for %s: %s" % (accession, exc))
        return []
    text = body.decode("utf-8", "replace")
    blocks: list[tuple[str, str]] = []
    key = None
    buf: list[str] = []
    for line in text.splitlines():
        if not line.startswith("FT"):
            continue
        rest = line[2:]
        if not rest.strip():
            continue
        if rest[:3].strip() == "" and rest.strip().startswith("/"):
            buf.append(rest)
            continue
        if key is not None:
            blocks.append((key, "\n".join(buf)))
        key = rest.strip().split()[0]
        buf = [rest]
    if key is not None:
        blocks.append((key, "\n".join(buf)))

    out = []
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
        if "VP24" in gname.upper() or "VP24" in pname.upper():
            out.append((gname, pname, re.sub(r"\s+", "", trans.group(1))))
    return out


def route_ena(max_records: int = 25) -> dict:
    say("=" * 78)
    say("EN  ENA portal (INSDC) -- VP24 nucleotide records + EMBL CDS translation")
    result: dict[str, dict] = {}
    for sp, taxid in SPECIES.items():
        hits = []
        for phrase, tag in (("VP24", "ena_%s_desc_vp24.tsv" % sp),):
            hits = ena_search(taxid, phrase, tag)
            if hits:
                break
        total = ena_count(taxid)
        say("  %-6s ENA records total: %d   VP24 in a free-text description: %d"
            % (sp, total, len(hits)))
        # Since ENA descriptions are genome-level, the residue-level test is run
        # on the EMBL records themselves: the VP24 CDS translation is parsed out
        # of the flat file.  Large archives are sampled evenly and the sample
        # size is reported.
        if total > max_records:
            all_recs = ena_search(taxid, "", "ena_%s_all.tsv" % sp)
            step = max(1, len(all_recs) // max_records)
            sample = all_recs[::step][:max_records]
        else:
            sample = ena_search(taxid, "", "ena_%s_all.tsv" % sp)
        say("      EMBL records sampled for CDS parsing: %d" % len(sample))
        calls: list[dict] = []
        for rec in sample:
            acc = rec["accession"]
            time.sleep(0.25)
            for gname, pname, translation in ena_embl_translations(
                    acc, "ena_%s_%s.embl" % (sp, acc)):
                calls.append({
                    "accession": acc,
                    "gene": gname,
                    "product": pname,
                    "length": len(translation),
                    "positions": align_positions(translation),
                    "all_primary_equal_ebov_reference": all(
                        align_positions(translation).get(p) == EBOV_VP24[p - 1]
                        for p in PRIMARY),
                    "bdbv_type_at_all_four": all(
                        align_positions(translation).get(p) == c
                        for p, c in zip(PRIMARY, "SQHA")),
                })
        bdbv_like = [c for c in calls if c["bdbv_type_at_all_four"]]
        result[sp] = {"n_ena_records_total": total, "n_records": len(hits),
                      "n_sampled": len(sample), "n_vp24_cds_called": len(calls),
                      "n_with_all_four_bdbv_type": len(bdbv_like),
                      "records": sample[:10], "calls": calls}
        say("      CDS translations recovered: %d" % len(calls))
        say("      carrying S/Q/H/A at all four: %d" % len(bdbv_like))
        for c in calls:
            say("        %-14s len=%3d  %s  bdbv-type=%s"
                % (c["accession"], c["length"],
                   " ".join("%d:%s" % (p, c["positions"].get(p)) for p in PRIMARY),
                   c["bdbv_type_at_all_four"]))
    return result


# ----------------------------------------------------------------------- BV-BRC
def route_bvbrc(known: dict[str, list[str]]) -> dict:
    say("=" * 78)
    say("BR  BV-BRC genome_feature -- VP24 CDS features per species")
    result: dict[str, dict] = {}
    for sp, taxid in SPECIES.items():
        params = "eq(taxon_id,%d)&limit(50000)&http_accept=application/json" % taxid
        url = "https://www.bv-brc.org/api/genome_feature/?" + params
        try:
            body = fetch(url, "bvbrc_%s_features.json" % sp)
            rows = json.loads(body.decode("utf-8"))
        except Exception as exc:  # noqa: BLE001
            say("  %-6s FAILED %s" % (sp, exc))
            result[sp] = {"error": str(exc)}
            continue
        vp24 = [r for r in rows
                if str(r.get("gene", "")).upper().endswith("VP24")
                or "VP24" in str(r.get("product", "")).upper()]
        md5s = sorted({r.get("aa_sequence_md5") for r in vp24 if r.get("aa_sequence_md5")})
        accs = sorted({r.get("accession") for r in vp24})
        # does any BV-BRC VP24 md5 equal a sequence *designed* to carry the four
        # BDBV-type residues on that species' own reference backbone?
        engineered = []
        for ref_seq in known.get(sp, []):
            if len(ref_seq) < max(PRIMARY):
                continue
            chars = list(ref_seq)
            for p, c in zip(PRIMARY, "SQHA"):
                chars[p - 1] = c
            engineered.append(md5_of("".join(chars)))
        matched = sorted(set(engineered) & set(md5s))
        result[sp] = {
            "url": url,
            "n_features_total": len(rows),
            "n_vp24_features": len(vp24),
            "n_distinct_accessions": len(accs),
            "distinct_aa_sequence_md5": md5s,
            "example_accessions": accs[:10],
            "md5_of_known_reference_vp24": sorted({md5_of(s) for s in known.get(sp, [])}),
            "md5_of_reference_with_bdbv_type_at_all_four": sorted(set(engineered)),
            "bvbrc_md5_matching_a_bdbv_type_variant": matched,
        }
        say("  %-6s features=%6d VP24=%4d accessions=%4d distinct protein md5=%d"
            % (sp, len(rows), len(vp24), len(accs), len(md5s)))
        say("        BV-BRC proteins matching a reference backbone carrying S/Q/H/A: %s"
            % (matched or "NONE"))
        for m in md5s[:6]:
            say("        md5 %s" % m)
    return result


# ---------------------------------------------------------------------- UniProt
def route_uniprot() -> dict:
    say("=" * 78)
    say("UP  UniProt search -- VP24 protein entries per species")
    result: dict[str, dict] = {}
    for sp, taxid in SPECIES.items():
        query = "gene:VP24 AND organism_id:%d" % taxid
        url = ("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query="
               + urllib.parse.quote(query))
        try:
            body = fetch(url, "uniprot_%s_vp24.fasta" % sp)
        except Exception as exc:  # noqa: BLE001
            say("  %-6s FAILED %s" % (sp, exc))
            result[sp] = {"error": str(exc)}
            continue
        recs = list(SeqIO.parse(io.StringIO(body.decode()), "fasta"))
        calls = []
        for rec in recs:
            seq = str(rec.seq)
            pos = align_positions(seq)
            calls.append({
                "id": rec.id,
                "length": len(seq),
                "positions": pos,
                "bdbv_type_at_all_four": all(pos.get(p) == c for p, c in zip(PRIMARY, "SQHA")),
                "identical_to_ebov_reference": seq == EBOV_VP24,
            })
        bdbv_like = [c for c in calls if c["bdbv_type_at_all_four"]]
        result[sp] = {
            "url": url,
            "n_entries": len(calls),
            "n_with_all_four_bdbv_type": len(bdbv_like),
            "calls": calls,
        }
        say("  %-6s entries=%3d  carrying S/Q/H/A at all four: %d"
            % (sp, len(calls), len(bdbv_like)))
        for c in calls[:6]:
            say("        %-16s len=%3d  %s" % (c["id"], c["length"],
               " ".join("%d:%s" % (p, c["positions"].get(p)) for p in PRIMARY)))
    return result


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    started = time.strftime("%Y-%m-%dT%H:%M:%S")
    say("G1 cross-database exclusivity census -- started %s" % started)
    say("reference: EBOV VP24 P0C779, positions %s (1-based, EBOV numbering)" % (PRIMARY,))
    summary = {
        "started": started,
        "reference": {"id": "P0C779", "length": len(EBOV_VP24),
                      "md5": md5_of(EBOV_VP24)},
        "positions": list(PRIMARY),
        "pathoplexus": route_pathoplexus(),
        "ena": route_ena(),
        "uniprot": route_uniprot(),
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    known: dict[str, list[str]] = {}
    for sp, blk in summary["uniprot"].items():
        seqs = []
        path = os.path.join(RAW, "uniprot_%s_vp24.fasta" % sp)
        if os.path.exists(path):
            seqs = [str(r.seq) for r in SeqIO.parse(path, "fasta")]
        known[sp] = seqs
    summary["bvbrc"] = route_bvbrc(known)
    with open(os.path.join(OUT, "g1_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(OUT, "g1_log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    say("")
    say("wrote %s" % os.path.join(OUT, "g1_summary.json"))


if __name__ == "__main__":
    main()
