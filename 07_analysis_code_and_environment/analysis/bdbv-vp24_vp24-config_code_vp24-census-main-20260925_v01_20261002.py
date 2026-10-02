"""Independent census of the BDBV VP24 four-position configuration.

Two claims are checked end to end from raw NCBI responses:

  Claim 1 (within-species conservation / cross-outbreak)
      Every full-length BDBV VP24 sequence in public databases carries
      S/Q/H/A at the positions corresponding to EBOV VP24 83/135/140/141.

  Claim 2 (exclusivity)
      No non-BDBV ebolavirus carries the exact four-residue combination
      S/Q/H/A at those positions.

The four positions are transferred from the EBOV VP24 reference by a real
global alignment (BLOSUM62); numbering is never assumed to be 1:1.

Three independent retrieval routes are used on the BDBV side, because the
answer depends on what is counted:

  route P   db=protein  "...[Organism] AND VP24[All Fields]"
            the record set that yields the widely quoted count
  route N   db=nuccore  "...[Organism] AND VP24[All Fields]"
            the corresponding annotated genome records
  route G   db=nuccore  "...[Organism] AND complete genome[Title]"
            includes unannotated 2026 surveillance genomes that the VP24
            text search misses; the VP24 CDS is located and translated here

Outputs land in analysis/vp24_census_main_20260925/.
"""

import csv
import json
import os
import time
import urllib.parse
import urllib.request
from io import StringIO

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "vp24_census_main_20260925")
RAW = os.path.join(OUT, "raw")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# EBOV VP24 reference protein (Zaire ebolavirus Mayinga-76), UniProt P0C779.
EBOV_VP24 = (
    "MAKATGRYNLISPKKDLEKGVVLSDLCNFLVSQTIQGWKVYWAGIEFDVTHKGMALLHRL"
    "KTNDFAPAWSMTRNLFPHLFQNPNSTIESPLWALRVILAAGIQDQLIDQSLIEPLAGALG"
    "LISDWLLTTNTNHFNMRTQRVKEQLSLKMLSLIRSNILKFINKLDALHVVNYNGLLSSIE"
    "IGTQNHTIIITRTNMGFLVELQEPDKSAMNRMKPGPAKFSLLHESTLKAFTQGSSTRMQS"
    "LILEFNSSLAI"
)

PRIMARY = (83, 135, 140, 141)   # the four engineered BDBV-type positions
SECONDARY = (184, 217)          # additional BDBV interface substitutions
TARGETS = PRIMARY + SECONDARY

BDBV = "Bundibugyo virus[Organism]"
BDBV_VP24_HITS = f"{BDBV} AND VP24[All Fields]"
BDBV_COMPLETE = f"{BDBV} AND complete genome[Title]"

# Orthoebolavirus (six species) plus two out-of-genus VP24 orthologues, so the
# exclusivity statement can be read at the right taxonomic level.
SPECIES = {
    "BDBV": BDBV,
    "EBOV": "Zaire ebolavirus[Organism]",
    "SUDV": "Sudan ebolavirus[Organism]",
    "TAFV": "Tai Forest ebolavirus[Organism]",
    "RESTV": "Reston ebolavirus[Organism]",
    "BOMV": "Bombali virus[Organism]",
    "MARV": "Marburg marburgvirus[Organism]",
    "LLOV": "Lloviu virus[Organism]",
}


# ------------------------------------------------------------------ network
def get(endpoint, params, tries=4):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    last = None
    for _ in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(4)
    raise RuntimeError(f"{endpoint} failed: {last}")


def esearch(term, db, retmax=2000):
    text = get(
        "esearch.fcgi",
        {"db": db, "term": term, "retmax": retmax, "retmode": "json"},
    )
    return json.loads(text)["esearchresult"]


def efetch(ids, db, rettype, retmode="text", chunk=100):
    parts = []
    for i in range(0, len(ids), chunk):
        parts.append(
            get(
                "efetch.fcgi",
                {
                    "db": db,
                    "id": ",".join(ids[i : i + chunk]),
                    "rettype": rettype,
                    "retmode": retmode,
                },
            )
        )
        time.sleep(0.45)
    return "\n".join(parts)


def parse_fasta(text):
    records, header, seq = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def tokens(defline):
    return defline.split()[0]


# ---------------------------------------------------------------- alignment
ALIGNER = PairwiseAligner()
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score = -11
ALIGNER.extend_gap_score = -1
ALIGNER.mode = "global"

VALID_AA = set("ARNDCQEGHILKMFPSTWYVBZX")


def clean_protein(seq):
    """Keep only residues BLOSUM62 knows; anything else becomes X, length-preserving."""
    seq = "".join(ch for ch in seq.upper() if ch.isalpha() or ch == "*")
    return "".join(ch if ch in VALID_AA else "X" for ch in seq.replace("*", ""))


def map_positions(seq, ref=EBOV_VP24):
    """Residues of `seq` sitting at each reference (EBOV) position in TARGETS."""
    seq = clean_protein(seq)
    aln = ALIGNER.align(ref, seq)[0]
    a, b = aln[0], aln[1]
    ref_i, seq_i = -1, -1
    want = {p: None for p in TARGETS}
    gap_cols = 0
    for ca, cb in zip(a, b):
        if ca != "-":
            ref_i += 1
        if cb != "-":
            seq_i += 1
        else:
            gap_cols += 1
        if ca != "-" and (ref_i + 1) in want and want[ref_i + 1] is None:
            want[ref_i + 1] = cb if cb != "-" else "?"
    return want, gap_cols


def config_string(res):
    return "".join(str(res[p]) if res[p] else "?" for p in PRIMARY)


def scan_proteins(records):
    """records: iterable of (accession, description, sequence)."""
    rows = []
    for acc, desc, seq in records:
        seq = clean_protein(seq)
        res, gap_cols = map_positions(seq)
        rows.append(
            {
                "accession": acc,
                "description": desc[:150],
                "length": len(seq),
                "full_length": len(seq) == len(EBOV_VP24),
                "alignment_gap_columns": gap_cols,
                **{f"pos{p}": res[p] for p in TARGETS},
                "four_position_config": config_string(res),
            }
        )
    return rows


# --------------------------------------------------------------- genbank io
def genbank_records(ids, db, chunk=50):
    text = efetch(ids, db, "gb", chunk=chunk)
    return list(SeqIO.parse(StringIO(text), "genbank")), text


SOURCE_KEYS = ("isolate", "strain", "collection_date", "geo_loc_name",
               "country", "host", "isolation_source", "note")


def source_meta(record):
    meta = {}
    for feature in record.features:
        if feature.type == "source":
            for key in SOURCE_KEYS:
                if key in feature.qualifiers:
                    meta[key] = feature.qualifiers[key][0]
            break
    return meta


def vp24_cds(record):
    """Find the VP24 CDS of a genome record; return (protein, product)."""
    for feature in record.features:
        if feature.type != "CDS":
            continue
        product = " ".join(feature.qualifiers.get("product", []))
        note = " ".join(feature.qualifiers.get("note", []))
        blob = f"{product} {note}".lower()
        if "vp24" in blob.replace(" ", "") or "membrane-associated" in blob \
                or "membrane associated" in blob:
            translation = feature.qualifiers.get("translation", [None])[0]
            if translation:
                return translation.replace(" ", "").rstrip("*"), (product or note).strip()
    return None, None


def write_tsv(path, rows):
    if not rows:
        return
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def cached_json(path, producer):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    value = producer()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(value, fh, indent=1)
    return value


def cached_text(path, producer):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return fh.read()
    value = producer()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(value)
    return value


def main():
    os.makedirs(RAW, exist_ok=True)
    summary = {
        "reference": {"protein": "EBOV VP24 / UniProt P0C779", "length": len(EBOV_VP24)},
        "queries": {"bdbv_protein": BDBV_VP24_HITS, "bdbv_nuccore": BDBV_VP24_HITS,
                    "bdbv_complete_genomes": BDBV_COMPLETE},
    }

    # ------------------------------------------------ route P: BDBV proteins
    srch = cached_json(
        os.path.join(RAW, "bdbv_protein_esearch.json"),
        lambda: esearch(BDBV_VP24_HITS, "protein", retmax=500),
    )
    ids = srch.get("idlist", [])
    summary["route_P"] = {"term": BDBV_VP24_HITS, "esearch_count": int(srch["count"]),
                          "ids_retrieved": len(ids)}

    gb_text = cached_text(
        os.path.join(RAW, "bdbv_protein.gb"),
        lambda: efetch(ids, "protein", "gb", chunk=16),
    )
    gb_records = list(SeqIO.parse(StringIO(gb_text), "genbank"))

    prot_rows = []
    for record in gb_records:
        seq = clean_protein(str(record.seq))
        acc = record.id
        res, gap_cols = map_positions(seq)
        meta = source_meta(record)
        coded_by = ""
        for feature in record.features:
            if feature.type == "CDS" and "coded_by" in feature.qualifiers:
                coded_by = feature.qualifiers["coded_by"][0]
                break
        prot_rows.append(
            {
                "accession": acc,
                "definition": record.description[:150],
                "length": len(seq),
                "full_length": len(seq) == len(EBOV_VP24),
                "alignment_gap_columns": gap_cols,
                **{f"pos{p}": res[p] for p in TARGETS},
                "four_position_config": config_string(res),
                "coded_by": coded_by,
                "isolate": meta.get("isolate", meta.get("strain", "")),
                "collection_date": meta.get("collection_date", ""),
                "geo_loc_name": meta.get("geo_loc_name", meta.get("country", "")),
                "host": meta.get("host", ""),
            }
        )
    write_tsv(os.path.join(OUT, "routeP_bdbv_vp24_proteins.tsv"), prot_rows)

    # convenience view: the full-length BDBV records, grouped by outbreak year
    def outbreak_year(row):
        date = row.get("collection_date") or ""
        for year in ("2007", "2012", "2026"):
            if year in date:
                return year
        return "undated"

    breakdown = [
        {
            "accession": r["accession"],
            "outbreak_year": outbreak_year(r),
            "collection_date": r["collection_date"],
            "geo_loc_name": r["geo_loc_name"],
            "isolate": r["isolate"],
            "host": r["host"],
            "length": r["length"],
            "pos83": r["pos83"],
            "pos135": r["pos135"],
            "pos140": r["pos140"],
            "pos141": r["pos141"],
            "pos184": r["pos184"],
            "pos217": r["pos217"],
            "four_position_config": r["four_position_config"],
        }
        for r in sorted(prot_rows, key=lambda x: (outbreak_year(x), x["accession"]))
        if r["full_length"]
    ]
    write_tsv(os.path.join(OUT, "bdbv_full_length_breakdown.tsv"), breakdown)

    full = [r for r in prot_rows if r["full_length"]]
    partial = [r for r in prot_rows if not r["full_length"]]
    cfg = {}
    for r in prot_rows:
        cfg[r["four_position_config"]] = cfg.get(r["four_position_config"], 0) + 1
    summary["route_P"].update(
        {
            "records": len(prot_rows),
            "full_length": len(full),
            "partial": len(partial),
            "config_counts_all": cfg,
            "config_counts_full_length_only": {},
            "partial_records": [
                {"accession": r["accession"], "length": r["length"],
                 "definition": r["definition"]} for r in partial
            ],
            "claim1_supported": bool(full) and all(
                r["four_position_config"] == "SQHA" for r in full
            ),
        }
    )
    for r in full:
        key = r["four_position_config"]
        d = summary["route_P"]["config_counts_full_length_only"]
        d[key] = d.get(key, 0) + 1

    # ------------------------------------- route N: BDBV VP24 genome records
    srch_n = cached_json(
        os.path.join(RAW, "bdbv_nuccore_esearch.json"),
        lambda: esearch(BDBV_VP24_HITS, "nuccore", retmax=500),
    )
    ids_n = srch_n.get("idlist", [])
    text_n = cached_text(
        os.path.join(RAW, "bdbv_nuccore.gb"),
        lambda: efetch(ids_n, "nuccore", "gb", chunk=16),
    )
    recs_n = list(SeqIO.parse(StringIO(text_n), "genbank"))
    nuccore_rows = []
    for record in recs_n:
        protein, product = vp24_cds(record)
        meta = source_meta(record)
        row = {
            "accession": record.id,
            "definition": record.description[:150],
            "genome_length": len(record.seq),
            "vp24_found": bool(protein),
            "vp24_length": len(protein) if protein else 0,
            "isolate": meta.get("isolate", meta.get("strain", "")),
            "collection_date": meta.get("collection_date", ""),
            "geo_loc_name": meta.get("geo_loc_name", meta.get("country", "")),
            "host": meta.get("host", ""),
        }
        if protein:
            protein = clean_protein(protein)
            res, gap_cols = map_positions(protein)
            row.update({f"pos{p}": res[p] for p in TARGETS})
            row["four_position_config"] = config_string(res)
            row["alignment_gap_columns"] = gap_cols
        nuccore_rows.append(row)
    write_tsv(os.path.join(OUT, "routeN_bdbv_nuccore_vp24.tsv"), nuccore_rows)
    summary["route_N"] = {
        "term": BDBV_VP24_HITS,
        "esearch_count": int(srch_n["count"]),
        "records": len(nuccore_rows),
        "vp24_translated": sum(1 for r in nuccore_rows if r["vp24_found"]),
        "config_counts": {},
    }
    for r in nuccore_rows:
        if r["vp24_found"]:
            key = r["four_position_config"]
            d = summary["route_N"]["config_counts"]
            d[key] = d.get(key, 0) + 1

    # ----------------------------------- route G: all BDBV complete genomes
    srch_g = cached_json(
        os.path.join(RAW, "bdbv_complete_esearch.json"),
        lambda: esearch(BDBV_COMPLETE, "nuccore", retmax=500),
    )
    ids_g = srch_g.get("idlist", [])
    text_g = cached_text(
        os.path.join(RAW, "bdbv_complete_genomes.gb"),
        lambda: efetch(ids_g, "nuccore", "gb", chunk=16),
    )
    recs_g = list(SeqIO.parse(StringIO(text_g), "genbank"))
    genome_rows = []
    for record in recs_g:
        protein, product = vp24_cds(record)
        meta = source_meta(record)
        row = {
            "accession": record.id,
            "definition": record.description[:150],
            "genome_length": len(record.seq),
            "vp24_found": bool(protein),
            "vp24_length": len(protein) if protein else 0,
            "isolate": meta.get("isolate", meta.get("strain", "")),
            "collection_date": meta.get("collection_date", ""),
            "geo_loc_name": meta.get("geo_loc_name", meta.get("country", "")),
            "host": meta.get("host", ""),
            "year_hint": "",
        }
        if protein:
            protein = clean_protein(protein)
            res, gap_cols = map_positions(protein)
            row.update({f"pos{p}": res[p] for p in TARGETS})
            row["four_position_config"] = config_string(res)
            row["alignment_gap_columns"] = gap_cols
        genome_rows.append(row)
    write_tsv(os.path.join(OUT, "routeG_bdbv_complete_genomes.tsv"), genome_rows)
    g_cfg, g_found = {}, 0
    for r in genome_rows:
        if r["vp24_found"]:
            g_found += 1
            key = r["four_position_config"]
            g_cfg[key] = g_cfg.get(key, 0) + 1
    summary["route_G"] = {
        "term": BDBV_COMPLETE,
        "esearch_count": int(srch_g["count"]),
        "records": len(genome_rows),
        "vp24_translated": g_found,
        "config_counts": g_cfg,
        "claim1_supported": all(k == "SQHA" for k in g_cfg),
    }

    # ------------------------------------------------------ pan-species scan
    species_rows, per_species = [], {}
    for code, organism in SPECIES.items():
        term = f"{organism} AND VP24[All Fields]"
        res = cached_json(
            os.path.join(RAW, f"{code}_vp24_esearch.json"),
            lambda term=term: esearch(term, "protein", retmax=1500),
        )
        ids_s = res.get("idlist", [])
        fasta = cached_text(
            os.path.join(RAW, f"{code}_vp24.fasta"),
            lambda ids_s=ids_s: efetch(ids_s, "protein", "fasta"),
        )
        records = [
            (tokens(d), d, s) for d, s in parse_fasta(fasta) if len(s) > 80
        ]
        rows = scan_proteins(records)
        for row in rows:
            row["species"] = code
        species_rows.extend(rows)
        per_species[code] = {
            "term": term,
            "esearch_count": int(res["count"]),
            "records_scanned": len(rows),
            "full_length": sum(1 for r in rows if r["full_length"]),
            "config_counts": {},
            "residues_by_position": {f"pos{p}": {} for p in TARGETS},
        }
        for row in rows:
            key = row["four_position_config"]
            d = per_species[code]["config_counts"]
            d[key] = d.get(key, 0) + 1
            for p in TARGETS:
                value = str(row[f"pos{p}"])
                e = per_species[code]["residues_by_position"][f"pos{p}"]
                e[value] = e.get(value, 0) + 1

    write_tsv(os.path.join(OUT, "pan_species_vp24_positions.tsv"), species_rows)
    non_bdbv = {}
    for row in species_rows:
        if row["species"] == "BDBV" or not row["full_length"]:
            continue
        key = row["four_position_config"]
        non_bdbv[key] = non_bdbv.get(key, 0) + 1
    summary["per_species"] = per_species
    summary["non_bdbv_full_length_config_counts"] = non_bdbv
    summary["claim2_supported"] = "SQHA" not in non_bdbv
    summary["claim2_note"] = (
        "Exclusivity is asserted only over the sequences actually scanned, and "
        "per-species depth differs by orders of magnitude; any manuscript "
        "sentence must carry those counts."
    )

    with open(os.path.join(OUT, "vp24_census_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)

    print("route P (BDBV protein, VP24 text hits)")
    print(f"  records={len(prot_rows)} full={len(full)} partial={len(partial)} "
          f"configs_all={cfg} configs_full={summary['route_P']['config_counts_full_length_only']}")
    for r in full:
        print(f"    {r['accession']:<14s} {r['four_position_config']}  "
              f"isolate={r['isolate']!r} date={r['collection_date']!r} geo={r['geo_loc_name']!r}")
    for r in partial:
        print(f"    PARTIAL {r['accession']:<12s} len={r['length']:>4d} {r['definition'][:80]}")
    print("\nroute N (BDBV nuccore VP24 hits)")
    print(f"  records={len(nuccore_rows)} vp24_translated={summary['route_N']['vp24_translated']} "
          f"configs={summary['route_N']['config_counts']}")
    print("\nroute G (all BDBV complete genomes)")
    print(f"  records={len(genome_rows)} vp24_translated={g_found} configs={g_cfg}")
    print("\npan-species")
    for code, info in per_species.items():
        print(f"  {code:6s} total_hits={info['esearch_count']:>6d} scanned={info['records_scanned']:>4d} "
              f"full={info['full_length']:>4d} configs={info['config_counts']}")
    print(f"\nnon-BDBV full-length configs: {non_bdbv}")
    print(f"claim1 (route P) = {summary['route_P']['claim1_supported']}   "
          f"claim2 = {summary['claim2_supported']}")


if __name__ == "__main__":
    main()
