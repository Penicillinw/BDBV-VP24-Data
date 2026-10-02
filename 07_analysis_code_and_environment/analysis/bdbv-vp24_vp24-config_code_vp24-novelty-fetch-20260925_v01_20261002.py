"""Second pass of the VP24 novelty audit: fetch the texts that decide it.

Targets
  * PMC5286879  Schwarz 2017 (J Virol) - the paper that first mapped the four
    BDBV-type positions onto the eVP24-KPNA5 surface
  * PPR1304528  Research Square 2026 preprint comparing all seven BDBV
    proteins across the 2007, 2012 and 2026 outbreaks (novelty threat)
  * the 2026 outbreak genomics / clinical papers and the RPE-VP24 paper

Outputs land in analysis/vp24_novelty_lit_20260925/.
"""

import json
import os
import re
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "vp24_novelty_lit_20260925")
RAW = os.path.join(OUT, "raw")
TEXT = os.path.join(OUT, "fulltext")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "Mozilla/5.0 (research novelty audit; contact: local)"}

# ids are Europe PMC ids: PMIDs for journal articles, PPR ids for preprints
EPMC_IDS = [
    "PPR1304528",   # 2026 whole-ORF comparative analysis of BDBV, 3 outbreaks
    "42246951",     # 2026 outbreak: role of genomics
    "42277463",     # 2026 index case, clinical + genomic
    "42708069",     # 2026 re-emergence review
    "41953650",     # human RPE cells vs ebolavirus VP24
    "42742179",     # Nat Med 2026, emergence of BDBV variant
    "42742180",     # Nat Med 2026, treated human case
    "42732232",     # bioRxiv 2026, BDBV NHP pathogenesis
    "27974555",     # Schwarz 2017
    "37647113",     # Khan 2023, eight filovirus VP24s
    "40872766",     # 2025 biophysical basis of VP24 nuclear transport
    "41141938",     # 2025 global distribution / genetic diversity of orthoebolaviruses
]

PMC_IDS = ["PMC5286879"]


def get(url, tries=3, params=None, raw=False):
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    last = None
    for _ in range(tries):
        try:
            request = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(request, timeout=120) as resp:
                data = resp.read()
                return data if raw else data.decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(4)
    raise RuntimeError(f"{url} failed: {last}")


def save_text(name, text):
    os.makedirs(RAW, exist_ok=True)
    with open(os.path.join(RAW, name), "w", encoding="utf-8") as fh:
        fh.write(text)


def epmc_record(identifier):
    data = json.loads(
        get(f"{EPMC}/search", params={"query": f"EXT_ID:{identifier} OR ID:{identifier}",
                                      "format": "json", "resultType": "core"})
    )
    save_text(f"epmc_record_{identifier}.json", json.dumps(data, ensure_ascii=False, indent=1))
    results = data.get("resultList", {}).get("result", [])
    return results[0] if results else None


def strip_tags(text):
    text = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", text, flags=re.S | re.I)
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def main():
    os.makedirs(TEXT, exist_ok=True)
    summary = []

    # ---- 1. PMC full texts (NCBI efetch works for author manuscripts too)
    for pmcid in PMC_IDS:
        try:
            xml = get(
                f"{EUTILS}/efetch.fcgi",
                params={"db": "pmc", "id": pmcid, "retmode": "xml"},
            )
            save_text(f"{pmcid}.xml", xml)
            flat = strip_tags(xml)
            with open(os.path.join(TEXT, f"{pmcid}.txt"), "w", encoding="utf-8") as fh:
                fh.write(flat)
            print(f"[pmc] {pmcid}: {len(xml)} bytes raw, {len(flat)} chars text")
        except Exception as exc:  # noqa: BLE001
            print(f"[pmc] {pmcid} FAILED: {exc}")
        time.sleep(0.5)

    # ---- 2. Europe PMC records (abstracts) for the deciding references
    for identifier in EPMC_IDS:
        try:
            record = epmc_record(identifier)
        except Exception as exc:  # noqa: BLE001
            print(f"[epmc] {identifier} FAILED: {exc}")
            continue
        if not record:
            print(f"[epmc] {identifier}: no record")
            continue
        abstract = record.get("abstractText") or ""
        title = record.get("title") or ""
        abstract_flat = strip_tags(abstract)
        if abstract_flat:
            with open(os.path.join(TEXT, f"{identifier}_abstract.txt"), "w",
                      encoding="utf-8") as fh:
                fh.write(f"{title}\n\n{abstract_flat}\n")
        flags = []
        for needle in ("VP24", "Bundibugyo", "83", "135", "140", "141", "P83S",
                       "conserve", "invariant", "signature", "interface"):
            if needle.lower() in abstract_flat.lower():
                flags.append(needle)
        summary.append(
            {
                "id": identifier,
                "source": record.get("source"),
                "pmid": record.get("pmid"),
                "pmcid": record.get("pmcid"),
                "doi": record.get("doi"),
                "year": record.get("pubYear"),
                "title": title,
                "open_access": record.get("isOpenAccess"),
                "in_epmc": record.get("inEPMC"),
                "abstract_chars": len(abstract_flat),
                "abstract_flags": flags,
            }
        )
        print(f"[epmc] {identifier:>10s} {record.get('pubYear')} "
              f"{'OA' if record.get('isOpenAccess') == 'Y' else '--'} "
              f"abs={len(abstract_flat):>5d} flags={flags}")
        print(f"        {title[:120]}")
        time.sleep(0.5)

    save_text("epmc_record_index.json", json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
