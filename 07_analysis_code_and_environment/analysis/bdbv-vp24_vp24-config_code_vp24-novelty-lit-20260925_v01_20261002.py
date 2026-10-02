"""Novelty audit for two BDBV VP24 claims (Claim A and Claim B).

Claim A  BDBV natively carries the complete four-position configuration
         P83S/N135Q/R140H/V141A, conserved across outbreaks.
Claim B  No other ebolavirus carries this exact combination, and no paper has
         stated the composite claim.

Searches Europe PMC (which indexes full text), PubMed and OpenAlex; saves every
raw response; then reports, per query, the hit count and the top records so the
novelty judgement can be re-checked without trusting this summary.

Outputs land in analysis/vp24_novelty_lit_20260925/.
"""

import json
import os
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "analysis", "vp24_novelty_lit_20260925")
RAW = os.path.join(OUT, "raw")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
OPENALEX = "https://api.openalex.org/works"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "Mozilla/5.0 (research novelty audit; contact: local)"}


def get_json(url, tries=3, params=None):
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    last = None
    for _ in range(tries):
        try:
            request = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(request, timeout=120) as resp:
                return json.loads(resp.read().decode("utf-8", "replace"))
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(4)
    raise RuntimeError(f"{url} failed: {last}")


def get_text(url, tries=3):
    last = None
    for _ in range(tries):
        try:
            request = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(request, timeout=120) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(4)
    raise RuntimeError(f"{url} failed: {last}")


def save(name, payload):
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, name)
    with open(path, "w", encoding="utf-8") as fh:
        if isinstance(payload, str):
            fh.write(payload)
        else:
            json.dump(payload, fh, ensure_ascii=False, indent=1)
    return path


# Europe PMC: default field is full text, so these catch statements anywhere.
EPMC_QUERIES = {
    "epmc_vp24_bundibugyo_title_abs": 'TITLE_ABS:"Bundibugyo" AND TITLE_ABS:"VP24"',
    "epmc_bundibugyo_vp24_anywhere": '"Bundibugyo" AND "VP24"',
    "epmc_p83s_exact": '"P83S"',
    "epmc_p83s_vp24": '"P83S" AND "VP24"',
    "epmc_4x_mutant": '("P83S/N135Q/R140H/V141A")',
    "epmc_n135q_r140h": '("N135Q" AND "R140H")',
    "epmc_v141a_bundibugyo": '"V141A" AND "Bundibugyo"',
    "epmc_vp24_kpna_species": 'TITLE_ABS:"VP24" AND (TITLE_ABS:"karyopherin" OR '
                              'TITLE_ABS:"importin" OR TITLE_ABS:"KPNA")',
    "epmc_vp24_species_specific": '"VP24" AND "species-specific" AND "ebolavirus"',
    "epmc_vp24_signature": '"VP24 signature" OR "VP24 motif"',
    "epmc_bundibugyo_signature": 'TITLE_ABS:"Bundibugyo" AND (TITLE_ABS:"signature" '
                                 'OR TITLE_ABS:"motif" OR TITLE_ABS:"residue")',
    "epmc_bundibugyo_conserved_vp24": '"Bundibugyo" AND "VP24" AND ("conserved" OR '
                                      '"invariant")',
    "epmc_bundibugyo_interface_residues": '"Bundibugyo" AND ("interface" AND "residues") '
                                          'AND "VP24"',
    "epmc_bdbv_natural_config": '"bVP24" AND ("naturally" OR "native" OR "endogenous")',
}

PUBMED_QUERIES = {
    "pubmed_bundibugyo_vp24": "Bundibugyo[All Fields] AND VP24[All Fields]",
    "pubmed_vp24_kpna": "VP24[All Fields] AND (karyopherin[All Fields] OR "
                        "importin[All Fields] OR KPNA[All Fields])",
    "pubmed_bundibugyo_residue": "Bundibugyo[All Fields] AND VP24[All Fields] AND "
                                 "(residue[All Fields] OR sequence[All Fields])",
}

OPENALEX_QUERIES = {
    "openalex_bundibugyo_vp24": "Bundibugyo VP24",
    "openalex_vp24_karyopherin": "VP24 karyopherin alpha interferon antagonism ebolavirus",
    "openalex_bundibugyo_signature": "Bundibugyo virus VP24 species-specific residues",
}

# Papers whose full text decides the novelty question.
KEY_PMIDS = {
    "25121748": "structure of eV24-KPNA5 complex (defining the interface)",
    "27974555": "BDBV/RESTV/EBOV VP24-KPNA affinity + VP24 stability",
    "37243162": "engineered eVP24 4x mutant with BDBV residues",
    "35069519": "filovirus VP24 promoter assays (includes BDBV)",
    "37647113": "eight filovirus VP24 IFN pathway inhibition",
    "20394957": "BDBV vs EBOV in human PBMC",
    "42742179": "2026 BDBV outbreak virus characterisation",
    "42742180": "2026 BDBV human case treated with mAbs",
    "42732232": "BDBV pathogenesis in nonhuman primates (2026 preprint)",
}


def epmc_search(name, query, page_size=25):
    data = get_json(
        f"{EPMC}/search",
        params={"query": query, "format": "json", "pageSize": page_size,
                "resultType": "core"},
    )
    save(f"{name}.json", data)
    hits = data.get("hitCount")
    results = data.get("resultList", {}).get("result", [])
    print(f"\n=== EPMC {name}: hitCount={hits}  (query: {query})")
    for item in results[:8]:
        pmid = item.get("pmid") or "-"
        pmcid = item.get("pmcid") or "-"
        year = item.get("pubYear") or "-"
        print(f"   {pmid:>9s} {pmcid:>12s} {year:>5s}  {(item.get('title') or '')[:110]}")
    return data


def pubmed_search(name, query, retmax=25):
    data = get_json(
        f"{EUTILS}/esearch.fcgi",
        params={"db": "pubmed", "term": query, "retmax": retmax, "retmode": "json"},
    )
    save(f"{name}.json", data)
    ids = data.get("esearchresult", {}).get("idlist", [])
    print(f"\n=== PUBMED {name}: count={data.get('esearchresult', {}).get('count')} ids={len(ids)}")
    if ids:
        summ = get_json(
            f"{EUTILS}/esummary.fcgi",
            params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"},
        )
        save(f"{name}_esummary.json", summ)
        for uid in summ.get("result", {}).get("uids", []):
            record = summ["result"][uid]
            print(f"   {uid:>9s} {(record.get('pubdate') or '')[:5]:>5s}  "
                  f"{(record.get('title') or '')[:110]}")
    return ids


def openalex_search(name, query, per_page=10):
    data = get_json(
        OPENALEX,
        params={"search": query, "per-page": per_page,
                "select": "id,doi,title,publication_year,primary_location,"
                          "cited_by_count,type"},
    )
    save(f"{name}.json", data)
    print(f"\n=== OPENALEX {name}: count={data.get('meta', {}).get('count')}")
    for item in data.get("results", [])[:8]:
        print(f"   {item.get('publication_year')}  cites={item.get('cited_by_count'):>5}"
              f"  {(item.get('title') or '')[:110]}")
    return data


def fetch_fulltext(pmid):
    """Try Europe PMC for an open-access full text; return (pmcid, text)."""
    search = get_json(
        f"{EPMC}/search",
        params={"query": f"EXT_ID:{pmid}", "format": "json", "resultType": "core"},
    )
    save(f"key_pmid_{pmid}_record.json", search)
    results = search.get("resultList", {}).get("result", [])
    if not results:
        return None, None
    record = results[0]
    pmcid = record.get("pmcid")
    text = None
    if pmcid and record.get("isOpenAccess") == "Y":
        try:
            text = get_text(f"{EPMC}/{pmcid}/fullTextXML")
            save(f"{pmid}_{pmcid}_fulltext.xml", text)
        except Exception as exc:  # noqa: BLE001
            print(f"   full text fetch failed for {pmcid}: {exc}")
    return pmcid, text


def main():
    os.makedirs(RAW, exist_ok=True)
    print("=" * 100)
    print("EUROPE PMC (full-text index)")
    for name, query in EPMC_QUERIES.items():
        try:
            epmc_search(name, query)
        except Exception as exc:  # noqa: BLE001
            print(f"   !! {name} failed: {exc}")
        time.sleep(0.6)

    print("\n" + "=" * 100)
    print("PUBMED")
    for name, query in PUBMED_QUERIES.items():
        try:
            pubmed_search(name, query)
        except Exception as exc:  # noqa: BLE001
            print(f"   !! {name} failed: {exc}")
        time.sleep(0.6)

    print("\n" + "=" * 100)
    print("OPENALEX")
    for name, query in OPENALEX_QUERIES.items():
        try:
            openalex_search(name, query)
        except Exception as exc:  # noqa: BLE001
            print(f"   !! {name} failed: {exc}")
        time.sleep(0.6)

    print("\n" + "=" * 100)
    print("KEY PAPERS: identifiers and open-access status")
    rows = []
    for pmid, why in KEY_PMIDS.items():
        try:
            pmcid, text = fetch_fulltext(pmid)
        except Exception as exc:  # noqa: BLE001
            print(f"   !! {pmid} failed: {exc}")
            continue
        rows.append({"pmid": pmid, "pmcid": pmcid or "", "fulltext": bool(text),
                     "why": why})
        print(f"   {pmid:>9s} pmcid={pmcid or '-':<12s} fulltext={'yes' if text else 'no'}  {why}")
        time.sleep(0.6)
    save("key_papers_availability.json", rows)


if __name__ == "__main__":
    main()
