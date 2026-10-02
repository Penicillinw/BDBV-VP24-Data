"""Third pass: the decisive checks for the composite novelty claim.

 1. Can the Schwarz 2017 full text be retrieved (it is the paper that first
    mapped the four residues onto the eVP24-KPNA5 surface)?
 2. Do the BDBV-native forms (S83, Q135, H140, A141) co-occur with
    "Bundibugyo" anywhere in the indexed literature?  This is the decisive
    test of whether the composite statement has been published.
 3. What do the neighbouring papers say (eight-filovirus VP24 panel,
    orthoebolavirus diversity survey, VP24 nuclear-transport biophysics)?

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
TEXT = os.path.join(OUT, "fulltext")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
UA = {"User-Agent": "Mozilla/5.0 (research novelty audit; contact: local)"}


def get(url, params=None, tries=3):
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
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


def save(name, text):
    os.makedirs(RAW, exist_ok=True)
    with open(os.path.join(RAW, name), "w", encoding="utf-8") as fh:
        fh.write(text)


def epmc_count(query):
    data = json.loads(get(f"{EPMC}/search",
                          {"query": query, "format": "json", "pageSize": 5,
                           "resultType": "lite"}))
    save(f"targeted_{abs(hash(query)) % 10**8}.json",
         json.dumps({"query": query, "hits": data.get("hitCount"),
                     "results": [
                         {k: r.get(k) for k in ("id", "source", "pmid", "pmcid",
                                                "pubYear", "title")}
                         for r in data.get("resultList", {}).get("result", [])
                     ]}, ensure_ascii=False, indent=1))
    return data


def main():
    os.makedirs(TEXT, exist_ok=True)

    # ---- 1. Schwarz 2017 full text through the Europe PMC endpoint
    print("=" * 90)
    print("1. Schwarz 2017 (PMC5286879) full text via Europe PMC")
    try:
        xml = get(f"{EPMC}/PMC5286879/fullTextXML")
        save("PMC5286879_epmc_fulltext.xml", xml)
        print(f"   retrieved {len(xml)} bytes")
        import re
        flat = re.sub(r"<[^>]+>", " ", xml)
        flat = re.sub(r"\s+", " ", flat).strip()
        with open(os.path.join(TEXT, "PMC5286879_epmc.txt"), "w", encoding="utf-8") as fh:
            fh.write(flat)
        for needle in ("P83", "N135", "R140", "V141", "83", "135", "140", "141",
                       "alignment", "Bundibugyo"):
            idx = flat.find(needle)
            print(f"   '{needle}': {'found' if idx >= 0 else 'absent'}")
    except Exception as exc:  # noqa: BLE001
        print(f"   FAILED: {exc}")

    # ---- 2. decisive composite-claim searches
    print("\n" + "=" * 90)
    print("2. Does the BDBV-native residue set ever co-occur with Bundibugyo?")
    composite_queries = {
        '"Bundibugyo" AND "S83"': '"Bundibugyo" AND "S83"',
        '"Bundibugyo" AND "Q135"': '"Bundibugyo" AND "Q135"',
        '"Bundibugyo" AND "H140"': '"Bundibugyo" AND "H140"',
        '"Bundibugyo" AND "A141"': '"Bundibugyo" AND "A141"',
        '"Bundibugyo" AND "SQHA"': '"Bundibugyo" AND "SQHA"',
        '"Bundibugyo" AND "four positions"': '"Bundibugyo" AND "four positions"',
        '"Bundibugyo" AND "four residues"': '"Bundibugyo" AND "four residues"',
        '"Bundibugyo" AND "signature residue"': '"Bundibugyo" AND "signature residue"',
        '"Bundibugyo" AND "natively"': '"Bundibugyo" AND "natively"',
        '"Bundibugyo" AND "exclusive" AND "VP24"': '"Bundibugyo" AND "exclusive" AND "VP24"',
        '"VP24" AND "BDBV-specific"': '"VP24" AND "BDBV-specific"',
        '"VP24" AND "species-specific signature"': '"VP24" AND "species-specific signature"',
    }
    report = {}
    for label, query in composite_queries.items():
        try:
            data = epmc_count(query)
            hits = data.get("hitCount")
            results = data.get("resultList", {}).get("result", [])
            report[label] = {
                "hits": hits,
                "top": [f"{r.get('pmid') or r.get('id')} {r.get('pubYear')} "
                        f"{(r.get('title') or '')[:90]}" for r in results],
            }
            print(f"\n   {label:<45s} hits={hits}")
            for r in results[:5]:
                print(f"      - {r.get('pmid') or r.get('id')} {r.get('pubYear')} "
                      f"{(r.get('title') or '')[:95]}")
        except Exception as exc:  # noqa: BLE001
            print(f"\n   {label:<45s} FAILED: {exc}")
        time.sleep(0.5)
    save("composite_claim_searches.json", json.dumps(report, ensure_ascii=False, indent=1))

    # ---- 3. neighbouring papers: abstracts to quote
    print("\n" + "=" * 90)
    print("3. Neighbouring abstracts")
    for identifier in ("37647113", "41141938", "40872766", "42277463"):
        try:
            data = json.loads(get(f"{EPMC}/search",
                                  {"query": f"EXT_ID:{identifier}", "format": "json",
                                   "resultType": "core"}))
            results = data.get("resultList", {}).get("result", [])
            if not results:
                print(f"   {identifier}: no record")
                continue
            record = results[0]
            abstract = (record.get("abstractText") or "")
            save(f"epmc_record_{identifier}.json",
                 json.dumps(data, ensure_ascii=False, indent=1))
            with open(os.path.join(TEXT, f"{identifier}_abstract.txt"), "w",
                      encoding="utf-8") as fh:
                fh.write(f"{record.get('title')}\n\n{abstract}\n")
            print(f"\n   --- {identifier} ({record.get('pubYear')}) {record.get('title')}")
            print(f"   {abstract[:900]}")
        except Exception as exc:  # noqa: BLE001
            print(f"   {identifier} FAILED: {exc}")
        time.sleep(0.5)


if __name__ == "__main__":
    main()
