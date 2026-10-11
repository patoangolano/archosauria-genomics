"""Step 10 - resolve bibliography DOIs against CrossRef and verify they exist.

Standalone (no platform helpers): queries api.crossref.org directly.
Writes data/bibliografia.json.
"""
import json, time, urllib.parse, urllib.request, os, sys

API = "https://api.crossref.org/works"
MAILTO = os.environ.get("CROSSREF_CONTACT_EMAIL")   # CrossRef's polite pool


def _get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "archosauria-genomics/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def lookup(ref_string):
    """Top CrossRef match for a free-text citation. Always inspect the returned
    title: a confident-looking hit can still be the wrong paper."""
    p = {"query.bibliographic": ref_string, "rows": 1}
    if MAILTO:
        p["mailto"] = MAILTO
    items = _get(f"{API}?{urllib.parse.urlencode(p)}")["message"]["items"]
    if not items:
        return None
    it = items[0]
    return {"doi": it.get("DOI"), "title": (it.get("title") or [""])[0],
            "journal": (it.get("container-title") or [""])[0],
            "year": (it.get("issued", {}).get("date-parts", [[None]])[0] or [None])[0],
            "authors": "; ".join(f"{a.get('given','')} {a.get('family','')}".strip()
                                 for a in (it.get("author") or [])[:6]),
            "score": it.get("score")}


def verify(doi):
    try:
        it = _get(f"{API}/{urllib.parse.quote(doi)}")["message"]
        return {"ok": True, "title": (it.get("title") or [""])[0]}
    except Exception as e:
        return {"ok": False, "error": str(e)}


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "references.tsv"
    out = {}
    for line in open(src):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, query = line.split("\t", 1)
        rec = lookup(query)
        if rec:
            rec.update(verify(rec["doi"]))
        out[key] = rec
        print(f"{key:16s} {rec['doi'] if rec else 'NOT FOUND'}")
        time.sleep(0.2)
    json.dump(out, open("../data/bibliografia.json", "w"), indent=1)
