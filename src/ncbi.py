"""Thin client for the NCBI Datasets v2 REST API and E-utilities."""
import json, time, urllib.parse, urllib.request, zipfile, io, re, os

BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CONTACT = os.environ.get("NCBI_CONTACT_EMAIL")   # optional; NCBI asks for one

TAXA = {
    "Gallus_gallus": 9031,
    "Alligator_mississippiensis": 8496,
    "Crocodylus_porosus": 8502,
    "Gavialis_gangeticus": 94835,
    "Chrysemys_picta": 8478,          # gene records sit under C. p. bellii
    "Anolis_carolinensis": 28377,
    "Homo_sapiens": 9606,
}
ASSEMBLY_TAXA = dict(TAXA, Chrysemys_picta=8479, Taeniopygia_guttata=59729)


def get(url, params=None, timeout=120):
    if params:
        url = url + "?" + urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def get_zip(url, timeout=300):
    req = urllib.request.Request(url, headers={"Accept": "application/zip"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return zipfile.ZipFile(io.BytesIO(r.read()))


def parse_fasta(txt):
    name, buf = None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:], []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def gene_id_of(header):
    m = re.search(r"\[GeneID=(\d+)\]", header)
    return m.group(1) if m else None


def efetch_region(acc, start, stop, timeout=180):
    """Genomic subsequence of a RefSeq accession, 1-based inclusive."""
    p = {"db": "nuccore", "id": acc, "rettype": "fasta", "retmode": "text",
         "seq_start": max(1, int(start)), "seq_stop": int(stop)}
    if CONTACT:
        p["email"] = CONTACT
    req = urllib.request.Request(EUTILS + "?" + urllib.parse.urlencode(p))
    with urllib.request.urlopen(req, timeout=timeout) as r:
        txt = r.read().decode()
    return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))


CODON_TABLE = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _b1 in enumerate(_bases):
    for _j, _b2 in enumerate(_bases):
        for _k, _b3 in enumerate(_bases):
            CODON_TABLE[_b1 + _b2 + _b3] = _aas[_i * 16 + _j * 4 + _k]


def translate(seq):
    s = seq.upper().replace("U", "T")
    return "".join(CODON_TABLE.get(s[i:i + 3], "X") for i in range(0, len(s) - len(s) % 3, 3))


def polite_sleep(seconds=0.34):
    time.sleep(seconds)
