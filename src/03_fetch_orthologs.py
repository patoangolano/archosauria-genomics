"""Step 3 - resolve 1:1 orthologs across seven amniotes and download CDS/protein.

Writes data/orthologs.json, data/usable_genes.json and data/sequences/<GENE>.{pep,cds}.fa.
Protein and CDS are paired by exact translation, so the codon back-translation in
step 4 is unambiguous.
"""
import json, os, urllib.parse
from ncbi import BASE, TAXA, get, get_zip, parse_fasta, gene_id_of, translate, polite_sleep

SYMBOLS = ["SHH", "BMP4", "FGF8", "MSX2", "EDA", "EDAR", "WNT5A", "TBX5", "PITX1", "GLI3",
           "SOX9", "RUNX2", "PAX6", "NKX2-5", "TBX20", "GATA4", "HAND2", "TBX1", "ALX4", "LEF1",
           "CTNNB1", "RAG1", "RAG2", "BDNF", "CXCR4", "ZEB2", "EGR1", "SALL1", "SIX1", "EYA1",
           "FOXD3", "TWIST1", "PRRX1", "DLX5", "BARX1", "SATB2", "NOTCH1", "JAG1", "HEY2", "NRG1",
           "MEF2C", "TBX2", "TBX3", "FOXP2", "SMAD4", "PTCH1", "SMO", "GREM1", "WNT7A", "LMX1B"]

os.makedirs("../data/sequences", exist_ok=True)

# --- chicken gene ids (batched: the symbol endpoint pages at 20) ---
chick = {}
for i in range(0, len(SYMBOLS), 8):
    batch = SYMBOLS[i:i + 8]
    d = get(f"{BASE}/gene/symbol/{','.join(batch)}/taxon/9031", {"page_size": 100})
    for r in d.get("reports", []):
        g = r.get("gene", {})
        if g.get("symbol", "").upper() in [b.upper() for b in batch] and \
           g.get("type") == "PROTEIN_CODING":
            chick[g["symbol"].upper()] = g["gene_id"]
    polite_sleep()

# --- orthologs restricted to the seven target taxa ---
orth = {}
for sym, gid in chick.items():
    try:
        o = get(f"{BASE}/gene/id/{gid}/orthologs",
                {"taxon_filter": list(TAXA.values()), "page_size": 50})
        orth[sym] = {int(r["gene"]["tax_id"]): r["gene"]["gene_id"] for r in o.get("reports", [])}
    except Exception:
        orth[sym] = {}
    polite_sleep()
complete = sorted(s for s, h in orth.items() if len(h) == len(TAXA))
json.dump({"chicken_ids": chick, "orthologs": orth, "complete": complete},
          open("../data/orthologs.json", "w"))

# --- download protein + CDS for every ortholog ---
all_ids = sorted({str(g) for s in complete for g in orth[s].values()})
prot, cds = {}, {}
for i in range(0, len(all_ids), 30):
    url = (f"{BASE}/gene/id/{','.join(all_ids[i:i + 30])}/download?"
           + urllib.parse.urlencode([("include_annotation_type", "FASTA_PROTEIN"),
                                     ("include_annotation_type", "FASTA_CDS")], doseq=True))
    z = get_zip(url)
    for nm in z.namelist():
        if nm.endswith("protein.faa"):
            for h, s in parse_fasta(z.read(nm).decode()):
                prot.setdefault(gene_id_of(h), []).append((len(s), s))
        elif nm.endswith("cds.fna"):
            for h, s in parse_fasta(z.read(nm).decode()):
                cds.setdefault(gene_id_of(h), []).append((len(s), s))
    polite_sleep(0.5)

# --- pair the longest protein with the CDS that translates to it exactly ---
selected = {}
for g, plist in prot.items():
    plist.sort(reverse=True)
    clist = cds.get(g, [])
    chosen = None
    for _, ps in plist:
        for _, cs in clist:
            if len(cs) % 3:
                continue
            if translate(cs).rstrip("*") == ps:
                chosen = (ps, cs)
                break
        if chosen:
            break
    if chosen is None and plist and clist:         # fallback: closest in-frame CDS
        ps = plist[0][1]
        cands = [(abs(len(cs) // 3 - len(ps)), cs) for _, cs in clist if len(cs) % 3 == 0]
        if cands:
            cs = min(cands)[1]
            chosen = (translate(cs).rstrip("*"), cs)
    if chosen:
        selected[g] = chosen

usable = sorted(s for s in complete
                if all(str(orth[s][t]) in selected for t in TAXA.values()))
for sym in usable:
    with open(f"../data/sequences/{sym}.pep.fa", "w") as fp, \
         open(f"../data/sequences/{sym}.cds.fa", "w") as fc:
        for sp, tid in TAXA.items():
            p, c = selected[str(orth[sym][tid])]
            fp.write(f">{sp}\n{p}\n")
            fc.write(f">{sp}\n{c}\n")
json.dump({"usable": usable, "taxa": TAXA}, open("../data/usable_genes.json", "w"))
print(f"symbols resolved: {len(chick)} | complete orthology: {len(complete)} | usable: {len(usable)}")
