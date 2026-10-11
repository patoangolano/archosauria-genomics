"""Step 7 - Hox cluster coordinates and genomic regions for the CNE scan.

Writes data/hox_locations.json, results/tabela5_hox_sintenia.csv and the
per-cluster FASTA under data/hox/ consumed by step 8.
"""
import json, os
import pandas as pd
from ncbi import BASE, get, efetch_region, get_zip, parse_fasta, gene_id_of, polite_sleep
import urllib.parse

FOCAL = {"Gallus_gallus": 9031, "Alligator_mississippiensis": 8496, "Crocodylus_porosus": 8502}
HOX = [f"HOX{c}{n}" for c in "ABCD" for n in range(1, 14)]
FLANK = 20_000
os.makedirs("../data/hox", exist_ok=True)

rows = []
for sp, tid in FOCAL.items():
    for i in range(0, len(HOX), 8):
        batch = HOX[i:i + 8]
        try:
            d = get(f"{BASE}/gene/symbol/{','.join(batch)}/taxon/{tid}", {"page_size": 100})
        except Exception:
            continue
        for r in d.get("reports", []):
            g = r.get("gene", {})
            if g.get("symbol", "").upper() not in [b.upper() for b in batch]:
                continue
            for a in g.get("annotations", []) or []:
                for loc in a.get("genomic_locations", []) or []:
                    gr = loc.get("genomic_range", {})
                    rows.append({"species": sp, "symbol": g["symbol"].upper(),
                                 "gene_id": g["gene_id"], "seq": loc.get("sequence_name"),
                                 "acc": loc.get("genomic_accession_version"),
                                 "begin": int(gr.get("begin", 0)), "end": int(gr.get("end", 0)),
                                 "strand": gr.get("orientation")})
        polite_sleep()

json.dump(rows, open("../data/hox_locations.json", "w"))
hox = pd.DataFrame(rows).drop_duplicates(subset=["species", "symbol"])
hox["cluster"] = hox.symbol.str[3]

spans = (hox.groupby(["species", "cluster", "acc", "seq"])
            .agg(start=("begin", "min"), stop=("end", "max"), n_genes=("symbol", "size"))
            .reset_index())
spans["span_kb"] = ((spans.stop - spans.start) / 1000).round(1)
spans.to_csv("../results/tabela5_hox_sintenia.csv", index=False)

for _, r in spans.iterrows():
    seq = efetch_region(r.acc, int(r.start) - FLANK, int(r.stop) + FLANK)
    with open(f"../data/hox/{r.species}_HOX{r.cluster}.fa", "w") as fh:
        fh.write(f">{r.species}_HOX{r.cluster} {r.acc}:{max(1, int(r.start) - FLANK)}-"
                 f"{int(r.stop) + FLANK}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")
    polite_sleep(0.4)

# Hox CDS, used in step 8 to mask coding exons inside each cluster
ids = sorted(set(hox.gene_id.astype(str)))
cds = {}
for i in range(0, len(ids), 30):
    url = (f"{BASE}/gene/id/{','.join(ids[i:i + 30])}/download?"
           + urllib.parse.urlencode([("include_annotation_type", "FASTA_CDS")]))
    z = get_zip(url)
    for nm in z.namelist():
        if nm.endswith("cds.fna"):
            for h, s in parse_fasta(z.read(nm).decode()):
                g = gene_id_of(h)
                if g and (g not in cds or len(s) > len(cds[g])):
                    cds[g] = s
    polite_sleep(0.4)

idx = hox.set_index(hox.gene_id.astype(str))
for sp in FOCAL:
    with open(f"../data/hox/{sp}_cds.fa", "w") as fh:
        for g, s in cds.items():
            if g in idx.index and idx.loc[g, "species"] == sp:
                fh.write(f">{idx.loc[g, 'symbol']}_{g}\n{s}\n")
print(spans.to_string(index=False))
