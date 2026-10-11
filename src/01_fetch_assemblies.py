"""Step 1 - reference assembly reports and per-chromosome sequence reports.

Writes data/ncbi_dataset_reports.json and data/ncbi_sequence_reports.json.gz
plus results/tabela1_estatisticas_genomicas.csv.
"""
import gzip, json, os
import pandas as pd
from ncbi import BASE, ASSEMBLY_TAXA, get, polite_sleep

os.makedirs("../data", exist_ok=True)
os.makedirs("../results", exist_ok=True)

reports = {}
for name, tid in ASSEMBLY_TAXA.items():
    reports[name] = get(f"{BASE}/genome/taxon/{tid}/dataset_report",
                        {"filters.reference_only": "true",
                         "filters.assembly_version": "current", "page_size": 20})
    polite_sleep(0.4)
json.dump(reports, open("../data/ncbi_dataset_reports.json", "w"))

rows = []
for name, d in reports.items():
    r = d["reports"][0]
    ai, st = r.get("assembly_info", {}), r.get("assembly_stats", {})
    gc = (r.get("annotation_info", {}) or {}).get("stats", {}).get("gene_counts", {})
    rows.append({
        "species": name.replace("_", " "), "accession": r["accession"],
        "assembly": ai.get("assembly_name"), "level": ai.get("assembly_level"),
        "release_date": ai.get("release_date"),
        "size_Gb": round(int(st.get("total_sequence_length", 0)) / 1e9, 3),
        "gc_percent": st.get("gc_percent"), "contig_n50": st.get("contig_n50"),
        "scaffold_n50": st.get("scaffold_n50"), "n_scaffolds": st.get("number_of_scaffolds"),
        "genes_total": gc.get("total"), "protein_coding": gc.get("protein_coding"),
        "non_coding": gc.get("non_coding"), "pseudogene": gc.get("pseudogene"),
        "ann_name": (r.get("annotation_info", {}) or {}).get("name")})
pd.DataFrame(rows).to_csv("../results/tabela1_estatisticas_genomicas.csv", index=False)

seqreps = {}
for name in ["Gallus_gallus", "Alligator_mississippiensis", "Chrysemys_picta",
             "Taeniopygia_guttata", "Anolis_carolinensis"]:
    acc = reports[name]["reports"][0]["accession"]
    out, token = [], None
    while True:
        p = {"page_size": 1000}
        if token:
            p["page_token"] = token
        d = get(f"{BASE}/genome/accession/{acc}/sequence_reports", p)
        out += d.get("reports", [])
        token = d.get("next_page_token")
        if not token:
            break
    seqreps[name.replace("_", " ")] = out
    polite_sleep(0.4)

with gzip.open("../data/ncbi_sequence_reports.json.gz", "wt") as fh:
    json.dump(seqreps, fh)
print(f"assemblies: {len(rows)} | sequence reports: {sum(len(v) for v in seqreps.values())}")
