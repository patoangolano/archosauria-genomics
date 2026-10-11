"""Step 6 - per-lineage substitution distances and Tajima's relative-rate test.

Writes results/tabela_taxas_lineagem.csv, results/tabela_tajima_rrt.csv
and data/rate_results.json.
"""
import json
import numpy as np, pandas as pd
from Bio import Phylo
from scipy.stats import chi2 as chi2_dist

ARCHOSAURS = ["Gallus_gallus", "Alligator_mississippiensis",
              "Crocodylus_porosus", "Gavialis_gangeticus"]
CROCS = ARCHOSAURS[1:]


def read_phylip(path):
    out = {}
    for line in open(path).read().split("\n")[1:]:
        if line.strip():
            name, seq = line.split()
            out[name] = seq
    return out


def distances_from_archosaur_mrca(treefile):
    """Branch lengths measured from the bird/crocodilian split: same elapsed
    time for every tip, so the two lineages are directly comparable."""
    t = Phylo.read(treefile, "newick")
    t.root_with_outgroup({"name": "Homo_sapiens"})
    mrca = t.common_ancestor({"name": "Gallus_gallus"},
                             {"name": "Alligator_mississippiensis"})
    return {s: t.distance(mrca, s) for s in ARCHOSAURS}


rows = []
for label, path in [("Proteina (LG+F+G4)", "../alignments/trees/pep_fix.treefile"),
                    ("Codon pos 1+2 (GTR+F+G4)", "../alignments/trees/p12_fix.treefile"),
                    ("Codon pos 3 (GTR+F+G4)", "../alignments/trees/p3_fix.treefile")]:
    d = distances_from_archosaur_mrca(path)
    croc = float(np.mean([d[s] for s in CROCS]))
    rows.append({"partition": label, **{s.split("_")[0]: round(d[s], 4) for s in ARCHOSAURS},
                 "crocodylia_mean": round(croc, 4),
                 "ratio_bird_croc": round(d["Gallus_gallus"] / croc, 2)})
rate_df = pd.DataFrame(rows)
rate_df.to_csv("../results/tabela_taxas_lineagem.csv", index=False)


def tajima_rrt(a, b, outgroup, seqs, site_filter=None):
    """Tajima (1993) 1D relative-rate test. m1/m2 count sites where exactly one
    ingroup differs from the outgroup; equal rates predict m1 == m2."""
    A, B, O = seqs[a], seqs[b], seqs[outgroup]
    m1 = m2 = 0
    for i in range(len(A)):
        if site_filter and not site_filter(i):
            continue
        x, y, z = A[i], B[i], O[i]
        if "-" in (x, y, z) or "N" in (x, y, z):
            continue
        if x != z and y == z:
            m1 += 1
        elif y != z and x == z:
            m2 += 1
    chi = (m1 - m2) ** 2 / (m1 + m2) if (m1 + m2) else 0.0
    return m1, m2, chi, float(chi2_dist.sf(chi, 1))


codon = read_phylip("../alignments/super.codon.phy")
taj = []
for out in ["Chrysemys_picta", "Anolis_carolinensis", "Homo_sapiens"]:
    for label, filt in [("all sites", None), ("3rd position", lambda i: i % 3 == 2)]:
        m1, m2, chi, p = tajima_rrt("Gallus_gallus", "Alligator_mississippiensis",
                                    out, codon, filt)
        taj.append({"outgroup": out.replace("_", " "), "sites": label,
                    "m1_Gallus": m1, "m2_Alligator": m2, "m1_over_m2": round(m1 / m2, 2),
                    "chi2": round(chi, 1), "P": f"{p:.2e}"})
taj_df = pd.DataFrame(taj)
taj_df.to_csv("../results/tabela_tajima_rrt.csv", index=False)

json.dump({"rates": rows, "tajima": taj}, open("../data/rate_results.json", "w"))
print(rate_df.to_string(index=False)); print(); print(taj_df.to_string(index=False))
