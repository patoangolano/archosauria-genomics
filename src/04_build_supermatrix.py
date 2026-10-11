"""MAFFT protein alignment -> codon back-translation -> trimmed supermatrices."""
import json, os, subprocess, sys
from collections import OrderedDict

SP = ["Gallus_gallus", "Alligator_mississippiensis", "Crocodylus_porosus",
      "Gavialis_gangeticus", "Chrysemys_picta", "Anolis_carolinensis", "Homo_sapiens"]
MIN_OCC = 6          # protein column kept if >= 6/7 taxa non-gap
os.makedirs("../alignments", exist_ok=True)

def read_fasta(p):
    d, name, buf = OrderedDict(), None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            if name: d[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else: buf.append(line)
    if name: d[name] = "".join(buf)
    return d

genes = sorted(json.load(open("../data/usable_genes.json"))["usable"])
pep_blocks, cod_blocks, parts, dropped = [], [], [], []

for gsym in genes:
    pep_in, cds = f"../data/sequences/{gsym}.pep.fa", read_fasta(f"../data/sequences/{gsym}.cds.fa")
    pep_aln_f = f"../alignments/{gsym}.pep.aln"
    with open(pep_aln_f, "w") as fh:
        subprocess.run(["mafft", "--auto", "--quiet", "--anysymbol", pep_in],
                       stdout=fh, check=True)
    aln = read_fasta(pep_aln_f)
    if set(aln) != set(SP):
        dropped.append((gsym, "taxa")); continue
    L = len(next(iter(aln.values())))
    keep = [i for i in range(L) if sum(aln[s][i] != "-" for s in SP) >= MIN_OCC]
    if len(keep) < 50:
        dropped.append((gsym, "short")); continue
    pep_t = {s: "".join(aln[s][i] for i in keep) for s in SP}
    # back-translate: walk each taxon's ungapped CDS codon by codon
    cod_t, bad = {}, False
    for s in SP:
        nt = cds[s].upper()
        if len(nt) % 3: nt = nt[:len(nt) - len(nt) % 3]
        codons = [nt[i:i+3] for i in range(0, len(nt), 3)]
        if codons and codons[-1] in ("TAA", "TAG", "TGA"): codons = codons[:-1]
        out, k = [], 0
        for ch in aln[s]:
            if ch == "-": out.append("---")
            else:
                out.append(codons[k] if k < len(codons) else "NNN"); k += 1
        if k > len(codons): bad = True
        cod_t[s] = "".join(out[i] for i in keep)
    if bad:
        dropped.append((gsym, "cds_mismatch")); continue
    pep_blocks.append((gsym, pep_t)); cod_blocks.append((gsym, cod_t))
    parts.append((gsym, len(keep)))

def write_phy(blocks, path, factor):
    cat = {s: "".join(b[s] for _, b in blocks) for s in SP}
    n = len(next(iter(cat.values())))
    with open(path, "w") as fh:
        fh.write(f"{len(SP)} {n}\n")
        for s in SP: fh.write(f"{s}  {cat[s]}\n")
    return n

n_pep = write_phy(pep_blocks, "../alignments/super.pep.phy", 1)
n_cod = write_phy(cod_blocks, "../alignments/super.codon.phy", 3)
# codon-position partitions
cat_cod = {s: "".join(b[s] for _, b in cod_blocks) for s in SP}
for tag, sel in [("p12", lambda i: i % 3 != 2), ("p3", lambda i: i % 3 == 2)]:
    sub = {s: "".join(c for i, c in enumerate(v) if sel(i)) for s, v in cat_cod.items()}
    with open(f"../alignments/super.{tag}.phy", "w") as fh:
        fh.write(f"{len(SP)} {len(next(iter(sub.values())))}\n")
        for s in SP: fh.write(f"{s}  {sub[s]}\n")

json.dump({"genes_used": [g for g, _ in parts], "partitions": parts,
           "n_aa": n_pep, "n_nt": n_cod, "dropped": dropped},
          open("../alignments/supermatrix_info.json", "w"))
print(f"genes used: {len(parts)} | dropped: {dropped}")
print(f"protein supermatrix: {n_pep} aa | codon supermatrix: {n_cod} nt")
