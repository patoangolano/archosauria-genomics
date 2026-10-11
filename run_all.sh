#!/usr/bin/env bash
# Full pipeline, from public databases to figures.
#
# conda env create -f environment.yml && conda activate archosaur && ./run_all.sh
#
# Steps 1, 3 and 7 hit the NCBI APIs and take several minutes. Set
# NCBI_CONTACT_EMAIL (and CROSSREF_CONTACT_EMAIL for step 10) to identify
# yourself to those services, as they request.
set -euo pipefail
cd "$(dirname "$0")/src"

python 01_fetch_assemblies.py      # assembly + per-chromosome reports  -> Table 1
python 02_karyotype_gc.py          # microchromosomes, GC compartment   -> Table 2
python 03_fetch_orthologs.py       # 1:1 orthologs, protein + CDS
python 04_build_supermatrix.py     # MAFFT, codon back-translation
bash   05_phylogeny.sh             # IQ-TREE: ML search + fixed topology
python 06_rates_tajima.py          # lineage rates, Tajima RRT          -> Tables 3, 4
python 07_fetch_hox.py             # Hox coordinates + regions          -> Table 5
python 08_cne_scan.py              # conserved element scan             -> Table 6
python 09_figures.py               # Figures 1-3
# python 10_bibliography.py references.tsv   # optional: re-resolve DOIs
