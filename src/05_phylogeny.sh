#!/usr/bin/env bash
# Step 5 - ML tree plus fixed-topology branch lengths for three partitions.
#
# KMP_AFFINITY=disabled is load-bearing: IQ-TREE's OpenMP runtime aborts with
# "pthread_setaffinity_np() failed / Operation not permitted" inside containers
# that deny CPU-affinity syscalls.
set -euo pipefail
export KMP_AFFINITY=disabled OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
cd "$(dirname "$0")/.."
mkdir -p alignments/trees
T=alignments/tree_fixed.nwk

iqtree -s alignments/super.pep.phy -m LG+F+G4 -bb 1000 -alrt 1000 \
       -nt "$OMP_NUM_THREADS" -pre alignments/trees/pep_ml -quiet -redo
for P in "pep:LG+F+G4" "p12:GTR+F+G4" "p3:GTR+F+G4"; do
  PART="${P%%:*}"; MODEL="${P##*:}"
  iqtree -s "alignments/super.${PART}.phy" -te "$T" -m "$MODEL" \
         -nt "$OMP_NUM_THREADS" -pre "alignments/trees/${PART}_fix" -quiet -redo
done
for f in alignments/trees/*.treefile; do echo "== $f"; cat "$f"; echo; done
