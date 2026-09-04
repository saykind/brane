#!/bin/bash
#SBATCH --job-name=brane-grid
#SBATCH --output=brane_%A_%a.out
#SBATCH --array=0-63            # one task per (N,p8) cell; adjust to N_cells-1
#SBATCH --cpus-per-task=32      # replicas per cell = cores on the node slice
#SBATCH --time=02:00:00
#SBATCH --mem=4G                # engine is memory-light
#
# SLURM array job: each task runs ONE (N, p8) cell using all its allocated
# cores as replicas. Submit with:  sbatch cloud/slurm_grid.sbatch
#
# The array index is mapped to a (N, p8) pair below. Edit NS/P8S and set
# --array=0-(len(NS)*len(P8S)-1) to match. Results land in the shared
# data/N<N>/p<p8>/data.dat tree, ready for tools/analyze.py --all locally.
set -euo pipefail

NS=(32 40 48 56 64 80 96 120)
P8S=(0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

nP=${#P8S[@]}
i=${SLURM_ARRAY_TASK_ID:-0}
N=${NS[$(( i / nP ))]}
P=${P8S[$(( i % nP ))]}
NT=${SLURM_CPUS_PER_TASK:-8}

module load gcc 2>/dev/null || true          # site-dependent; ignore if absent
make CC=gcc

out=$(printf 'data/N%s/p%.2f/data.dat' "$N" "$P")
mkdir -p "$(dirname "$out")"
echo ">>> task $i: N=$N p8=$P nt=$NT on $(hostname)"
./brane "N=$N" "p8=$P" "nt=$NT" therm=100 sweeps=4000 eps=0.005 minsweeps=100 "out=$out"
