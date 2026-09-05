#!/usr/bin/env bash
# simcloud_task.sh -- runs INSIDE one Simcloud batch job.
#
# Simcloud posts a batch of identical jobs (see simcloud_submit.sh); each job
# gets a unique index in $SC_BATCH_ID. We map that index to a single (N, p8)
# cell of the grid -- exactly like the SLURM array job (slurm_grid.sbatch) --
# build the engine, run that one cell using all allocated cores as replicas,
# and write the result under $OUTDIR so it can be shipped back as an output
# bundle. Error on each G(q) point falls ~ 1/sqrt(nt * sweeps).
#
# The repo is delivered as a bundle unpacked at /brane (see simcloud_submit.sh).
set -euo pipefail

# --- grid definition (must match simcloud_submit.sh) ------------------------
NS="${NS:-32,40,48,56,64,80,96,120}"           # lattice sizes (L=2N+1)
P8S="${P8S:-0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}"  # couplings (q8~p8)
THERM="${THERM:-300}"                           # thermalization sweeps (legacy MTH=300)
SWEEPS="${SWEEPS:-4000}"                         # measurement-sweep cap per cell
EPS="${EPS:-0.005}"                             # Delta2 rel-error convergence target
MINSW="${MINSW:-100}"                           # floor before convergence can trip
OUTDIR="${OUTDIR:-/out}"                         # collected via --output-to-bundle
NT="${NT:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc)}"  # replicas = cores
IT="${IT:-1}"                                    # inner threads/replica (large-N: nt=1 it=cores)
SRCDIR="${SRCDIR:-/brane}"                       # where the source bundle unpacks

# --- toolchain (a package bundle usually provides gcc/make already) ---------
if ! command -v gcc >/dev/null 2>&1; then
  # mr2 x86 containers have direct internet; fall back to apt if no bundle.
  if command -v apt-get >/dev/null 2>&1; then
    apt-get update -y && apt-get install -y build-essential
  fi
fi

# --- build ------------------------------------------------------------------
cd "$SRCDIR"
make -B CC=gcc    # -B: force rebuild (a stale macOS binary ships in the bundle)

# --- map batch index -> (N, p8) cell ---------------------------------------
IFS=',' read -ra NARR <<< "$NS"
IFS=',' read -ra PARR <<< "$P8S"
nP=${#PARR[@]}
i="${SC_BATCH_ID:-${SC_BATCH_INDEX:-0}}"        # CLI exposes the index as SC_BATCH_ID
N="${NARR[$(( i / nP ))]}"
P="${PARR[$(( i % nP ))]}"

# Let the engine build the descriptive path under $OUTDIR:
#   $OUTDIR/N<N>/p<p8>/<stop>/therm..._nt..._it..._seed....dat
echo ">>> batch index $i -> N=$N p8=$P nt=$NT it=$IT on $(hostname) ($(nproc) cores)"
./brane "N=$N" "p8=$P" "nt=$NT" "it=$IT" "therm=$THERM" "sweeps=$SWEEPS" \
        "eps=$EPS" "minsweeps=$MINSW" "outdir=$OUTDIR"
echo ">>> wrote under $OUTDIR/N$N/p$P/"
