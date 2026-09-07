#!/bin/bash
# run.sh -- multi-size finite-size sweep (thesis-style): run a range of lattice
# sizes N at one fixed coupling p8, then pool them to see the finite-size
# collapse of the inverse Green function G^-1(q). Each size runs in well under a
# minute on a laptop.
#
# Everything is overridable via env vars, e.g.:
#   P8=0.6 NT=8 SIZES="24 32 40" ./run.sh
set -e
P8=${P8:-0.4}
NT=${NT:-12}                       # throughput peaks near the P-core count (M4 Max)
THERM=${THERM:-80}
SWEEPS=${SWEEPS:-120}
SIZES=${SIZES:-"20 24 28 32 36"}

make -s                            # build ./brane if missing / out of date

for N in $SIZES; do
    # outdir=data lets the engine build its descriptive path
    # data/N<N>/p<p8>/fixed<sweeps>/therm..._nt..._it..._seed....dat
    # (one file per cell; never overwrites a different config).
    ./brane N=$N p8=$P8 nt=$NT therm=$THERM sweeps=$SWEEPS outdir=data
done

PDIR=$(printf 'p%.2f' "$P8")       # engine formats the coupling dir as p%.2f
echo
echo "=== finite-size collapse: every size on disk at $PDIR pooled ==="
# --combined pools every mode from all sizes into one G^-1(q) cloud + running
# exponent (writes plots/combined/), colored by N so the collapse is visible.
# The glob matches the descriptive layout, so any earlier runs / fetched cloud
# cells at this p8 are included too -- more sizes = a cleaner collapse.
uv run tools/analyze.py --combined "data/N*/${PDIR}/*/*.dat"
