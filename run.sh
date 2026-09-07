#!/bin/bash
# run.sh -- multi-size sweep: run a range of lattice sizes N at one fixed p8.
# Results land in the descriptive layout under data/. Analyze separately, e.g.
#   uv run tools/analyze.py --combined 'data/N*/p0.40/*/*.dat'
#
# Override any knob via env vars:  P8=0.6 NT=8 SIZES="24 32 40" ./run.sh
set -e
P8=${P8:-0.4}
NT=${NT:-12}
THERM=${THERM:-80}
SWEEPS=${SWEEPS:-120}
SIZES=${SIZES:-"20 24 28 32 36"}

make -s
for N in $SIZES; do
    ./brane N=$N p8=$P8 nt=$NT therm=$THERM sweeps=$SWEEPS outdir=data
done
