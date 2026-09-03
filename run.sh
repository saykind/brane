#!/bin/bash
# run.sh -- multi-size sweep, mirroring the thesis workflow (legacy/save.sh).
# Each size runs in well under a minute; results land in data/N=<N>.dat and
# can be overlaid to check finite-size collapse of the inverse Green function.
set -e
P8=${P8:-0.4}
NT=${NT:-$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)}
mkdir -p data
for N in 20 24 28 32 36; do
    ./brane N=$N p8=$P8 nt=$NT therm=80 sweeps=120 out=data/N=$N.dat
done
echo
echo "Analyze the largest size:"
python3 tools/analyze.py data/N=36.dat
