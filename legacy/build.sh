#!/bin/bash
# Minimal fix so the original code builds on macOS (Apple clang needs libomp).
# The only change vs. the thesis-era script is the OpenMP incantation and
# matching -O3 so speed comparisons reflect the algorithm, not compiler flags.
# (The original used "-g" with no optimization.)
set -e
cd "$(dirname "$0")"

OMP_CFLAGS="-fopenmp"
OMP_LDLIBS="-fopenmp"
if [ "$(uname -s)" = "Darwin" ]; then
    LIBOMP="$(brew --prefix libomp 2>/dev/null)"
    if [ -n "$LIBOMP" ]; then
        OMP_CFLAGS="-Xpreprocessor -fopenmp -I$LIBOMP/include"
        OMP_LDLIBS="-L$LIBOMP/lib -lomp"
    else
        echo "warning: libomp not found; run 'brew install libomp'" >&2
    fi
fi

OPT="-O3 -march=native -ffast-math -funroll-loops"

gcc -std=gnu11 -Wall $OPT $OMP_CFLAGS main.c storage.c simulate.c -lm $OMP_LDLIBS -o a.out
gcc -std=gnu11 -Wall $OPT $OMP_CFLAGS main.c storage.c simexact.c -lm $OMP_LDLIBS -o e.out
echo "built legacy/a.out and legacy/e.out"
