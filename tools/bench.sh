#!/bin/bash
# bench.sh -- compare the legacy single-chain code against the new
# replica-parallel engine, in single-mode-updates per second.
#
# The two use OpenMP differently: legacy parallelizes the inner O(L^2) loop of
# ONE chain (a fork/join per move); new runs one independent chain per thread.
# We therefore report raw update throughput, which is directly comparable.
set -e
cd "$(dirname "$0")/.."

N=${N:-40}
SWEEPS=${SWEEPS:-10}
P8=${P8:-0.4}
MAXNT=${MAXNT:-$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)}
L=$((2*N+1)); LL=$((L*L))

[ -x ./brane ] || make >/dev/null
[ -x ./legacy/a.out ] || ./legacy/build.sh >/dev/null

legacy_time () { # nt -> seconds of the measurement loop
  ./legacy/a.out N=$N p8=$P8 nt=$1 MTH=0 M=$SWEEPS 2>&1 \
    | grep 'time =' | tail -1 | grep -oE '[0-9.]+' | awk '{print $1*60}'
}
new_time () { # nt -> seconds of the sampling loop
  ./brane N=$N p8=$P8 nt=$1 therm=0 sweeps=$SWEEPS out=data/bench.dat 2>&1 \
    | grep '^time' | grep -oE 'time = [0-9.]+' | grep -oE '[0-9.]+'
}

echo "N=$N (L=$L), $SWEEPS sweeps, updates/sweep/chain=$LL, cores=$MAXNT"
echo
printf "%-16s %10s %8s %16s\n" "run" "time[s]" "chains" "updates/sec"
t=$(legacy_time 1);      printf "%-16s %10.2f %8d %16.0f\n" "legacy nt=1"      "$t" 1 "$(echo "$SWEEPS*$LL/$t"|bc -l)"
t=$(legacy_time $MAXNT); printf "%-16s %10.2f %8d %16.0f\n" "legacy nt=$MAXNT" "$t" 1 "$(echo "$SWEEPS*$LL/$t"|bc -l)"
t=$(new_time 1);         printf "%-16s %10.2f %8d %16.0f\n" "new nt=1"         "$t" 1 "$(echo "$SWEEPS*$LL/$t"|bc -l)"
t=$(new_time $MAXNT);    printf "%-16s %10.2f %8d %16.0f\n" "new nt=$MAXNT"    "$t" "$MAXNT" "$(echo "$MAXNT*$SWEEPS*$LL/$t"|bc -l)"
rm -f data/bench.dat
