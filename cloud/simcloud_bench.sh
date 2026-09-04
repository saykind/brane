#!/usr/bin/env bash
# simcloud_bench.sh -- runs INSIDE one Simcloud job (or locally). Measures how
# long one cell takes and how it scales with cores (replicas).
#
# The engine is replica-parallel: nt independent Markov chains run as OpenMP
# threads, each doing exactly `sweeps` sweeps (we set eps=0 to disable the
# adaptive stop, so the work is fixed and timings are comparable). We sweep:
#   * N  -> single-chain cost grows ~ N^4 (sets per-cell wall time)
#   * nt -> replicas; wall time should stay ~flat up to core count, so extra
#           replicas (= smaller error bars) are almost free until nt > cores.
#
# Output: a human table + machine-readable CSV (lines prefixed "CSV,").
set -euo pipefail

SRCDIR="${SRCDIR:-/brane}"
cd "$SRCDIR" 2>/dev/null || cd "$(cd "$(dirname "$0")/.." && pwd)"

command -v gcc >/dev/null 2>&1 || { command -v apt-get >/dev/null 2>&1 && \
  { apt-get update -y && apt-get install -y build-essential; }; } || true
make -B CC="${CC:-gcc}" >/dev/null    # -B: force rebuild (stale macOS binary ships in bundle)

NS="${NS:-32,48,64}"          # lattice sizes to probe
NTS="${NTS:-1,2,4,8,16}"      # replica/thread counts to probe
SWEEPS="${SWEEPS:-400}"        # fixed measurement sweeps (eps=0 -> no early stop)
THERM="${THERM:-40}"
P8="${P8:-0.4}"

cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc)
echo "=== brane scaling benchmark ==="
echo "host=$(hostname) cores=$cores  sweeps=$SWEEPS therm=$THERM p8=$P8"
printf '%-5s %-4s %-9s %-12s %-8s\n' N nt wall_s samples/s eff
echo "CSV,N,nt,wall_s,samples_per_s,efficiency"

IFS=',' read -ra NARR <<< "$NS"
IFS=',' read -ra TARR <<< "$NTS"
total=$(( ${#NARR[@]} * ${#TARR[@]} )); step=0

for N in "${NARR[@]}"; do
  base_tp=""                       # throughput at nt=1 for this N (for efficiency)
  for nt in "${TARR[@]}"; do
    tmp=$(mktemp)
    s=$(date +%s.%N)
    if ! ./brane "N=$N" "p8=$P8" "nt=$nt" "therm=$THERM" "sweeps=$SWEEPS" eps=0 \
                 out="$tmp" >/dev/null 2>&1; then
      echo "  run failed: N=$N nt=$nt" >&2; rm -f "$tmp"; continue
    fi
    e=$(date +%s.%N); rm -f "$tmp"
    line=$(awk -v s="$s" -v e="$e" -v nt="$nt" -v sw="$SWEEPS" -v base="$base_tp" 'BEGIN{
        w=e-s; if(w<=0)w=1e-6; tp=nt*sw/w;
        eff=(base=="")?1:tp/(base*nt);
        printf "%.3f %.1f %.2f\n", w, tp, eff }')
    read -r wall tp eff <<< "$line"
    [ -z "$base_tp" ] && base_tp="$tp"
    printf '%-5s %-4s %-9s %-12s %-8s\n' "$N" "$nt" "$wall" "$tp" "$eff"
    echo "CSV,$N,$nt,$wall,$tp,$eff"
    step=$(( step + 1 ))
    printf '  ... %d/%d\r' "$step" "$total" >&2
  done
done
echo >&2
echo "=== done ($total runs) ==="
