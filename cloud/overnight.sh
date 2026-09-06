#!/usr/bin/env bash
# overnight.sh -- launch the sparse-grid production run on M2 Ultra (mr2-as).
# Each cell is one job; a cell's output bundle is exported only when the job
# COMPLETES, so the timeout must exceed the run time (the engine also
# checkpoints data.dat every 60s as a backup).
#
# Grid: N in {140,160,180,200} x p8 in {0.3,0.4,0.5,0.6,0.7} = 20 cells, one job
# each. Layout: NT=4 replicas x IT=4 inner threads = 16 cores/cell. IT>1
# (intra-chain) speeds a single chain on Apple Silicon (macOS/libomp); the
# tradeoff vs NT=16,IT=1 is FEWER replicas (4 -> larger per-mode error bars) in
# exchange for finishing the large-N cells sooner. If more statistics are needed,
# raise SWEEPS or run more seeds.
#
# Sizing (M2 Ultra, ~N^4): N=160 measured ~58 s/sweep single-thread; with
# therm=300 + sweeps=800 = 1100 sweeps the single-thread wall is N=160 ~17.8h,
# N=180 ~28.5h, N=200 ~43.5h. IT=4 cuts each by ~2-2.5x (macOS only). TIMEOUT=48h
# covers the worst cell (N=200) EVEN IF intra-chain parallelism yields no speedup,
# so no cell is lost to a timeout. (If quota max_timeout < 48h, lower TIMEOUT and
# drop N=200, or split the grid.)
#
# Usage:
#   bash cloud/overnight.sh                 # launch with defaults below
#   NS=140,160 IT=1 NT=16 bash cloud/overnight.sh   # override any knob
# Afterwards, pull results:
#   CLUSTER=mr2-as bash cloud/simcloud_fetch.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
scuser="$(simcloud -q -c mr2 user info 2>/dev/null | awk -F': *' '/^Username:/{print $2; exit}')"
[ -z "$scuser" ] && scuser="${USER}"

CLUSTER="${CLUSTER:-mr2-as}" \
OWNER="${OWNER:-hwt:atg:sph:$scuser}" \
NET="${NET:-e57cff0a-d781-4250-8ca5-065e283c8da1}" \
TOOLCHAIN="${TOOLCHAIN:-0}" \
CPUS="${CPUS:-16}" MEMORY="${MEMORY:-16}" DISK="${DISK:-30}" TIMEOUT="${TIMEOUT:-48h}" \
NS="${NS:-140,160,180,200}" P8S="${P8S:-0.3,0.4,0.5,0.6,0.7}" \
THERM="${THERM:-300}" SWEEPS="${SWEEPS:-800}" EPS="${EPS:-0}" MINSW="${MINSW:-100}" \
NT="${NT:-4}" IT="${IT:-4}" \
TAG="${TAG:-brane-overnight}" \
bash "$here/simcloud_submit.sh"
