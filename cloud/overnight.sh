#!/usr/bin/env bash
# overnight.sh -- launch the sparse-grid production run on M2 Ultra (mr2-as).
# Each cell is one job; a cell's output bundle is exported only when the job
# COMPLETES, so the timeout must exceed the run time (the engine also
# checkpoints data.dat every 60s as a backup).
#
# Grid: N in {180,200} x p8 in {0.3,0.4,0.5,0.6,0.7} = 10 cells, one job each.
# (N=140,160 were already run -- see data/; add them back via NS= if needed.)
# Layout: NT=16 replicas x IT=1 = 16 cores/cell. IMPORTANT: mr2-as jobs run
# ARM *Linux* (SMI ubuntu22.04-v1, libgomp), NOT macOS -- and intra-chain
# parallelism (IT>1) saturates at ~1.6x by IT=2 on this cloud (measured), so it
# never beats more replicas; keep IT=1. IT>1 helps more on macOS/libomp (the
# local M-series Mac). 16 replicas maximizes statistics (error ~ 1/sqrt(nt*sweeps)).
#
# Sizing (M2 Ultra ARM, ~N^4): N=160 measured ~58 s/sweep (1 thread); with
# therm=300 + sweeps=800 = 1100 sweeps the wall (= one chain's single-thread time,
# all 16 chains run in parallel) is N=160 ~17.8h, N=180 ~28.5h, N=200 ~43.5h.
# TIMEOUT=48h covers the slowest cell (N=200) with margin. (If quota max_timeout
# < 48h, lower TIMEOUT and drop N=200, or split the grid.)
#
# Usage:
#   bash cloud/overnight.sh                 # launch with defaults below
#   NS=140,160,180,200 bash cloud/overnight.sh   # override any knob
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
NS="${NS:-180,200}" P8S="${P8S:-0.3,0.4,0.5,0.6,0.7}" \
THERM="${THERM:-300}" SWEEPS="${SWEEPS:-800}" EPS="${EPS:-0}" MINSW="${MINSW:-100}" \
NT="${NT:-16}" IT="${IT:-1}" \
TAG="${TAG:-brane-overnight}" \
bash "$here/simcloud_submit.sh"
