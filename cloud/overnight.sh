#!/usr/bin/env bash
# overnight.sh -- launch the sparse-grid production run on M2 Ultra (mr2-as),
# sized to finish within an 8h window so every cell COMPLETES and every output
# bundle is saved (Simcloud exports --output-to-bundle on completion). The
# engine additionally checkpoints each cell's data.dat every 60s as a backup,
# so even a killed cell keeps its latest data.
#
# Grid: N in {140,160} x p8 in {0.3,0.4,0.5,0.6,0.7} = 10 cells, one job each,
# 16 replicas/cell. Calibrated on M2 Ultra: N=100 ~ 9.7 s/sweep, scaling ~N^4,
# so N=140 ~ 37 s/sweep and N=160 ~ 64 s/sweep. With therm=300 sweeps=800:
# N=140 ~ 11.4h, N=160 ~ 19.5h wall (the slowest cell). 24h timeout leaves
# margin so every cell COMPLETES (Simcloud exports the output bundle only on
# completion; a timeout-killed cell loses its data). Total wall ~ slowest cell.
#
# Usage:
#   bash cloud/overnight.sh                 # launch with defaults below
#   SWEEPS=1200 bash cloud/overnight.sh     # override any knob
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
CPUS="${CPUS:-16}" MEMORY="${MEMORY:-16}" DISK="${DISK:-30}" TIMEOUT="${TIMEOUT:-24h}" \
NS="${NS:-140,160}" P8S="${P8S:-0.3,0.4,0.5,0.6,0.7}" \
THERM="${THERM:-300}" SWEEPS="${SWEEPS:-800}" EPS="${EPS:-0}" MINSW="${MINSW:-100}" \
NT="${NT:-16}" IT="${IT:-1}" \
TAG="${TAG:-brane-overnight}" \
bash "$here/simcloud_submit.sh"
