#!/usr/bin/env bash
# overnight.sh -- launch the sparse-grid production run on M2 Ultra (mr2-as),
# sized to finish within an 8h window so every cell COMPLETES and every output
# bundle is saved (Simcloud exports --output-to-bundle on completion). The
# engine additionally checkpoints each cell's data.dat every 60s as a backup,
# so even a killed cell keeps its latest data.
#
# Grid: N in {80,90,100,110,120} x p8 in {0.3,0.4,0.5} = 15 cells, one job each,
# 16 replicas/cell. Calibrated on M2 Ultra: N=100 ~ 9.7 s/sweep, N=120 ~ 20
# s/sweep, so therm=100 sweeps=800 => N=120 ~ 6h wall (the slowest cell);
# smaller N finish sooner. Total wall ~ slowest cell (cells run in parallel).
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
CPUS="${CPUS:-16}" MEMORY="${MEMORY:-16}" DISK="${DISK:-30}" TIMEOUT="${TIMEOUT:-8h}" \
NS="${NS:-80,90,100,110,120}" P8S="${P8S:-0.3,0.4,0.5}" \
THERM="${THERM:-300}" SWEEPS="${SWEEPS:-800}" EPS="${EPS:-0}" MINSW="${MINSW:-100}" \
NT="${NT:-16}" IT="${IT:-1}" \
TAG="${TAG:-brane-overnight}" \
bash "$here/simcloud_submit.sh"
