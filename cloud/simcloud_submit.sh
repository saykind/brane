#!/usr/bin/env bash
# simcloud_submit.sh -- package the repo + submit the whole (N, p8) grid as a
# Simcloud batch, one job per cell (mirrors slurm_grid.sbatch, but on ACS).
#
# Run this on your laptop (needs the `simcloud` CLI + AppleConnect). It:
#   1. builds a source bundle from the repo (no data/plots/binaries),
#   2. ensures a toolchain (build-essential) package bundle for the SMI,
#   3. posts a batch of len(NS)*len(P8S) jobs running cloud/simcloud_task.sh,
#      each cell using --cpus cores as replicas and saving its data.dat to an
#      output bundle (tagged so simcloud_fetch.sh can pull them back).
#
# Everything is env-configurable (defaults below). Typical use:
#   # x86 (mr2): slow cores, big node -- fine for statistics, poor for large N
#   OWNER=hw:others:$USER CPUS=16 bash cloud/simcloud_submit.sh
#
#   # M2 Ultra (mr2-as): ~3.6x faster/core, near-linear scaling -- RECOMMENDED.
#   # gcc ships in the SMI (TOOLCHAIN=0); needs group quota + Denali VPC net id.
#   CLUSTER=mr2-as OWNER=hwt:atg:sph:$USER TOOLCHAIN=0 CPUS=16 \
#     NET=e57cff0a-d781-4250-8ca5-065e283c8da1 bash cloud/simcloud_submit.sh
#
# See cloud/SIMCLOUD.md for cluster comparison, VPC net ids, and scaling data.
set -euo pipefail

# --- cluster / resources ----------------------------------------------------
CLUSTER="${CLUSTER:-mr2}"                        # mr2=x86 (internet); mr2-as=M2 Ultra
export SIMCLOUD_CLUSTER="$CLUSTER"
SMI="${SMI:-ubuntu22.04-v1}"                     # Linux image (see: simcloud cluster smis)
CPUS="${CPUS:-8}"                                # cores per cell == replicas per cell
MEMORY="${MEMORY:-8}"                            # GB; engine is memory-light
DISK="${DISK:-20}"                              # GB
TIMEOUT="${TIMEOUT:-4h}"                         # per-cell wall-clock cap
OWNER="${OWNER:-}"                              # group quota fqn, e.g. hwt:atg:sph:$USER
NET="${NET:-}"                                  # Denali VPC network id (REQUIRED on -as)
TOOLCHAIN="${TOOLCHAIN:-1}"                      # 1=ship build-essential bundle; 0=SMI has gcc (AS)

# --- grid definition (must match simcloud_task.sh) --------------------------
NS="${NS:-32,40,48,56,64,80,96,120}"
P8S="${P8S:-0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}"
export NS P8S
export THERM="${THERM:-300}" SWEEPS="${SWEEPS:-4000}" EPS="${EPS:-0.005}" MINSW="${MINSW:-100}"
# replicas (NT) x inner-threads (IT) per cell. Default: NT=CPUS replicas, IT=1
# (statistics). For large-N reach set NT=1 IT=$CPUS (fewer replicas, one fast
# chain -- intra-chain parallelism, engine 'it=' knob, wins for N>=~70).
NT="${NT:-$CPUS}"; IT="${IT:-1}"

TAG="${TAG:-brane-grid}"                         # output bundles tagged with this

# --- derive cell count ------------------------------------------------------
IFS=',' read -ra NARR <<< "$NS"
IFS=',' read -ra PARR <<< "$P8S"
COUNT=$(( ${#NARR[@]} * ${#PARR[@]} ))

repo="$(cd "$(dirname "$0")/.." && pwd)"
echo "=== brane -> Simcloud ($CLUSTER, $SMI) ==="
echo "grid: ${#NARR[@]} N x ${#PARR[@]} p8 = $COUNT cells; $CPUS cpus/cell (=${CPUS} replicas)"
echo "params: therm=$THERM sweeps=$SWEEPS eps=$EPS minsweeps=$MINSW timeout=$TIMEOUT"
echo "concurrent CPU demand if all run at once: $(( COUNT * CPUS )) (check your quota)"
[ -n "$OWNER" ] && echo "owner (quota): $OWNER"

# --- 1. source bundle -------------------------------------------------------
# Ship only what the build needs; exclude big dirs. NOTE: --exclude takes ONE
# regex matched against the ABSOLUTE path (which contains the repo dir name
# "brane"), so we anchor on big data/output dirs, not on "brane" itself. The
# small stale binaries are harmless -- the task rebuilds with `make`.
echo "--- creating source bundle (unpacks at /brane) ---"
srcbundle=$(simcloud -q bundle create -C "$repo" . --root /brane \
  --tag "brane-src" \
  --exclude '(^|/)(data|plots|eta|example_data|legacy|\.git|\.venv|__pycache__)(/|$)')
echo "source bundle: $srcbundle"

# --- 2. toolchain package bundle (built once, reused via tag) ---------------
bundles="$srcbundle"
if [ "$TOOLCHAIN" = "1" ]; then
  tcname="brane-toolchain-${SMI}"
  toolbundle=$(simcloud -q bundle find --tag "$tcname" 2>/dev/null || true)
  if [ -z "${toolbundle:-}" ]; then
    echo "--- building toolchain bundle (build-essential for $SMI) ---"
    toolbundle=$(simcloud -q bundle package --apt build-essential --smi "$SMI" --tag "$tcname")
  fi
  echo "toolchain bundle: $toolbundle"
  bundles="$toolbundle,$srcbundle"
else
  echo "toolchain: (skipped -- SMI provides gcc)"
fi

# --- 3. post the batch ------------------------------------------------------
# There is no --env flag; pass grid params by prefixing the command (values
# have no spaces, so this is safe). Each job derives its cell from $SC_BATCH_ID.
owner_flag=(); [ -n "$OWNER" ] && owner_flag=(--owner "$OWNER")
vpc_flag=();   [ -n "$NET" ]   && vpc_flag=(--denali-vpc "ipv4_network_id=$NET")
envprefix="NS=$NS P8S=$P8S THERM=$THERM SWEEPS=$SWEEPS EPS=$EPS MINSW=$MINSW NT=$NT IT=$IT OUTDIR=/out"
echo "--- posting batch of $COUNT jobs ---"
batch=$(simcloud -q batch post \
  --count "$COUNT" \
  --smi "$SMI" \
  --cpus "$CPUS" --memory "$MEMORY" --disk "$DISK" \
  --timeout "$TIMEOUT" \
  --bundles "$bundles" \
  --output /out --output-to-bundle \
  --output-bundle-tags "$TAG" \
  "${owner_flag[@]+"${owner_flag[@]}"}" \
  "${vpc_flag[@]+"${vpc_flag[@]}"}" \
  --command "$envprefix bash /brane/cloud/simcloud_task.sh")

echo
echo "=== batch submitted: $batch ==="
echo "monitor : simcloud -c $CLUSTER job wait --batch $batch --summary --poll-interval 30s"
echo "list    : simcloud -c $CLUSTER job list --batch $batch"
echo "fetch   : CLUSTER=$CLUSTER BATCH=$batch bash cloud/simcloud_fetch.sh"
echo
# stash for the fetch script's convenience
echo "$batch" > "$repo/cloud/.last_batch"
