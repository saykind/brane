#!/usr/bin/env bash
# simcloud_bench_run.sh -- launch the scaling benchmark on Simcloud.
# Interrupt-safe: uses `job post` (not `job run`) and a trap that cancels the
# remote job if this script is interrupted, so a Ctrl-C can't leave a box
# burning cores. Streams the console live while it runs.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
. "$here/simcloud_lib.sh"

CLUSTER="${CLUSTER:-mr2}"; export SIMCLOUD_CLUSTER="$CLUSTER"
SMI="${SMI:-ubuntu22.04-v1}"
CPUS="${CPUS:-32}"                 # cores on the box (mr2 node limit is 72); also caps nt
MEMORY="${MEMORY:-16}"
TIMEOUT="${TIMEOUT:-5m}"           # cap only; job exits as soon as the bench finishes
NET="${NET:-}"                     # Denali VPC network id (required on Apple Silicon)
TOOLCHAIN="${TOOLCHAIN:-1}"        # 1=ship build-essential bundle; 0=SMI already has gcc (AS)
# group quota is required for >20 cpus; derive the simcloud username (not $USER)
scuser="$(simcloud -q -c mr2 user info 2>/dev/null | awk -F': *' '/^Username:/{print $2; exit}')"
[ -z "$scuser" ] && scuser="${USER}"
OWNER="${OWNER:-hw:others:$scuser}"

# what to probe (passed through to simcloud_bench.sh)
NS="${NS:-24,32}"
NTS="${NTS:-1,2,4,8,16,32}"
SWEEPS="${SWEEPS:-100}"
THERM="${THERM:-20}"
P8="${P8:-0.4}"

repo="$(cd "$here/.." && pwd)"
echo "=== brane scaling benchmark -> Simcloud ($CLUSTER, $SMI, ${CPUS} cores) ==="
echo "probing N={$NS} nt={$NTS} at sweeps=$SWEEPS  (timeout cap $TIMEOUT, owner $OWNER)"

echo "--- source bundle ---"; src=$(sc_src_bundle "$repo"); echo "  $src"
bundles="$src"
if [ "$TOOLCHAIN" = "1" ]; then
  tool=$(sc_toolchain "$SMI"); echo "toolchain: $tool"; bundles="$tool,$src"
else
  echo "toolchain: (skipped -- SMI provides gcc)"
fi

owner_flag=(); [ -n "$OWNER" ] && owner_flag=(--owner "$OWNER")
vpc_flag=();   [ -n "$NET" ]   && vpc_flag=(--denali-vpc "ipv4_network_id=$NET")
cmd="NS=$NS NTS=$NTS SWEEPS=$SWEEPS THERM=$THERM P8=$P8 bash /brane/cloud/simcloud_bench.sh"

echo "--- posting job ---"
job=$(simcloud -q job post \
  --smi "$SMI" --cpus "$CPUS" --memory "$MEMORY" --timeout "$TIMEOUT" \
  --bundles "$bundles" \
  "${owner_flag[@]+"${owner_flag[@]}"}" \
  "${vpc_flag[@]+"${vpc_flag[@]}"}" \
  --command "$cmd")
echo "job: $job"

# cancel the remote job on any interrupt/exit-before-completion
cleanup() { simcloud -q job cancel "$job" >/dev/null 2>&1 || true; }
trap 'echo; echo "interrupted -- cancelling $job"; cleanup; exit 130' INT TERM

echo "--- waiting for start, then streaming ---"
simcloud job wait "$job" --status running --poll-interval 3s >&2 || true
simcloud job console --live "$job" 2>&1 \
  | grep -vE 'string_fortified|builtin|__glibc|__dest|\^~~|In file|In function|inlined|warning:'
trap - INT TERM
echo "=== job $job finished (status: $(simcloud -q job info -f '{{.Status}}' "$job" 2>/dev/null)) ==="

