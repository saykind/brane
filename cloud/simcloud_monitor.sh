#!/usr/bin/env bash
# simcloud_monitor.sh -- watch a brane batch with a live progress bar.
#
# A friendlier replacement for `simcloud job wait --summary`: polls the batch
# and renders a single updating line showing how many jobs are done, running,
# and failed, plus a progress bar, until every job reaches a terminal state.
#
# Usage:
#   CLUSTER=mr2-as bash cloud/simcloud_monitor.sh                 # uses cloud/.last_batch
#   CLUSTER=mr2-as BATCH=<batch-id> bash cloud/simcloud_monitor.sh
#   CLUSTER=mr2-as bash cloud/simcloud_monitor.sh <batch-id> [poll-seconds]
set -euo pipefail

CLUSTER="${CLUSTER:-mr2-as}"
export SIMCLOUD_CLUSTER="$CLUSTER"
repo="$(cd "$(dirname "$0")" && pwd)/.."
BATCH="${1:-${BATCH:-$(cat "$repo/cloud/.last_batch" 2>/dev/null || true)}}"
POLL="${2:-${POLL:-30}}"

if [ -z "${BATCH:-}" ]; then
  echo "error: no batch id (pass one, set BATCH=, or run simcloud_submit.sh first)" >&2
  exit 1
fi

# Render a single updating progress-bar line to stderr.
bar() {
  local cur=$1 tot=$2 label="${3:-}" width="${SC_BAR_WIDTH:-30}"
  local filled=0 pct=0
  if [ "$tot" -gt 0 ]; then filled=$(( cur * width / tot )); pct=$(( cur * 100 / tot )); fi
  [ "$filled" -gt "$width" ] && filled=$width
  local b e
  printf -v b '%*s' "$filled" ''; b=${b// /#}
  printf -v e '%*s' "$(( width - filled ))" ''; e=${e// /.}
  printf '\r[%s%s] %3d%% %s\033[K' "$b" "$e" "$pct" "$label" >&2
}

echo "=== monitoring batch $BATCH on $CLUSTER (poll ${POLL}s) ===" >&2
while :; do
  # One `job list` call; tally states so the poll cost is a single request.
  statuses=$(simcloud -q job list --batch "$BATCH" -f '{{.Status}}' 2>/dev/null || true)
  total=$(printf '%s\n' "$statuses" | grep -c . || true)
  done_n=$(printf '%s\n' "$statuses" | grep -ciE 'complete|error|cancel|fail|timeout' || true)
  fail_n=$(printf '%s\n' "$statuses" | grep -ciE 'error|cancel|fail|timeout' || true)
  run_n=$(printf '%s\n' "$statuses" | grep -ciE 'run|start|active' || true)
  total=${total:-0}; done_n=${done_n:-0}; fail_n=${fail_n:-0}; run_n=${run_n:-0}
  bar "$done_n" "$total" "done $done_n/$total  running $run_n  failed $fail_n"
  { [ "$total" -gt 0 ] && [ "$done_n" -ge "$total" ]; } && { echo >&2; break; }
  sleep "$POLL"
done

echo "=== batch $BATCH: all $total jobs terminal ($fail_n failed) ===" >&2
echo "next: CLUSTER=$CLUSTER BATCH=$BATCH bash cloud/simcloud_fetch.sh" >&2
