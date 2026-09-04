#!/usr/bin/env bash
# simcloud_fetch.sh -- wait for a brane grid batch, download every cell's
# output bundle, and reassemble the local data/N<N>/p<p8>/data.dat tree that
# tools/analyze.py --all and tools/heatmap.py expect.
#
# Usage:
#   CLUSTER=mr2 BATCH=<batch-id> bash cloud/simcloud_fetch.sh
#   bash cloud/simcloud_fetch.sh            # uses cloud/.last_batch from submit
#
set -euo pipefail

CLUSTER="${CLUSTER:-mr2}"
export SIMCLOUD_CLUSTER="$CLUSTER"
repo="$(cd "$(dirname "$0")/.." && pwd)"
BATCH="${BATCH:-$(cat "$repo/cloud/.last_batch" 2>/dev/null || true)}"
OUTDIR="${OUTDIR:-$repo/data}"

if [ -z "${BATCH:-}" ]; then
  echo "error: no batch id (set BATCH=... or run simcloud_submit.sh first)" >&2
  exit 1
fi

echo "=== waiting for batch $BATCH on $CLUSTER ==="
simcloud job wait --batch "$BATCH" --summary --poll-interval 30s || true

dl="$repo/cloud/_dl/$BATCH"
mkdir -p "$dl"
echo "=== downloading results -> $dl ==="
# Pull consoles + output bundles for every job in the batch.
simcloud job list --batch "$BATCH" -f id \
  | simcloud job download - --to "$dl"

echo "=== merging into $OUTDIR ==="
mkdir -p "$OUTDIR"
n=0
# Output arrives as job-output-bundle.tgz per job; extract them first.
while IFS= read -r tgz; do
  tar -xzf "$tgz" -C "$(dirname "$tgz")" 2>/dev/null || true
done < <(find "$dl" -type f -name 'job-output-bundle.tgz')
# Each bundle carries an out/N<N>/p<p8>/data.dat; hoist them into data/.
while IFS= read -r f; do
  rel="${f#*out/}"           # -> N<N>/p<p8>/data.dat
  dst="$OUTDIR/$rel"
  mkdir -p "$(dirname "$dst")"
  cp "$f" "$dst"
  n=$((n + 1))
done < <(find "$dl" -type f -path '*out/N*/p*/data.dat')

echo "=== merged $n cells into $OUTDIR ==="
if [ "$n" -eq 0 ]; then
  echo "warning: found no data.dat under $dl -- inspect it and job consoles:" >&2
  echo "  simcloud -c $CLUSTER job list --batch $BATCH" >&2
else
  echo "next: uv run tools/analyze.py --all && uv run tools/heatmap.py --replot-all"
fi
