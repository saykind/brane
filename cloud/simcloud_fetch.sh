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
echo "=== downloading output bundles -> $dl ==="
# Fetch the per-cell OUTPUT BUNDLES through the bundle service (reliable), NOT
# `job download` -- the latter also pulls console logs directly from each VM,
# which needs DC-VPN/direct reachability and times out on Apple Silicon. Output
# bundles are auto-tagged with the batch id, so filter on that.
ids=$(simcloud -q bundle list --tags "$BATCH" -f '{{.ID}}' 2>/dev/null | grep -E '^bundle-' || true)
if [ -z "$ids" ]; then
  echo "error: no output bundles tagged $BATCH (jobs not done, or none produced output)" >&2
  echo "  check: simcloud -c $CLUSTER job list --batch $BATCH" >&2
  exit 1
fi
for b in $ids; do
  simcloud -q bundle download "$b" --to "$dl" >/dev/null 2>&1 || echo "warn: download failed for $b" >&2
done

echo "=== merging into $OUTDIR ==="
mkdir -p "$OUTDIR"
n=0
# Each downloaded bundle is a gzip tarball (named bundle-XXX) of the job's /out
# tree; extract them, then hoist out/N<N>/p<p8>/data.dat into data/.
while IFS= read -r tgz; do
  tar -xzf "$tgz" -C "$dl" 2>/dev/null || true
done < <(find "$dl" -maxdepth 1 -type f -name 'bundle-*')
while IFS= read -r f; do
  rel="${f#*out/}"           # -> N<N>/p<p8>/data.dat
  dst="$OUTDIR/$rel"
  mkdir -p "$(dirname "$dst")"
  cp "$f" "$dst"
  n=$((n + 1))
done < <(find "$dl" -type f -path '*out/N*/p*/*/*.dat')

echo "=== merged $n cells into $OUTDIR ==="
if [ "$n" -eq 0 ]; then
  echo "warning: found no data.dat under $dl -- inspect it:" >&2
  echo "  simcloud -c $CLUSTER job list --batch $BATCH" >&2
else
  echo "next: uv run tools/analyze.py --all && uv run tools/heatmap.py --replot-all"
fi
