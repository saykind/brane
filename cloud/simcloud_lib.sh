#!/usr/bin/env bash
# simcloud_lib.sh -- shared helpers for the Simcloud driver scripts.
# Source this: . "$(dirname "$0")/simcloud_lib.sh"

# --- progress bar (renders to stderr, single updating line) -----------------
sc_bar() {
  local cur=$1 tot=$2 label="${3:-}" width="${SC_BAR_WIDTH:-30}"
  local filled=0
  [ "$tot" -gt 0 ] && filled=$(( cur * width / tot ))
  [ "$filled" -gt "$width" ] && filled=$width
  local pct=0; [ "$tot" -gt 0 ] && pct=$(( cur * 100 / tot ))
  local b e
  printf -v b '%*s' "$filled" ''; b=${b// /#}
  printf -v e '%*s' "$(( width - filled ))" ''; e=${e// /.}
  printf '\r[%s%s] %3d%% %d/%d %s\033[K' "$b" "$e" "$pct" "$cur" "$tot" "$label" >&2
}

# --- create the source bundle (unpacks at /brane), echo its id --------------
# args: <repo-dir>   env: none
sc_src_bundle() {
  local repo="$1"
  simcloud -q bundle create -C "$repo" . --root /brane \
    --tag "brane-src" \
    --exclude '(^|/)(data|plots|eta|example_data|legacy|\.git|\.venv|__pycache__)(/|$)'
}

# --- find-or-build the build-essential toolchain bundle for an SMI ----------
# args: <smi>
sc_toolchain() {
  local smi="$1" name="brane-toolchain-$1" id
  id=$(simcloud -q bundle find --tag "$name" 2>/dev/null || true)
  if [ -z "${id:-}" ]; then
    echo "--- building toolchain bundle (build-essential for $smi) ---" >&2
    id=$(simcloud -q bundle package --apt build-essential --smi "$smi" --tag "$name")
  fi
  echo "$id"
}

# --- poll a batch and draw a progress bar until every job is terminal -------
# args: <batch-id> <total> [poll-seconds]
sc_watch_batch() {
  local batch="$1" total="$2" poll="${3:-10}" done_n
  while :; do
    # count jobs in a terminal state (complete/error/cancelled/failed/timeout)
    done_n=$(simcloud -q job list --batch "$batch" -f '{{.Status}}' 2>/dev/null \
      | grep -ciE 'complete|error|cancel|fail|timeout' || true)
    done_n=${done_n:-0}
    sc_bar "$done_n" "$total" "jobs done"
    [ "$done_n" -ge "$total" ] && { echo >&2; break; }
    sleep "$poll"
  done
}
