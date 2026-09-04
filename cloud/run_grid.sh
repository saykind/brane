#!/usr/bin/env bash
# run_grid.sh -- portable "build + grind the whole (N, p8) grid" driver.
#
# Works on any Linux box (or macOS). Installs a C toolchain if missing, builds
# the engine, and runs every (N, p8) cell with ONE replica per core (nt=nproc),
# which is exactly the knob that fixes our statistics bottleneck: error on each
# G(q) point falls ~ 1/sqrt(nt * sweeps), so a 96-core box gives ~sqrt(96/12)
# ~ 2.8x smaller error bars than the 12-core laptop at the same wall-time.
#
# Everything is configurable via environment variables (defaults below). Results
# land in data/N<N>/p<p8>/data.dat -- the exact layout tools/analyze.py expects.
#
# Typical use on a fresh cloud VM (code already rsync'd or cloned into ./):
#   NT=96 SWEEPS=4000 EPS=0.005 bash cloud/run_grid.sh
#
set -euo pipefail

NS="${NS:-32,40,48,56,64,80,96,120}"          # lattice sizes (L=2N+1)
P8S="${P8S:-0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}"  # couplings (q8~p8)
NT="${NT:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc)}"  # replicas = cores
THERM="${THERM:-100}"                          # thermalization sweeps
SWEEPS="${SWEEPS:-4000}"                        # measurement-sweep cap per cell
EPS="${EPS:-0.005}"                             # Delta2 rel-error convergence target
MINSW="${MINSW:-100}"                           # floor before convergence can trip
OUTDIR="${OUTDIR:-data}"
CC_BIN="${CC:-gcc}"                             # GCC has native -fopenmp on Linux
S3_BUCKET="${S3_BUCKET:-}"                       # optional: aws s3 upload target
GCS_BUCKET="${GCS_BUCKET:-}"                      # optional: gsutil upload target
SHUTDOWN="${SHUTDOWN:-0}"                         # 1 -> poweroff when done (spot)

echo "=== brane grid run: NT=$NT NS=$NS P8S=$P8S SWEEPS=$SWEEPS EPS=$EPS ==="

# --- 1. toolchain (best-effort; skip if already present) --------------------
if ! command -v "$CC_BIN" >/dev/null 2>&1; then
  if   command -v apt-get >/dev/null 2>&1; then sudo apt-get update -y && sudo apt-get install -y build-essential git; fi
  if   command -v dnf     >/dev/null 2>&1; then sudo dnf groupinstall -y "Development Tools" && sudo dnf install -y git; fi
  if   command -v yum     >/dev/null 2>&1; then sudo yum groupinstall -y "Development Tools" && sudo yum install -y git; fi
fi

# --- 2. build ---------------------------------------------------------------
make CC="$CC_BIN"
./brane -h >/dev/null 2>&1 || true

# --- 3. run every cell ------------------------------------------------------
IFS=',' read -ra NARR <<< "$NS"
IFS=',' read -ra PARR <<< "$P8S"
t0=$SECONDS
for N in "${NARR[@]}"; do
  for P in "${PARR[@]}"; do
    pdir=$(printf 'p%.2f' "$P")
    out="$OUTDIR/N$N/$pdir/data.dat"
    mkdir -p "$(dirname "$out")"
    echo ">>> N=$N p8=$P  ($(($SECONDS - t0))s elapsed)"
    ./brane "N=$N" "p8=$P" "nt=$NT" "therm=$THERM" "sweeps=$SWEEPS" \
            "eps=$EPS" "minsweeps=$MINSW" "out=$out"
  done
done
echo "=== grid done in $(($SECONDS - t0))s ==="

# --- 4. package + optional upload ------------------------------------------
tar czf "brane_data_$(hostname -s).tgz" "$OUTDIR"
echo "packaged -> brane_data_$(hostname -s).tgz"
[ -n "$S3_BUCKET" ]  && aws s3 cp   "brane_data_$(hostname -s).tgz" "s3://$S3_BUCKET/"  && echo "uploaded to s3://$S3_BUCKET/"
[ -n "$GCS_BUCKET" ] && gsutil cp   "brane_data_$(hostname -s).tgz" "gs://$GCS_BUCKET/" && echo "uploaded to gs://$GCS_BUCKET/"

# --- 5. optional self-terminate (spot instances) ---------------------------
[ "$SHUTDOWN" = "1" ] && sudo poweroff
