# Running brane on Simcloud (Apple ACS batch compute)

Notes and tooling for running the Fourier Monte Carlo engine on Apple's
Simcloud. Complements the workflow overview in [`README.md`](README.md).
**TL;DR: use `mr2-as` (M2 Ultra), not `mr2` (x86).**

## Scripts

| Script | Runs where | Purpose |
|---|---|---|
| [`overnight.sh`](overnight.sh) | laptop | **entrypoint** — sets production defaults, then execs `simcloud_submit.sh` |
| [`simcloud_submit.sh`](simcloud_submit.sh) | laptop | package repo + submit the whole grid as a batch |
| [`simcloud_task.sh`](simcloud_task.sh) | in-job | one batch task = one `(N,p8)` cell (maps `$SC_BATCH_ID`) |
| [`simcloud_monitor.sh`](simcloud_monitor.sh) | laptop | live progress-bar monitor of a batch |
| [`simcloud_fetch.sh`](simcloud_fetch.sh) | laptop | wait for batch, download output bundles, rebuild `data/` tree |

## Quick start

```sh
# 1. Submit a small grid slice first to validate the round-trip.
CLUSTER=mr2-as OWNER=hwt:atg:sph:$USER NET=e57cff0a-d781-4250-8ca5-065e283c8da1 \
  TOOLCHAIN=0 CPUS=16 NS=32,40,48 P8S=0.4,0.7 bash cloud/simcloud_submit.sh
#    -> prints the exact monitor/fetch commands with the batch id filled in

# 2. Watch it (live progress bar; uses cloud/.last_batch).
CLUSTER=mr2-as bash cloud/simcloud_monitor.sh

# 3. Pull results back into data/ once jobs finish.
CLUSTER=mr2-as bash cloud/simcloud_fetch.sh

# 4. Analyze locally.
uv run tools/analyze.py --all

# For the real production grid, just use the wrapper (edit knobs in its header):
bash cloud/overnight.sh
```

## How it maps onto Simcloud

The engine is **replica-parallel**: `nt` independent Markov chains run as OpenMP
threads, each doing `therm + sweeps` sweeps; results are averaged. Two facts
drive the design:

- **Statistics** scale as `1/sqrt(nt * sweeps)` — embarrassingly parallel, so
  more cores = more replicas = smaller error bars.
- **Reach** (large `N`) is bounded by *single-core clock*: one sweep costs
  `~N^4` and is **sequential**. Cores don't help a single big-`N` chain.

So each `(N, p8)` **cell is an independent Simcloud job** (`simcloud batch post`,
one job per cell, indexed by `$SC_BATCH_ID`), and within a cell `--cpus` = the
number of replicas. All cells run in parallel across the fleet, so grid wall
time ≈ the *slowest single cell*, not the sum.

## Cluster comparison (measured 2026-09-04)

Benchmark: `N=32`, `therm=20 sweeps=100` (120 "sweep-units"), single chain.

| Machine | single-thread wall | vs x86 | parallel efficiency |
|---|---|---|---|
| M4 P-core (laptop) | 5.9 s | 5.6× | ~80% @ 12 cores |
| **M2 Ultra (`mr2-as`)** | **9.0 s** | **3.6×** | **96–98% @ 8 cores** |
| x86 server (`mr2`) | 32.7 s | 1.0× | 73% @ 32 cores |

**The x86 fleet is ~5.6× slower single-thread than the M4 laptop** — a big
surprise, and fatal for large-`N` reach (a single N=120 chain would take days).
**M2 Ultra is only ~1.5× slower than the M4 and scales near-linearly**, making
it the right cluster for both goals.

### Projected per-cell wall (for a `therm=100 sweeps=4000` = 4100-sweep run)

| N | x86 `mr2` (nt=32) | M2 Ultra `mr2-as` (nt=16) |
|---|---|---|
| 32 | 26 min | **5 min** |
| 48 | 2.2 h | **26 min** |
| 64 | 6.8 h | **82 min** |
| 80 | 16.7 h | **3.3 h** |
| 96 | 34.6 h | **6.9 h** |
| 120 | 84.4 h | **16.9 h** |

`N^4` scaling was confirmed locally (N=48/N=32 wall ratio 5.25 ≈ 1.5⁴=5.06).
On x86, cap the grid at **N ≤ 64** and run big-`N` cells locally. On M2 Ultra the
whole grid is feasible: production now runs **N up to 180–200** at the shorter
`therm=300 sweeps=800` (=1100-sweep) setting — the per-cell wall then scales down
by ~4100/1100 from the table above, giving N=180 ~28.5 h and N=200 ~43.5 h (see
[`overnight.sh`](overnight.sh) for the current sizing). `TIMEOUT=72h` leaves
margin over the slowest cell; the `hwt:atg:sph` quota max timeout is 7 d.

### Intra-chain parallelism for large N (`it=`) — macOS only, avoid on Linux

By default cores = independent **replicas** (`nt=`), which buys *statistics*,
not per-cell speed. The engine also has an **inner-thread** knob `it=` that
parallelizes the `O(L²)` Metropolis step loop ([membrane.c](../src/membrane.c)) —
this is the *legacy* strategy (`legacy/simulate.c` ran one chain, inner loop
parallel, no replicas). Run `nt=1 it=<cores>` for one chain, or hybrid `nt*it=cores`.

**It helps only on macOS (LLVM libomp); it *regresses* on Linux (GCC libgomp)** —
which is what every Simcloud cluster runs. Measured:

| platform / N | it=1 | it=8 | it=16 |
|---|---|---|---|
| M4 laptop, macOS/libomp, N=120 | 64.9 s | **30.0 s** (2.2×) | — |
| M4 laptop, macOS/libomp, N=96 | 26.8 s | **15.6 s** (1.7×) | — |
| M2 Ultra, Linux/libgomp, N=96 | 29.4 s | 54.9 s ✗ | — |
| x86, Linux/libgomp, N=72 | 16.5 s | 20.8 s ✗ | 48.5 s ✗ |

libgomp's per-step fork/join cost grows with thread count and swamps the gain
(the step loop is re-entered `~L²` times per sweep). So:

- **On Simcloud (Linux): always use replicas, `it=1`** (the default). Don't set `IT>1`.
- `it=` is useful only for large-N runs **on a Mac** (`make` there links libomp).
- A persistent-team rewrite (fork once per sweep-run, `#pragma omp for` inside)
  would remove the per-step fork/join and likely make it pay on libgomp too —
  future work. The knob is default-off and gated (`LL≥20000`), so it never
  affects the replica path.

Correctness is unaffected — `it=1` vs `it=8` output is bitwise identical (RNG
stays serial; only the reduction sum reorders).

## Cluster-specific setup

### `mr2` — x86 (Sparks, NV)
- Publicly available; **direct internet** (apt works normally).
- Node limit **72 cpu**; requesting 64 may be rejected if no node has it free —
  32 is reliable.
- Needs a **toolchain bundle** (`TOOLCHAIN=1`, default): `bundle package --apt
  build-essential`, built once and reused by tag.
- Quota: personal `d_saykin` = 20 cpu; group `hw:others` = 5000 cpu (use
  `OWNER=hw:others:$USER` for the full grid).

### `mr2-as` — Apple Silicon M2 Ultra (Sparks, NV) — **recommended**
- **M2 Ultra**: max **16 cpu**, 168 GB, 7 d timeout. (M5 nodes also exist: only
  4 cpu but newest single-core; select by resource request.)
- **gcc preinstalled** in the ubuntu SMIs → set `TOOLCHAIN=0`, no bundle build.
- **No internet**, but apt reaches Apple-internal mirrors (works in-job).
- Requires **group quota** (personal quotas don't work) and a **Denali VPC**
  network id. Ours: `hwt:atg:sph` → closest ancestor `hwt` → MR2-AS net id
  **`e57cff0a-d781-4250-8ca5-065e283c8da1`**.
  (Discover with `simcloud -c mr2-as vpc network ls --owner hwt:atg:sph`, needs DC VPN.)
- **DC VPN**: docs say it's needed to read console / download output for M2
  Ultra jobs. In practice console + output worked without it for us — if a
  future run can't fetch results, request "Data Center VPN" in accessmanager.

VPC network ids for other groups/DCs: see the table at
`https://cloud.apple.com/docs/simcloud/tutorials/m2ultra/`.

## Gotchas learned the hard way

- **`simcloud job run` orphans jobs on interrupt.** If you Ctrl-C the local
  stream, the *remote* job keeps burning cores until its timeout. Prefer
  `job post` (as `simcloud_submit.sh` does) and monitor separately with
  `simcloud_monitor.sh`; always check `simcloud job list --status inprogress`
  after interrupts.
- **The repo's stale `brane` binary ships in the bundle** and is a macOS Mach-O.
  `make` thinks it's up to date and skips rebuilding → `exit 126` (can't exec)
  on Linux. Fix: build with **`make -B`** (force) in-job.
- **`bundle create --exclude` takes ONE regex**, matched against the *absolute*
  path (which contains the repo dir name `brane`). Multiple `--exclude` flags:
  only the last wins. Excluding a literal `brane` would nuke everything. Anchor
  on data dirs: `(^|/)(data|plots|...)(/|$)`. (Cut the source bundle 125 MB → 3 MB.)
- **No `--env` flag** on `job/batch post`. Pass params by prefixing the command
  string (safe when values have no spaces).
- **Batch index is `$SC_BATCH_ID`** (not `$SC_BATCH_INDEX`, despite some docs).
- **`error` status with `ExitCode: 0`** just means the job hit its `--timeout`
  *cap* after the task finished — harmless; give real runs generous timeouts.
- **`$USER` ≠ simcloud username** (`saykind` vs `d_saykin`). Derive it from
  `simcloud user info` when building an `--owner` fqn.
