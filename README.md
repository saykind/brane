# brane

**Fourier Monte Carlo simulation of a 2D crystalline membrane (graphene) in
the flat phase.**

`brane` measures the anomalous elasticity exponent **η** and the **Poisson
ratio ν** of a thermally fluctuating 2D membrane by Metropolis sampling of the
height field's Fourier modes `h_q`. It is a modernized, multi-threaded
rewrite of the code behind the M.Sc. thesis (`η = 0.78 ± 0.02`,
`ν = −0.76 ± 0.05`).

For the physics, the algorithm, full literature sources, and the acceleration
roadmap, see **[docs/model.md](docs/model.md)**.

---

## Quick start

```bash
# 1. One-time dependency (OpenMP runtime for Apple clang):
brew install libomp

# 2. Build:
make

# 3. Run the default simulation (N=36, p8=0.4, ~52 s on 16 cores):
./brane -v

# 4. Set up the Python analysis env (once) and extract eta + plot:
uv sync
uv run tools/analyze.py data/N=36.dat
```

Expected default output: acceptance ≈ 50 %, and `analyze.py` reporting
`eta ≈ 0.75` plus a PNG of the inverse Green function.

### Dependencies

| Tool | Needed for | Install |
|------|-----------|---------|
| `libomp` | multithreading (Apple clang) | `brew install libomp` |
| `uv` + `pyproject.toml` | Python analysis env (numpy + matplotlib) | `uv sync` (uv via `brew install uv`) |
| `gnuplot` | plot fallback if matplotlib is absent | optional, already installed |

Without `uv` you can still fit η with the system `python3`+`numpy`
(`python3 tools/analyze.py …`); plots then fall back to gnuplot.
On Linux with GCC no `libomp` is needed: `make CC=gcc`.

---

## Usage

```
./brane [options]

  N=<int>        half lattice size, L = 2N+1        (default 36)
  n=<int>        half move-zone size, l = 2n+1      (default N)
  p8=<float>     interaction strength, 0 < p8 < pi  (default 0.4)
  nt=<int>       independent replicas / threads     (default 12)
  therm=<int>    thermalization sweeps per replica  (default 80)
  sweeps=<int>   measurement sweeps per replica     (default 80)
  meas=<int>     measure every M sweeps             (default 1)
  d0=<float>     base Metropolis step size          (default 2.6)
  seed=<int>     base RNG seed (reproducible)       (default 12345)
  out=<path>     output file (default data/N=<N>.dat)
  -v             per-replica progress
  -h             help
```

One **sweep** = `l·l` attempted single-mode updates. Total statistics =
`nt × sweeps` measurements. Each replica is an independent Markov chain seeded
by `(seed, replica_index)`, so runs are fully reproducible.

### Reaching the anomalous regime on a small lattice

`η` lives in the window `3·a ≲ q ≲ p8`. On a small lattice this window is
narrow; widen it by raising `p8` (physical range up to ~π):

```bash
./brane N=36 p8=0.8 nt=16 therm=80 sweeps=80   # η ≈ 0.74 in ~50 s
```

### Multi-size sweep (finite-size study, thesis-style)

```bash
./run.sh            # N = 20…36, each well under a minute
```

---

## Output format

`data/N=<N>.dat` is a text table (one row per non-zero mode):

```
# Fourier MC membrane   N=36 L=73 p8=0.4000 samples=1280 nu=0.036500
# q1 q2 qx qy qmag G Ginv
q1  q2  qx  qy  |q|  G(q)=<|h_q|^2>  1/G(q)
```

The header carries `N, L, p8, samples, nu`; `tools/analyze.py` parses it
automatically.

---

## Testing

```bash
make test
```
`tests/test_correctness.c` runs the chain, then recomputes `S_q` from scratch
and checks it against the incrementally-maintained array (agreement `< 4e-14`)
plus the reality condition `h_{−q} = conj(h_q)`.

## Benchmark

```bash
./tools/bench.sh          # legacy vs new, single-mode updates/sec
```

Representative result at `N=40` on the 16-core M4 Max (raw update throughput,
directly comparable across parallelization strategies):

| run | time (s) | chains | updates/sec |
|-----|---------:|-------:|------------:|
| legacy nt=1  | 3.60 | 1  | 18,225 |
| legacy nt=16 | 9.60 | 1  | 6,834  |
| new nt=1     | 1.75 | 1  | 37,491 |
| new nt=16    | 5.44 | 16 | **192,971** |

- **Legacy gets *slower* with more threads** (nt=16 is ~0.4× nt=1): the
  original opens an OpenMP `parallel for` inside the per-move routine, ~L²
  times per sweep, so fork/join overhead dwarfs the tiny inner loop.
- **New vs legacy at the *same* 16 threads: ~28–30×.**
- **New vs legacy serial (1 thread each): ~2×** (contiguous arrays, PCG32,
  precomputed tables).
- **New's own scaling from 1→16 threads: ~5.4×** (~34 % of the ideal 16×).
  This M4 Max has 12 performance + 4 efficiency cores and each replica is
  memory-bandwidth heavy, so replica scaling is sublinear but still large.

The legacy code is built (with a minimal OpenMP fix) by `legacy/build.sh`;
the original only failed to compile because Apple clang needs libomp flags.

### Core scaling

```bash
uv run tools/scaling.py --N 40      # table + data/scaling.png
```

Replica parallelism at N=40 on the M4 Max (12 P + 4 E cores): linear to ~2
cores, then memory-bandwidth limited (each replica streams its own ~0.4 MB of
arrays), plateauing around **12 threads**; the 4 E-cores add nothing and, due
to the end-of-run barrier on the slowest replica, 16 threads is slightly worse
than 12. Practical sweet spot for this size: **~10-12 threads**.

| cores | speedup | efficiency |
|------:|--------:|-----------:|
| 2  | 2.00× | 100% |
| 4  | 3.62× | 91% |
| 8  | 4.79× | 60% |
| 12 | 5.45× | 45% |
| 16 | 5.08× | 32% |

---

## Project layout

```
src/
  membrane.h        API, Config/Geometry/Replica/Result structs
  membrane.c        core engine: geometry, incremental Metropolis, observables
  main.c            CLI + replica-parallel driver + I/O
  pcg.h             per-stream PCG32 RNG (reproducible, thread-safe)
tests/
  test_correctness.c  incremental-vs-exact S_q validation
tools/
  analyze.py        eta fit (numpy) + plot (matplotlib/gnuplot)
  bench.sh          legacy-vs-new throughput benchmark
  scaling.py        core-scaling benchmark (table + plot)
docs/
  model.md          physics, algorithm, sources, acceleration roadmap
legacy/             original thesis code (build.sh fixed for macOS libomp)
example_data/       reference .dat files from the thesis runs (old format)
lib/                method papers (Tröster; Los et al.)
Makefile            OpenMP autodetection (macOS libomp / Linux gcc)
pyproject.toml      Python analysis env (uv sync)
run.sh              multi-size sweep
```

## Notes

- **η** is reliable at these sizes; **ν** is genuinely hard and needs large
  `L` / long runs to stabilize (see docs/model.md §1.3, §3).
- The original code lives in `legacy/` and is preserved for provenance; it
  required `libomp` and used a single chain with an inner-loop `parallel for`.
