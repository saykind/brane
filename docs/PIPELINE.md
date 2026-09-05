# brane pipeline: data format, analysis, and workflow

This documents the modernized simulate→store→analyze pipeline on the
`modernize-fmc` branch. Physics background is in [model.md](model.md);
running on Apple Simcloud is in [../cloud/SIMCLOUD.md](../cloud/SIMCLOUD.md).

## Engine (`./brane`)

Replica-parallel Fourier Monte Carlo. Key options (`./brane -h`):

| option | meaning |
|---|---|
| `N`, `p8`, `n` | half lattice size (L=2N+1), coupling, move-zone |
| `nt` | replicas (independent chains, averaged) — statistics |
| `it` | inner threads per replica (intra-chain) — large-N latency; **helps only on macOS/libomp, regresses on Linux**, keep `it=1` on the cloud |
| `therm`, `sweeps` | thermalization + max measurement sweeps |
| `eps` | adaptive stop on rel. error of Δ₂=Σ_q⟨|h_q|²⟩ (`eps=0` → fixed length) |
| `minsweeps`, `block`, `meas`, `d0`, `seed` | convergence floor, block size, measure stride, step, RNG seed |
| `outdir` | base dir; engine builds the descriptive path (below) |
| `out` | explicit output path (overrides `outdir`) |
| `series` | write per-sweep instantaneous Δ₂ (replica 0) for τ measurement |

Defaults: `therm=300` (matches legacy), `eps=0.005`, `minsweeps=200`.

**Robustness:** output is written atomically (temp+rename) and **checkpointed
every 60 s**, so a killed run keeps its latest data. A per-block convergence
trace is written to `<out>.trace` (sweeps, Δ₂, rel_err, wall_s).

## Output layout & format

Descriptive paths keep different configs from overwriting each other:

```
data/N<N>/p<p8>/<stop>/therm<T>_nt<NT>_it<IT>_seed<S>.dat
                 └ stop = eps<eps> (adaptive) | fixed<sweeps> (fixed length)
```

Each `.dat` has a `key=value` metadata header (parsed by `analyze.py`) then
columns `q1 q2 qx qy qmag G Gerr Ginv`:

```
# N=100 L=201 n=100 p8=0.4000 N8=12 Y=... d0=2.60 seed=12345
# nt=16 it=1 cores=16
# therm=300 sweeps=800 sweeps_cap=800 min_sweeps=200 block=20 meas_every=1 steps_per_sweep=...
# eps=0.005 rel_err=0.0034 converged=0
# samples=12800 accept_rate=0.499 wall_s=... nu=0.048 nu_err=0.018
# engine_sha=<git>
```

`example_data/` (legacy large-N reference, N=100–200, p8=0.3) has been
reformatted to this same format via `tools/reformat_legacy.py`; unknown legacy
params are `NA` and `Gerr=0` (single chain → no per-mode error). The untouched
originals are in `legacy/example_data/` (gitignored) and in git history.

## Analysis (`tools/analyze.py`) — no angular averaging

Every mode `(q_x,q_y)` is a point `(q_r=|q|, G^-1)`; we do **not** average G over
angle (that would hide the membrane's anisotropy). Reported:

- **windowed η** — unbinned weighted log-log fit over the q-window; `η=4−slope`.
  Biased low when the window reaches the crossover `q8~p8`.
- **running η_eff(q_r)** — sliding local log-log slope over the raw cloud.
- **low-q plateau** — the flat part of `η_eff` at low q *is* the physical q→0
  exponent. `plateau_lowq()` is a **first, imperfect heuristic** (widest
  flat window) — see "Open items".

Plots color every mode by polar angle (anisotropy visible). `--combined`
pools all modes at fixed p8 across N into one fit (2-panel: cloud + running
exponent with plateau). `--all [GLOB]` batch-plots to `plots/<mirror>/`.

## Tools

| tool | purpose |
|---|---|
| `tools/analyze.py` | per-file & combined η analysis + plots |
| `tools/autocorr.py` | integrated autocorrelation time τ from a `series=` file |
| `tools/study_convergence.py` | error/thermalization vs sweeps from a `.trace` or log |
| `tools/reformat_legacy.py` | legacy dump → modern format |
| `tools/heatmap.py`, `green_map.py`, `scaling.py` | grids / maps / finite-size |

## Cloud

`cloud/` runs the grid on Apple Simcloud (M2 Ultra recommended). See
`cloud/SIMCLOUD.md`. `cloud/overnight.sh` launches a production grid; results
come back via `cloud/simcloud_fetch.sh` (pulls output bundles through the
bundle service, then merges into `data/`).

## Convergence findings (preliminary)

From a single N=100 study (`tools/study_convergence.py`):
- Error early-decay is consistent with 1/√(samples), but the tail is dominated
  by the noisy 4-replica error *estimator* — **not confirmed**. A proper
  multi-seed study (spread of the observable across independent seeds) is needed.
- Thermalization (running-mean Δ₂) settles by ~150 sweeps → `therm=300` is safe.
- `eps=0.005` reaches in ~500–800 sweeps at nt=16 — a reasonable production goal.

## Open items / TODO

- **Refine `plateau_lowq`** — the current heuristic can be pulled by a low-q
  finite-size shelf or the high-q crossover approach. Needs a principled,
  crossover-aware estimator. Do **not** tune it to reproduce an expected η.
- **Over-relaxation** — microcanonical moves to cut the autocorrelation time τ
  (biggest expected decorrelation win). Measure τ before/after with `autocorr.py`.
- **Anomalous warm start** — seed `|h_q|² ~ q^{−(4−η)}` (η≈0.7) instead of the
  harmonic `q^{−4}` to shorten thermalization.
- **Proper error-vs-samples study** — many independent seeds, per-observable
  spread, to actually establish the scaling law and finalize eps/therm.
