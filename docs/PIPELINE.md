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
| `overrelax` | over-relaxation sweeps interleaved per Metropolis sweep (`0`=off) — decorrelation, see below |
| `outdir` | base dir; engine builds the descriptive path (below) |
| `out` | explicit output path (overrides `outdir`) |
| `series` | write per-sweep instantaneous Δ₂ (replica 0) for τ measurement |

Defaults: `therm=300` (matches legacy), `eps=0.005`, `minsweeps=200`.

**Robustness:** output is written atomically (temp+rename) and **checkpointed
every 60 s**, so a killed run keeps its latest data. A per-block convergence
trace is written to `<out>.trace` (sweeps, Δ₂, rel_err, **accept**, wall_s) and
a per-mode Metropolis acceptance map to `<out>.accept` (`q1 q2 qmag proposed
accepted rate`).

### Step size & acceptance (Tröster OFMC)

The per-mode trial step is momentum-dependent,
`d[k] = d0/q² · (1 + Y/q²)^{−0.13}` (in `geometry_make`). Since the harmonic
amplitude is `⟨|h_q|⟩ ~ q⁻²`, the `d0/q²` factor is exactly Tröster's
`d[k] ~ ⟨|h_k|⟩` heuristic (2013, §4: tune the step per wave vector so each mode
has ~50% acceptance → uniform τ, no critical slowing down); the `(1+Y/q²)^{−0.13}`
softening bends `d` toward the anomalous amplitude `q^{−(4−η)/2}` below the
crossover.
- *Measured (N=64, p8=0.4):* overall acceptance 49.9%, **stable across sweeps**
  (0.497–0.510 from the first block on — it does not drift as the chain
  thermalizes). Acceptance is **flat at 50% for q above the crossover**, but
  drops to ~30% for the lowest-q modes (`q ≲ p8`) — the step slightly *over*-steps
  there (our `q^{−1.74}` vs the true `q^{−1.65}`). So the momentum-dependent step
  is close to optimal OFMC but **not fully tuned at the slowest modes** — the
  natural next step is a warm-up phase that adapts `d[k]` per mode to hit 50%
  (true OFMC). Tuning the *step* to 50% acceptance is a legitimate MC-efficiency
  criterion (detailed balance holds for any step; only efficiency changes) — not
  to be confused with tuning an *estimator* to a target η.

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
| `tools/plot_acceptance.py` | acceptance vs sweep + vs \|q\| from a run's `.trace`/`.accept` |
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

## Decorrelation & warm start (implemented)

Two levers were added to the engine and validated for correctness; both are
**large-N tools** whose benefit is marginal in the small-N regime tested so far.

**Over-relaxation** (`overrelax=R`). Holding all other modes fixed, a single
mode's conditional energy is exactly quadratic apart from a tiny quartic
self-coupling (the `q=±2k` terms). The move builds that quadratic (gradient
`g`, Hessian `M`) in one O(L²) pass and reflects the mode about the conditional
minimum, `u = −2M⁻¹g` (a volume-preserving involution), then accepts with the
*exact* ΔE — which corrects the quartic residual and keeps detailed balance
exact. `R` such sweeps are interleaved per Metropolis sweep (OR alone is
near-microcanonical, hence non-ergodic). Reuses the validated incremental-S
kernel (`mode_dS_energy`).
- *Correctness:* the S-consistency + reality test covers the OR path
  (`make test`), and a cross-replica comparison (nt=8) found R=1 vs R=0
  statistically identical — median |z|=0.66, `|z|>2` in 5.1% of modes, `|z|>3`
  in 0.6% (Gaussian expects 4.6%, 0.3%). OR reflects with ~100% acceptance.
- *Benefit (N=36, p8=0.4):* τ_int(Δ₂) drops 1.25→0.71→0.66 sweeps for R=0,1,2,
  and low-q error bars shrink markedly. **But** cost is ~2.8×/5.2× per sweep, so
  net efficiency (independent samples per wall-second) is ×0.62/×0.36 vs R=0 —
  net-negative here because τ≈1 already. OR is expected to pay off only in the
  critical-slowing-down regime (large N, large τ); **not yet demonstrated to win
  on wall-clock** — needs a large-N test.

**Anomalous warm start** — *tried and removed.* Seeding `|h_q|² ~ q^{−(4−η)}`
below `q_c≈p8` instead of harmonic `q⁻⁴` began closer to equilibrium (ensemble
Δ₂ excess +48% vs +158% at N=64) and shaved a few sweeps off thermalization, but
the transient is only ~10–15 sweeps for either start (≪ the conservative
`therm=300`), so it saved ~10 sweeps at most. Not worth the extra code path;
reverted to cold (harmonic) start only.

## Open items / TODO

- **Large-N decorrelation test** — measure whether over-relaxation's τ reduction
  beats its ~3× per-sweep cost at large N (where τ is large). Needs a cloud run.
- **Refine `plateau_lowq`** — the current heuristic can be pulled by a low-q
  finite-size shelf or the high-q crossover approach. Needs a principled,
  crossover-aware estimator. Do **not** tune it to reproduce an expected η.
- **Proper error-vs-samples study** — many independent seeds, per-observable
  spread, to actually establish the scaling law and finalize eps/therm.
