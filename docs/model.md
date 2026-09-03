# The Model, the Algorithm, and How We Accelerate It

This document summarizes the physics simulated by `brane`, the Fourier Monte
Carlo (FMC) algorithm it uses, and the concrete plan for modernizing and
accelerating the code. Every non-obvious claim carries a source; full
references are collected at the end.

---

## 1. Physical model: a 2D crystalline membrane in the flat phase

`brane` simulates a **two-dimensional crystalline (tethered) membrane** — the
continuum model of graphene and similar atomically-thin crystals — in its
**flat phase** at finite temperature. The membrane is parametrized by
`r(x) = (x + u(x), h(x))`, with in-plane phonons `u` and out-of-plane height
`h` over a 2D base lattice `x` [thesis `chapters/model.tex`; Nelson, Piran &
Weinberg 2004, ch. 6].

The elastic free energy, after expanding in gradients and rescaling, is
[thesis `model.tex`, eq. `free_energy_r`]:

```
F[r] = 1/2 ∫ d²x [ κ (∇²h)²
                   + (μ/2)(∂_α r · ∂_β r − δ_αβ)²
                   + (λ/4)(∂_α r · ∂_α r − D)² ].
```

**Why it is interesting.** Thermal fluctuations couple the soft bending modes
(`κ q⁴`) to in-plane stretching. Integrating out the phonons leaves an
effective, *strongly interacting* theory for `h` in which the bending rigidity
is anomalously renormalized: `κ_eff(q) ∼ q^(−η)` as `q → 0`. This "anomalous
Hooke's law" makes membranes far stiffer against bending at long wavelengths
than the harmonic theory predicts, and gives a **negative Poisson ratio** in
the universal regime [thesis `motivation.tex`, `problem.tex`; Le Doussal &
Radzihovsky 2018 review].

The central critical exponent is **η** (the anomalous elasticity exponent),
with the accepted value `η ≈ 0.78–0.85`; and the **Poisson ratio `ν`** in the
zero-stress universal regime, reported around `−0.3` to `−0.76`.

### 1.1 Effective Fourier action actually simulated

Integrating out the phonons (Gaussian) leaves a height-only action that is
**local in Fourier space but long-ranged in real space** — the ideal setting
for FMC. With external stress set to `σ = 0`, in dimensionless form [thesis
`chapters/mcmc.tex`, eq. `free_energy_final`]:

```
F  = Σ_{q≠0} { 1/2 (2π q / L)⁴ |h_q|²  +  (2π/3)(p₈²/L²) |S_q|² }

S_q = Σ_{k≠0} [ (2π/L) k × q̂ ]² ( h_k · h_{k+q} )
```

- The first term is the **harmonic** bending energy (`∝ q⁴`).
- The second is the **anharmonic** coupling; `S_q` is a *convolution* of the
  height field with itself weighted by the transverse projector `[k × q̂]²`.
- **All momenta are replaced by their lattice sine values**,
  `(2π/L) k₁ ↦ sin(2π k₁/L)` and `q² ↦ 4 sin²(π q₁/L) + 4 sin²(π q₂/L)`, so
  the model lives on a periodic `L × L` Brillouin zone (`L = 2N+1`) with the
  reality condition `h_{−q} = conj(h_q)` [thesis `mcmc.tex`].

### 1.2 Dimensionless coupling `p₈`

The single control parameter is the dimensionless interaction strength
[thesis `overview.tex`]:

```
p₈ = p* a = a √( (3/16π) Y₀ T / κ² )
```

For graphene (`a = 2.46 Å`, `Y₀ ≈ 22 eV·Å⁻²`, `κ ≈ 1.1 eV`) at `T = 300 K`
this gives **`p₈ ≈ 0.41`** — which is why `brane`'s default is `p8=0.4`.
Larger `p₈` pushes the crossover scale `q₈ ∼ p₈` to larger `q`, widening the
anomalous scaling window that is accessible on a small lattice.

### 1.3 Observables

- **Green function** `G(q) = ⟨|h_q|²⟩`. In the harmonic theory
  `G⁻¹(q) ∼ q⁴`; anomalous elasticity gives `G⁻¹(q) ∼ q^(4−η)` for `q ≪ p₈`.
  A straight-line fit of `log G⁻¹` vs `log q` in that window has slope
  `(4 − η)`, so **`η = 4 − slope`** (`tools/analyze.py`).
- **Poisson ratio** [thesis `mcmc.tex`]:

  ```
  ν = − Σ_{|q|<q₈} [⟨Kx Ky⟩ − ⟨Kx⟩⟨Ky⟩] / Σ_{|q|<q₈} [⟨Kx Kx⟩ − ⟨Kx⟩²],
  Kα_q = sin²(2π kα/L) |h_q|².
  ```

The thesis reports `η = 0.78 ± 0.02` and `ν = −0.76 ± 0.05` using lattices up
to `L = 360` (N=180) [thesis `mcmc.tex`, `overview.tex`]. These are **large-N,
long-run** numbers; on a sub-minute `N ≈ 36` lattice `brane` already recovers
`η ≈ 0.74–0.75`, while a clean `ν` needs the larger lattices.

---

## 2. The Fourier Monte Carlo algorithm

### 2.1 The move and the cost problem

FMC treats the **Fourier amplitudes `h_q` themselves as the MC degrees of
freedom** (real space is never touched). A single move picks one mode `k₀`
and shifts it together with its conjugate partner to preserve reality
[Tröster 2008 CPC, eq. 8; thesis `mcmc.tex`]:

```
h_{k₀} → h_{k₀} + ε,   h_{−k₀} → h_{−k₀} + conj(ε),   |ε| < r_ε.
```

The harmonic energy is diagonal, so its change is `O(1)`. The trouble is the
anharmonic term `Σ_q |S_q|²`: recomputing every `S_q` from scratch is a
double lattice sum — **`O(L⁴)` per move**, hopeless.

### 2.2 The acceleration trick (the heart of the method)

Store `S_q` and update only its **increment** `dS_q` each move. Writing the
quartic term as `Σ_x (h²(x))²`, which is diagonal in the Fourier amplitudes of
the *squared field*, Tröster derives the exact increment for a shift `ε` at
`k₀` [**Tröster 2008 CPC, eq. 14** — the formula the 2013 papers cite but omit]:

```
dS_q = (1/√N)[ 2ε h_{q−k₀} + 2 conj(ε) h_{q+k₀}
              + ε²        δ_{q, 2k₀}
              + 2|ε|²     δ_{q, 0}
              + conj(ε)²  δ_{q, −2k₀} ].
```

The first two terms touch every `q`; the last three are single-`q`
corrections. The Metropolis energy change follows in `O(L²)`
[Tröster 2008 CPC, eqs. 11–13]:

```
ΔE_harm = 2 D̃(k₀) Re[ ε conj(h_{k₀}) + ... ]
ΔE_anh  = (Y/L²) Σ_q Re[ (2 S_q + dS_q) conj(dS_q) ] / 4 .
```

This drops the per-move cost from **`O(L⁴)` to `O(L²)`** — the entire reason
the simulation is feasible. On acceptance, `h_{±k₀}` and every `S_q += dS_q`
are committed (`O(L²)`). `brane`'s `src/membrane.c` is a faithful, contiguous
port of this scheme; the increment is validated term-by-term against a
brute-force convolution to `< 4e-14` in `tests/test_correctness.c`.

### 2.3 Step-size tuning (OFMC)

With a single step radius, the acceptance rate is ~100% for small `q` but
collapses for large `q`, causing critical slowing down. **Optimized FMC**
treats each `k₀` as its own move type with a per-wavevector radius `r_ε(k₀)`
tuned to ≈ 50% acceptance; detailed balance holds because each `ε` is drawn
symmetrically [Tröster et al. 2013, PRB 87:144101 & EPL]. `brane` uses the
heuristic momentum-dependent step `d(k) ∼ 1/q²` softened by the coupling
(from the reference code), which already yields a flat **≈ 50% acceptance**
across all modes — matching the OFMC target out of the box.

### 2.4 What the method is good for

FMC pays off precisely when the Hamiltonian is **simple in Fourier space but
nonlocal in real space**, and one wants **long-wavelength/critical** behavior
(restrict the active modes to a cutoff and "stretch" the accessible `L`). The
requirement is that the field enter *at most quadratically* after introducing
the squared field [Tröster 2008 PRL 100:140602]. The same machinery was used
for the compressible φ⁴ model (Fisher renormalization) and, with Ewald
tabulation of the interaction, the long-range antiferromagnetic Ising model
[Tröster 2008 PRL; Tröster 2010 PRB 81:012406].

---

## 3. Numerical context (what others found)

Selected reported values [thesis `overview.tex`, Table 1; and the papers in
`lib/`]:

| Year | Author | Method | p₈ | L² | η | ν(σ=0) |
|---|---|---|---|---|---|---|
| 2009 | Los et al. | Molecular dynamics | – | 190² | ≈0.85 | |
| 2013 | Tröster | Fourier MC | 0.44 | 640²(40²) | 0.795 ± 0.01 | |
| 2016 | Los et al. | Real-space MC | 0.4 | 195² | ≈0.84 | +0.275 |
| 2019 | Saykin (thesis) | Fourier Green | 0.3 | 360² | 0.78 ± 0.02 | −0.76 ± 0.05 |

The spread across methods (`η ≈ 0.78–0.85`, `ν` from `−0.76` to `+0.28`)
reflects genuine finite-size and methodological difficulty — especially for
the Poisson ratio.

---

## 4. How we modernize and accelerate the code

### 4.1 What was wrong with the original

- **Did not build** on this machine: Apple `clang` needs Homebrew `libomp`
  for `-fopenmp` (now installed and auto-detected in the `Makefile`).
- **OpenMP anti-pattern:** a `parallel for` was opened *inside* the per-move
  routine, which runs `L²` times per sweep over a tiny `L²` loop — fork/join
  overhead dominated and made it *slower* than serial at small `L`. The one
  form of parallelism that MC actually offers cheaply (independent chains) was
  unused.
- Array-of-pointers `double**` (pointer chasing, no vectorization),
  pervasive `%L` in hot loops, global non-reproducible `rand()`, a `data/` vs
  `example_data/` path mismatch, fragile gnuplot post-processing, no tests.

### 4.2 The idea we implemented: **replica parallelism**

A single Metropolis chain is inherently sequential (each step depends on the
previous accept/reject), so parallelizing *within* a chain is where the
original went wrong. Equilibrium averages, however, are identical across
independent chains. So `brane` runs **one independent Markov chain per core**
(16 on this M4 Max), each with its own PCG32 stream, and **averages the Green
function and Poisson correlators across chains** at the end. This

- scales almost linearly with cores (no per-step synchronization),
- is fully reproducible from `(seed, stream)`,
- and gives `16 ×` the statistics for the same wall-clock.

This is *not* described in the Tröster papers (they mention only "a parallel
version... described elsewhere") — for equilibrium sampling it is the natural,
robust choice. Combined with **contiguous flat arrays**, a **precomputed
modulo/step/reciprocal table**, and a fast per-stream RNG, `brane` runs
`N=36, p₈=0.4` in ~52 s on 16 cores and recovers `η ≈ 0.75`.

### 4.3 A genuinely new idea worth trying next (documented, not yet built)

None of the papers in `lib/` use an FFT, and all state the per-move cost as
`O(L²)`. But `S_q` is a set of three cross-correlations:

```
S_q = q̂₂² A_q − 2 q̂₁q̂₂ B_q + q̂₁² C_q,
A_q = Σ_k (k₁² h_k) conj(h_{k+q}), etc.
```

Each `A, B, C` is a convolution, so the **full `S_q` array can be built in
`O(L² log L)` with three FFTs** instead of `O(L⁴)`. Single-mode Metropolis
does not benefit (the incremental update is already `O(L²)` per move), but
this unlocks two things the literature has not exploited here:

1. **Global moves via Hybrid/Hamiltonian Monte Carlo (HMC) or
   Fourier-accelerated Langevin.** The force `∂F/∂h_q` needs `S_q` and its
   derivative, both `O(L² log L)` by FFT. One HMC trajectory then updates
   *all* `L²` modes at once for `O(L² log L)` per step and decorrelates the
   whole lattice — potentially a large reduction in autocorrelation time
   versus single-mode sweeps.
2. **`O(L² log L)` re-thermalization / consistency checks**, replacing the
   `O(L⁴)` `calcS`.

This is the recommended path to *much* larger `L` (the regime where `η` and
especially `ν` sharpen). It is a research-grade reimplementation and is left
as future work; the current release deliberately preserves the validated
single-mode scheme so results match the thesis 1:1.

---

## 5. Sources

**Thesis chapters** (`/Users/saykind/GitHub/thesis_msc/chapters/`):
`overview.tex`, `mcmc.tex`, `model.tex`, `problem.tex`, `motivation.tex`,
`numerics.tex`, and reference code `thesis_msc/code/{s.c, simulate.c}`.

**Method papers** (`/Users/saykind/GitHub/brane/lib/`):

- A. Tröster, *"Monte Carlo simulation in Fourier space,"* Computer Physics
  Communications **179** (2008) 30–33. — **Primary algorithm source**;
  incremental update eq. 14, energy differences eqs. 11–13, `O(N)` per move.
- A. Tröster, *"Fourier Monte Carlo Simulation of Crystalline Membranes in the
  Flat Phase"* (2013). — Membrane action, discretization, OFMC step tuning,
  observables, finite-size scaling.
- A. Tröster, *"High-precision Fourier Monte Carlo simulation of crystalline
  membranes"* (2013). — High-precision `η`, error analysis (jackknife,
  autocorrelation, blocking), 50% ± 5% acceptance target.
- A. Tröster, *"Evidence for Fisher Renormalization in the Compressible φ⁴
  Model,"* Phys. Rev. Lett. **100**, 140602 (2008). — FMC applied to
  compressible φ⁴; requirement that the field be at most quadratic after
  squared-field substitution.
- A. Tröster, *"Short-range character of the antiferromagnetic Ising model
  with 1/rᵖ interaction,"* Phys. Rev. B **81**, 012406 (2010). — FMC with
  Ewald tabulation of a long-range interaction.
- J.H. Los et al., *"Scaling behavior and strain dependence of in-plane
  elastic properties of graphene,"* Phys. Rev. Lett. **116**, 015901 (2016).
  — Real-space atomistic MC benchmark (`η ≈ 0.84`, `ν ≈ +0.275`).

**External background** (cited in the thesis, not in `lib/`):
Le Doussal & Radzihovsky, review (2018); Nelson, Piran & Weinberg,
*Statistical Mechanics of Membranes and Surfaces* (2004);
Los et al., Phys. Rev. B **80**, 121405 (2009).
