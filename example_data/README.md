# `example_data/` — legacy large-N Fourier-MC output

These `N=<N>.dat` files are output of the **legacy** brane code
([`legacy/`](../legacy)), one file per linear size `N`, all produced at a single
physical coupling (`p8 = 0.3`, the legacy default). They are kept as a reference
dataset because they reach much larger lattices (`N = 100 … 200`,
`L = 2N+1 = 201 … 401`) than the quick modern runs, so the anomalous-elasticity
window is well resolved. See [`../docs/model.md`](../docs/model.md) for the physics.

## What a file contains

The lattice has `L = 2N+1` Fourier modes per axis. The writer
([`legacy/storage.c`](../legacy/storage.c), `dump()`) loops
`q1 = 0 … L-1`, `q2 = 0 … L-1` in **row-major order** and emits **three lines per
mode**:

```
c0   c1          # two ints: c0 = accepted-move counter, c1 = measurement count
Re(h_q)  Im(h_q) # last height mode value (complex); not needed for G(q)
g                # g = sum over measurements of |h_q|^2  (accumulator, not averaged)
```

After all `L*L` modes, one trailing line closes the file:

```
C    px0   px1   # global sample count and the two Poisson-ratio accumulators
```

So the file is `3 * L*L + 1` lines; as a flat stream of whitespace-separated
numbers it is `5 * L*L + 3` values (5 per mode: `c0 c1 Re Im g`, then `C px0 px1`).

There is **no header and no `p8` stored** — the size is only in the filename, and
the coupling is implied by how the run was launched (default `p8 = 0.3`).

## Reconstructing the Green function

`g` is the *sum* of `|h_q|^2` over `c1` measurements, so the height-height
correlator and its inverse are

```
G(q)    = <|h_q|^2> = g / c1
G^-1(q) = c1 / g
```

This is exactly the legacy definition (`storage.c:86`, where
`gy = c[1][0][i] / g[0][i]` is `G^-1` on the q-axis).

The mode indices are recovered from the row-major position `i = q1*L + q2`.
Frequencies are signed (`q1s = q1` if `q1 <= N`, else `q1 - L`), and with lattice
spacing `a = 2*pi/L` the continuum wavevector magnitude is

```
|q| = a * sqrt(q1s^2 + q2s^2)
```

matching the legacy axis convention `x = i*a` (and `i*a*sqrt(2)` on the diagonal).
The `(0,0)` zero mode has `g = 0` and is skipped.

### Difference from the modern format

The modern engine writes a single-line-per-mode, self-describing table
`q1 q2 qx qy qmag G Gerr Ginv` with a `#`-comment header carrying `N L p8 samples
sweeps nu nu_err rel_err converged`, and it stores an *averaged* `G` plus a
per-mode statistical error `Gerr` (from the replica spread). The legacy files
store raw accumulators (`g`, `c1`) with no error estimate, so a combined fit on
them falls back to the in-shell scatter for error bars.

## Using them

`tools/analyze.py` reads this format directly:

```sh
# combined multi-N pooled fit over all example_data files (one coupling):
uv run tools/analyze.py --legacy                 # -> plots/combined_legacy/
uv run tools/analyze.py --legacy --p8 0.3        # set the coupling (q8~p8) explicitly
```

Because every file shares the same physical coupling, the anomalous law
`G^-1 ~ q^(4-eta)` is size-independent, so the radially-averaged `G^-1(q)` points
from all N collapse onto one curve and can be pooled into a single log-log slope.
The larger lattices reach smaller `q`, widening the fit range and tightening
`eta`.

**Fit window matters.** The anomalous law holds only for `q << q8 ~ p8`; fitting
all the way up to `q8` itself includes the crossover roll-off and biases the slope
low (~0.68). The legacy analysis ([`legacy/plot.gp`](../legacy/plot.gp)) fits the
model `x^4 * (a*p8/x)^eta` (a pure power law `q^(4-eta)`) over the window
`[0.055, 0.11] = [q8/5.5, q8/2.7]`, well inside the plateau, and this is what our
`--legacy` mode reproduces:

```
eta = 0.786 +/- 0.006   (11 sizes N=100-200, 68 pts, window q in [0.055,0.11])
```

consistent with the thesis value `eta = 0.78 +/- 0.02`.

