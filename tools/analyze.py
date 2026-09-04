#!/usr/bin/env python3
"""
analyze.py -- Extract the anomalous elasticity exponent eta (and echo the
Poisson ratio) from a brane Fourier-Monte-Carlo output file, and plot the
inverse Green function.

Physics
-------
The height-height correlator is G(q) = <|h_q|^2>. Harmonic theory gives
G^-1 ~ q^4; anomalous elasticity gives G^-1 ~ q^(4-eta) for q << q8 ~ p8.
The two regimes join at the crossover scale q8.

Because the theory is isotropic (G depends only on q_r = |q|), we map every
mode (q_x, q_y) to q_r and *rotationally average* G over log-spaced |q| shells
(the theta average) -- using all L^2 modes, not just axes/diagonals.

Choosing the fit window
-----------------------
eta is *defined* as 4 - d ln G^-1/d ln q in the q->0 limit, so the principled,
assumption-light estimator is the plateau of the effective exponent

    eta_eff(q) = 4 - d ln G^-1 / d ln q

read over the scaling window (flat part of eta_eff, below the crossover q8 and
above finite-size effects). If eta_eff has no plateau, this lattice simply
cannot determine eta -- no single number should be forced. Three numbers are
reported:

  * plateau eta   -- PRIMARY: mean of eta_eff over the window, with its spread.
  * windowed slope-- straight log-log fit over the same window (for its error
    bar); equivalent to the plateau when the window is flat.
  * crossover fit -- CROSS-CHECK ONLY. Fit to the phenomenological ansatz
    G^-1 = C q^4 (1 + (q8/q)^eta), i.e. kappa(q)=kappa0[1+(q8/q)^eta] from the
    membrane literature. It captures both asymptotes but the crossover *shape*
    is not derived from theory, so with a limited q-range the eta it returns is
    model-dependent and can be biased -- treat it as indicative, not truth.

Only numpy is required. matplotlib is used if present (gnuplot fallback else).

Usage
-----
    python3 tools/analyze.py data/N=40.dat
    python3 tools/analyze.py data/N=40.dat --qmin 0.23 --qmax 0.40
"""
import argparse
import os
import re
import sys
import numpy as np


# ----------------------------------------------------------------------------
def load(path):
    header = {}
    saw_col_header = False
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "qmag" in line and "Ginv" in line:
                    saw_col_header = True
                for tok in line[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        header[k] = v
                continue
            break
    if not saw_col_header:
        sys.exit(
            f"error: '{path}' is not a brane (new-format) output file.\n"
            "  It looks like the legacy multi-line .dat format. Regenerate with\n"
            f"      ./brane N=40 p8=0.4 out={path}\n"
        )
    data = np.loadtxt(path, comments="#")
    # columns: q1 q2 qx qy qmag G Gerr Ginv   (Gerr added; older files lack it)
    if data.shape[1] >= 8:
        return data[:, 4], data[:, 5], data[:, 6], data[:, 7], header  # qmag,G,Gerr,Ginv
    # backward compat: old 7-column files (no Gerr)
    return data[:, 4], data[:, 5], np.zeros(len(data)), data[:, 6], header


def load_legacy(path):
    """Read a LEGACY brane file (example_data/N=<N>.dat).

    Legacy layout: L*L modes in row-major (q1=0..L-1, q2=0..L-1) order, each
    written over three lines -- "c0 c1", "Re(h) Im(h)", "g" -- followed by a
    trailing "C px0 px1" line. Here g = sum_measurements |h_q|^2 and c1 is the
    measurement count, so (legacy storage.c:86) the inverse Green function is
    G^-1(q) = c1 / g, i.e. G(q) = g / c1. The mode index (q1,q2) is recovered
    from the row-major position; signed frequency is q1 if q1<=N else q1-L, and
    |q| = a*sqrt(q1s^2 + q2s^2) with a = 2*pi/L (continuum convention, matching
    the legacy x = i*a axis).

    Returns (qmag, G, N, L, a) for the usable modes (g>0, c1>0, q>0).
    """
    N = int(re.search(r"N=(\d+)", path).group(1)); L = 2 * N + 1; a = 2 * np.pi / L
    toks = np.fromstring(open(path).read().replace("\t", " "), sep=" ")
    body = toks[:L * L * 5].reshape(L * L, 5)      # [c0, c1, re, im, g]
    c1, g = body[:, 1], body[:, 4]
    idx = np.arange(L * L); q1, q2 = idx // L, idx % L
    s1 = np.where(q1 <= N, q1, q1 - L)
    s2 = np.where(q2 <= N, q2, q2 - L)
    qmag = a * np.sqrt(s1.astype(float) ** 2 + s2.astype(float) ** 2)
    good = (g > 0) & (c1 > 0) & (qmag > 0)
    G = g[good] / c1[good]
    return qmag[good], G, N, L, a


def radial_average(qmag, G, nbins, Gerr=None):
    """Rotationally average G over log-spaced |q| shells.

    Returns (qr, Gr, Ginv_r, cnt, Ginv_err). Ginv_err is the propagated
    statistical error of 1/<G> in each shell: if per-mode errors Gerr (from the
    replica spread) are given, the shell error of the mean is
    sqrt(sum Gerr_i^2)/n, else it falls back to the in-shell std / sqrt(n).
    """
    mask = (qmag > 0) & (G > 0) & np.isfinite(G)
    q, g = qmag[mask], G[mask]
    ge = (Gerr[mask] if Gerr is not None else np.zeros_like(g))
    edges = np.logspace(np.log10(q.min()), np.log10(q.max()), nbins + 1)
    idx = np.digitize(q, edges)
    qr, gr, cnt, gerr = [], [], [], []
    for b in range(1, nbins + 1):
        sel = idx == b
        n = int(sel.sum())
        if n < 1:
            continue
        gm = g[sel].mean()
        # statistical error of the shell mean
        if Gerr is not None and np.any(ge[sel] > 0):
            em = np.sqrt(np.sum(ge[sel] ** 2)) / n
        else:
            em = (g[sel].std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        qr.append(q[sel].mean()); gr.append(gm); cnt.append(n); gerr.append(em)
    qr, gr, cnt, gerr = map(np.array, (qr, gr, cnt, gerr))
    Ginv_r = 1.0 / gr
    Ginv_err = gerr / gr ** 2          # error of 1/<G>
    return qr, gr, Ginv_r, cnt, Ginv_err


def effective_exponent(qr, Ginv_r):
    """eta_eff(q) = 4 - d ln(Ginv) / d ln(q) on the radial averages."""
    lq, lG = np.log(qr), np.log(Ginv_r)
    return qr, 4.0 - np.gradient(lG, lq)


def plateau_eta(qr, Ginv_r, cnt, qlo, qhi):
    """Primary, assumption-light estimator: the mean of eta_eff(q) over the
    scaling window [qlo, qhi] (which should be the flat part of eta_eff, below
    the crossover and above finite-size effects).

    eta is *defined* as 4 - d ln G^-1/d ln q in the q->0 limit, so averaging
    eta_eff across the plateau is the direct estimate. Returns (eta, spread,
    nshells); spread is the std of eta_eff across the window -- if it is large,
    there is no real plateau and eta is not determined by this lattice.
    """
    qe, ee = effective_exponent(qr, Ginv_r)
    m = (qe >= qlo) & (qe <= qhi) & np.isfinite(ee)
    if m.sum() < 3:
        return None, None, int(m.sum())
    w = np.sqrt(cnt[m])
    eta = float(np.average(ee[m], weights=w))
    spread = float(np.sqrt(np.average((ee[m] - eta) ** 2, weights=w)))
    return eta, spread, int(m.sum())


def fit_eta_window(qr, Ginv_r, cnt, qmin, qmax, Ginv_err=None):
    """Straight log-log fit over [qmin, qmax]. If per-shell statistical errors
    Ginv_err are given, weight by inverse variance in log space so the returned
    slope error is a genuine propagated statistical error; else weight by
    sqrt(shell count)."""
    m = (qr >= qmin) & (qr <= qmax) & (Ginv_r > 0) & np.isfinite(Ginv_r)
    if m.sum() < 3:
        return None, None, int(m.sum())
    x, y = np.log(qr[m]), np.log(Ginv_r[m])
    if Ginv_err is not None and np.all(Ginv_err[m] > 0):
        sigma = Ginv_err[m] / Ginv_r[m]          # error of log(Ginv)
        w = 1.0 / sigma
    else:
        w = np.sqrt(cnt[m])
    coef, cov = np.polyfit(x, y, 1, w=w, cov=True)
    return 4.0 - coef[0], float(np.sqrt(cov[0, 0])), int(m.sum())


def _crossover_logGinv(logq, logC, log_q8, eta):
    q, q8 = np.exp(logq), np.exp(log_q8)
    return logC + 4.0 * logq + np.log1p((q8 / q) ** eta)


def fit_eta_crossover(qr, Ginv_r, cnt, qlo, qhi, eta0=0.75, q8_0=None):
    """Fit G^-1 = C q^4 (1 + (q8/q)^eta) via pure-numpy Levenberg-Marquardt.

    Returns dict(eta, eta_err, q8, C, npts) or None.
    """
    m = (qr >= qlo) & (qr <= qhi) & (Ginv_r > 0) & np.isfinite(Ginv_r)
    if m.sum() < 5:
        return None
    x, y, w = np.log(qr[m]), np.log(Ginv_r[m]), np.sqrt(cnt[m])
    if q8_0 is None:
        q8_0 = float(np.sqrt(qr[m].min() * qr[m].max()))
    theta = np.array([0.0, np.log(q8_0), eta0])

    def resid(t):
        return w * (y - _crossover_logGinv(x, *t))

    def jac(t):
        J = np.zeros((len(x), 3)); eps = 1e-6
        f0 = _crossover_logGinv(x, *t)
        for k in range(3):
            tp = t.copy(); tp[k] += eps
            J[:, k] = -w * (_crossover_logGinv(x, *tp) - f0) / eps
        return J

    lam, cost = 1e-2, float(resid(theta) @ resid(theta))
    for _ in range(300):
        r, J = resid(theta), jac(theta)
        JTJ, JTr = J.T @ J, J.T @ r
        step = np.linalg.solve(JTJ + lam * np.diag(np.diag(JTJ)), -JTr)
        newt = theta + step
        newcost = float(resid(newt) @ resid(newt))
        if newcost < cost:
            theta, lam = newt, lam * 0.5
            if cost - newcost < 1e-12:
                cost = newcost; break
            cost = newcost
        else:
            lam *= 3.0
            if lam > 1e8:
                break
    J = jac(theta)
    dof = max(1, len(x) - 3)
    cov = np.linalg.inv(J.T @ J) * (cost / dof)
    return dict(eta=float(theta[2]), eta_err=float(np.sqrt(cov[2, 2])),
                q8=float(np.exp(theta[1])), C=float(np.exp(theta[0])),
                npts=int(m.sum()))


# ----------------------------------------------------------------------------
# Combined (multi-N) fit at fixed p8
# ----------------------------------------------------------------------------
# Different N have grid spacing a = 2*pi/L = 2*pi/(2N+1), so a larger N reaches
# *smaller* q. At a FIXED p8 (same coupling) the anomalous law G^-1 ~ q^(4-eta)
# is a bulk, N-independent property, so pooling the radial-averaged points from
# every N -- each kept only inside its own scaling window [3a_N, p8] -- gives a
# far denser, wider q-range than any single lattice and a tighter slope. We only
# ever pool cells with the SAME p8 (mixing couplings mixes crossover scales).
def collect_pooled(p8_target, pattern="data/N*/p*/data.dat", nbins=60,
                   tol=1e-3):
    """Pool radial-averaged (q_r, G^-1, err, N) points from every cell whose p8
    matches p8_target, each restricted to its own window [3 a_N, p8]."""
    import glob
    q, gi, ge, Ns = [], [], [], []
    for f in sorted(glob.glob(pattern)):
        qmag, G, Gerr, Ginv, header = load(f)
        p8 = float(header["p8"])
        if abs(p8 - p8_target) > tol:
            continue
        N = int(header["N"]); L = float(header["L"]); a = 2 * np.pi / L
        qr, Gr, Ginv_r, cnt, Ginv_err = radial_average(qmag, G, nbins, Gerr)
        m = (qr >= 3 * a) & (qr <= p8) & (Ginv_r > 0) & np.isfinite(Ginv_r)
        q.extend(qr[m]); gi.extend(Ginv_r[m])
        ge.extend(Ginv_err[m]); Ns.extend([N] * int(m.sum()))
    return (np.array(q), np.array(gi), np.array(ge), np.array(Ns, int))


def fit_pooled(q, gi, ge):
    """Inverse-variance-weighted log-log slope of the pooled points.
    Returns (eta, err, npts); eta is None if fewer than 3 usable points."""
    m = (q > 0) & (gi > 0) & np.isfinite(gi)
    if m.sum() < 3:
        return None, None, int(m.sum())
    x, y = np.log(q[m]), np.log(gi[m])
    if np.all(ge[m] > 0):
        sigma = ge[m] / gi[m]          # error of log(G^-1)
        w = 1.0 / sigma
    else:
        w = np.ones_like(x)
    coef, cov = np.polyfit(x, y, 1, w=w, cov=True)
    return 4.0 - coef[0], float(np.sqrt(cov[0, 0])), int(m.sum())


def plot_combined(q, gi, ge, Ns, eta, err, p8, png, bg=None, coupling=None):
    """Overlay every N's pooled points (colored by N) with the single combined
    slope fit q^(4-eta).

    p8 is the fit-window ceiling (points are fit over [q.min, p8]); coupling, if
    given, is the physical p8/q8 shown in the title (differs from the ceiling
    when the fit is capped below the crossover, as for the legacy data).

    If bg=(qbg, gibg) is given (full-Brillouin-zone radial averages, all q), it
    is drawn faintly so the whole G^-1(q) curve is visible -- the windowed fit
    region is shaded, and the q^4 / q^(4-eta) reference lines span the full q
    range (dashed outside the fit window to show the fit is only claimed there).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib missing; skipping combined plot.")
        return
    fig, ax = plt.subplots(figsize=(7.5, 6))
    if bg is not None:
        qbg, gibg = bg
        ax.loglog(qbg, gibg, ".", ms=2, color="0.75", alpha=0.5, zorder=0,
                  label="full BZ (all q)")
    uN = sorted(set(Ns.tolist()))
    cmap = plt.get_cmap("viridis")
    for i, n in enumerate(uN):
        s = Ns == n
        c = cmap(i / max(1, len(uN) - 1))
        ax.errorbar(q[s], gi[s], yerr=ge[s], fmt="o", ms=3.5, color=c,
                    ecolor=c, elinewidth=0.7, capsize=1.5, alpha=0.9,
                    label=f"N={n}")
    # reference lines: span the full plotted q-range if bg present, else window
    qmax_ref = float(max(q.max(), bg[0].max())) if bg is not None else float(q.max())
    qmin_ref = float(min(q.min(), bg[0].min())) if bg is not None else float(q.min())
    qq = np.logspace(np.log10(qmin_ref), np.log10(qmax_ref), 300)
    q0 = float(np.exp(np.mean(np.log(q)))); g0 = float(np.exp(np.mean(np.log(gi))))
    ax.loglog(qq, g0 * (qq / q0) ** 4, "--", color="C3", lw=1.2,
              label=r"$q^4$ (harmonic)")
    if eta is not None:
        # solid inside the fit window, dashed (extrapolation) outside
        inw = (qq >= q.min()) & (qq <= p8)
        gline = g0 * (qq / q0) ** (4 - eta)
        ax.loglog(qq[~inw], gline[~inw], ":", color="k", lw=1.2, alpha=0.6)
        ax.loglog(qq[inw], gline[inw], "-", color="k", lw=2,
                  label=rf"combined fit $\eta={eta:.3f}\pm{err:.3f}$")
        ax.axvspan(q.min(), p8, color="C1", alpha=0.10, label="fit window")
        ax.axvline(p8, ls=":", color="C1", alpha=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$q_r=|q|$"); ax.set_ylabel(r"$G^{-1}(q_r)$")
    if coupling is not None:
        ttl = (rf"Combined multi-$N$ fit, $q_8\!\sim\!p_8={coupling:.2f}$, "
               rf"window $q\leq{p8:.3f}$  ({len(uN)} sizes, {len(q)} pts)")
    else:
        ttl = (rf"Combined multi-$N$ fit at $p_8={p8:.2f}$  "
               rf"({len(uN)} sizes, {len(q)} pts)")
    ax.set_title(ttl)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    print(f"[plot] wrote {png}")


def combined_all(pattern="data/N*/p*/data.dat", nbins=60, outdir="plots/combined"):
    """Run the pooled fit for every distinct p8 present; write one overlay per
    p8, a summary eta(p8) curve, and a CSV. Returns the list of result dicts."""
    import glob
    p8s = sorted({round(float(load(f)[4]["p8"]), 4)
                  for f in glob.glob(pattern)})
    if not p8s:
        sys.exit(f"no files match {pattern}")
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for p8 in p8s:
        q, gi, ge, Ns = collect_pooled(p8, pattern, nbins)
        if len(q) < 3:
            print(f"  p8={p8:.2f}: only {len(q)} pooled pts, skipping")
            continue
        eta, err, npts = fit_pooled(q, gi, ge)
        nN = len(set(Ns.tolist()))
        rows.append(dict(p8=p8, eta=eta, err=err, npts=npts, nN=nN))
        estr = f"{eta:.3f} +/- {err:.3f}" if eta is not None else "n/a"
        print(f"  p8={p8:.2f}: eta={estr}  ({nN} sizes, {npts} pts)")
        plot_combined(q, gi, ge, Ns, eta, err, p8,
                      os.path.join(outdir, f"p{p8:.2f}.png"))
    # summary eta(p8)
    csv = os.path.join(outdir, "combined_eta.csv")
    with open(csv, "w") as f:
        f.write("p8,eta,err,nN,npts\n")
        for r in rows:
            f.write(f"{r['p8']:.2f},{r['eta']:.4f},{r['err']:.4f},"
                    f"{r['nN']},{r['npts']}\n")
    print(f"[csv] wrote {csv}")
    _plot_eta_of_p8(rows, os.path.join(outdir, "eta_vs_p8.png"))
    return rows


def combined_legacy(pattern="example_data/N=*.dat", p8=0.3,
                    lo_frac=0.055 / 0.3, hi_frac=0.11 / 0.3,
                    nbins=40, outdir="plots/combined_legacy"):
    """Combined multi-N pooled fit on the LEGACY large-N data (all files share
    one physical coupling, default p8=0.3, so q8~p8=0.3).

    The anomalous law G^-1 ~ q^(4-eta) holds only for q << q8. Fitting up to q8
    itself includes the crossover roll-off and biases the slope LOW (~0.68).
    The legacy analysis (legacy/plot.gp) fits the model x^4*(a*p8/x)^eta -- a
    pure power law q^(4-eta) -- over the window [0.055, 0.11] = [q8/5.5, q8/2.7],
    well inside the plateau, and reports eta ~ 0.78. We reproduce that: pool the
    radial-averaged G^-1(q) from every N over [lo_frac*q8, hi_frac*q8] (default
    the legacy [0.055,0.11]) and fit one count-weighted log-log slope. Returns
    (eta, err, npts, nN)."""
    import glob
    files = sorted(glob.glob(pattern),
                   key=lambda p: int(re.search(r"N=(\d+)", p).group(1)))
    if not files:
        sys.exit(f"no legacy files match {pattern}")
    os.makedirs(outdir, exist_ok=True)
    lo, hi = lo_frac * p8, hi_frac * p8
    q, gi, ge, Ns = [], [], [], []
    qbg, gibg = [], []                       # full-BZ background for the plot
    for f in files:
        qmag, G, N, L, a = load_legacy(f)
        qr, Gr, Ginv_r, cnt, Ginv_err = radial_average(qmag, G, nbins)
        qbg.extend(qr); gibg.extend(Ginv_r)
        m = (qr >= lo) & (qr <= hi) & (Ginv_r > 0) & np.isfinite(Ginv_r)
        q.extend(qr[m]); gi.extend(Ginv_r[m])
        ge.extend(Ginv_err[m]); Ns.extend([N] * int(m.sum()))
        print(f"  N={N:3d}  a={a:.4f}  window=[{lo:.3f},{hi:.3f}]  "
              f"pts={int(m.sum())}")
    q, gi, ge, Ns = map(np.array, (q, gi, ge, Ns))
    qbg, gibg = np.array(qbg), np.array(gibg)
    eta, err, npts = fit_pooled(q, gi, ge)
    nN = len(set(Ns.tolist()))
    print(f"\nCOMBINED (legacy, q8~p8={p8}, window=[{lo:.3f},{hi:.3f}]"
          f"=[q8/{p8/lo:.1f}, q8/{p8/hi:.1f}]): "
          f"eta = {eta:.3f} +/- {err:.3f}   "
          f"({nN} sizes N={Ns.min()}-{Ns.max()}, {npts} pts)")
    plot_combined(q, gi, ge, Ns, eta, err, hi,
                  os.path.join(outdir, "combined_legacy.png"),
                  bg=(qbg, gibg), coupling=p8)
    return eta, err, npts, nN


def _plot_eta_of_p8(rows, png):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    good = [r for r in rows if r["eta"] is not None]
    if not good:
        return
    p = np.array([r["p8"] for r in good])
    e = np.array([r["eta"] for r in good])
    er = np.array([r["err"] for r in good])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(p, e, yerr=er, fmt="o-", color="C0", capsize=3,
                label=r"combined multi-$N$ fit")
    ax.set_xlabel(r"$p_8$"); ax.set_ylabel(r"combined-fit $\eta$")
    ax.set_title(r"Pooled-over-$N$ exponent $\eta(p_8)$")
    ax.grid(alpha=0.3); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    print(f"[plot] wrote {png}")


# ----------------------------------------------------------------------------
def plot(qr, Ginv_r, cnt, qmag, Ginv, eta_w, cross, p8, qmin, qmax,
         out_png, out_gp, datfile, Ginv_err=None):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        eta = cross["eta"] if cross else (eta_w or 0.78)
        with open(out_gp, "w") as f:
            f.write(f"""set terminal pngcairo size 900,650
set output '{out_png}'
set logscale xy
set xlabel 'q_r'; set ylabel 'G^{{-1}}'
plot '{datfile}' using 5:7 pt 7 ps 0.3 lc rgb '#cccccc' t 'all modes', \\
     x**4 lw 2 t 'q^4', x**(4-{eta}) lw 2 t 'q^{{4-eta}}'
""")
        print(f"[plot] matplotlib not found; wrote gnuplot script {out_gp}")
        return

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5))
    # left: inverse Green
    axL.loglog(qmag, Ginv, ".", ms=2, alpha=0.18, color="0.7",
               label="all modes (q_x, q_y)")
    if Ginv_err is not None and np.any(Ginv_err > 0):
        axL.errorbar(qr, Ginv_r, yerr=Ginv_err, fmt="o", ms=2.5, color="C0",
                     ecolor="k", elinewidth=0.8, capsize=2, zorder=3,
                     label=r"radial avg $G^{-1}(q_r)\pm$SE")
    else:
        axL.scatter(qr, Ginv_r, s=12 + 40 * cnt / cnt.max(), color="C0",
                    zorder=3, label=r"radial avg $G^{-1}(q_r)$")
    qq = np.logspace(np.log10(qmag[qmag > 0].min()), np.log10(qmag.max()), 200)
    axL.loglog(qq, qq**4, "--", color="C3", lw=1.3, label=r"$q^4$ (harmonic)")
    if cross:
        eta, q8, C = cross["eta"], cross["q8"], cross["C"]
        axL.loglog(qq, C * qq**4 * (1 + (q8 / qq) ** eta), "-", color="C2",
                   lw=2, label=rf"crossover ansatz $\eta={eta:.2f}$ (heuristic)")
        axL.axvline(q8, ls=":", color="C2", alpha=0.7,
                    label=rf"$q_8={q8:.2f}$")
    axL.axvspan(qmin, qmax, color="C1", alpha=0.10, label="windowed-fit range")
    axL.set_xlabel(r"$q_r=|q|$"); axL.set_ylabel(r"$G^{-1}(q_r)$")
    axL.set_title(rf"Inverse Green ($p_8={p8}$)")
    axL.legend(frameon=False, fontsize=8); axL.grid(True, which="both", alpha=0.3)

    # right: effective exponent
    qe, ee = effective_exponent(qr, Ginv_r)
    axR.semilogx(qe, ee, "o-", ms=4, color="C0", label=r"$\eta_{eff}(q)$")
    if cross:
        axR.axhline(cross["eta"], ls="-", color="C2",
                    label=rf"crossover $\eta={cross['eta']:.2f}$")
        axR.axvline(cross["q8"], ls=":", color="C2", alpha=0.7)
    if eta_w is not None:
        axR.axhline(eta_w, ls="--", color="C1", label=rf"windowed $\eta={eta_w:.2f}$")
    axR.axhline(0.0, color="0.8", lw=0.8)
    axR.axvspan(qmin, qmax, color="C1", alpha=0.10)
    axR.set_ylim(-0.5, 4.2)
    axR.set_xlabel(r"$q_r$"); axR.set_ylabel(r"$\eta_{eff}=4-d\ln G^{-1}/d\ln q$")
    axR.set_title("Effective exponent (plateau = scaling window)")
    axR.legend(frameon=False, fontsize=8); axR.grid(True, which="both", alpha=0.3)

    fig.tight_layout(); fig.savefig(out_png, dpi=140)
    plt.close(fig)
    print(f"[plot] wrote {out_png}")


# ----------------------------------------------------------------------------
def analyze_file(datfile, qmin_arg=None, qmax_arg=None, nbins=60, quv=1.0,
                 png=None, quiet=False):
    qmag, G, Gerr, Ginv, header = load(datfile)
    p8 = float(header.get("p8", 0.3))
    a = 2 * np.pi / float(header.get("L", 81))

    # NOTE: default window stops at the crossover q8~p8 (not above it). The
    # right-hand eta_eff(q) panel is there so you can judge whether even p8 is
    # too high a ceiling and tighten --qmax below the crossover.
    qmin = qmin_arg if qmin_arg is not None else 3 * a
    qmax = qmax_arg if qmax_arg is not None else p8

    qr, Gr, Ginv_r, cnt, Ginv_err = radial_average(qmag, G, nbins, Gerr)
    eta_p, spread_p, nsh_p = plateau_eta(qr, Ginv_r, cnt, qmin, qmax)
    eta_w, err_w, nsh = fit_eta_window(qr, Ginv_r, cnt, qmin, qmax, Ginv_err)
    cross = fit_eta_crossover(qr, Ginv_r, cnt, 2 * a, quv,
                              eta0=eta_w if eta_w else 0.75, q8_0=p8)

    if not quiet:
        print(f"file           : {datfile}")
        print(f"N, L           : {header.get('N','?')}, {header.get('L','?')}")
        print(f"p8             : {p8}")
        print(f"window         : q_r in [{qmin:.3f}, {qmax:.3f}]  ({nsh} shells)")
        if eta_p is not None:
            flat = "flat" if spread_p < 0.25 else "NOT flat -> eta ill-defined here"
            print(f"plateau eta    : {eta_p:.3f}  (eta_eff spread {spread_p:.3f}, {flat})  [PRIMARY]")
        if eta_w is not None:
            print(f"windowed slope : {eta_w:.3f} +/- {err_w:.3f}  (stat. error)")
        if cross:
            print(f"crossover fit  : {cross['eta']:.3f} +/- {cross['eta_err']:.3f}"
                  f"   q8={cross['q8']:.3f}   [heuristic ansatz, cross-check only]")
        print(f"Poisson ratio  : {header.get('nu','?')} +/- {header.get('nu_err','?')}")
        conv = header.get("converged", "?")
        print(f"run            : sweeps={header.get('sweeps','?')} "
              f"Delta2 rel.err={header.get('rel_err','?')} converged={conv}")

    # mirror the data subpath: data/N40/p0.40/data.dat -> plots/N40/p0.40/
    d = os.path.dirname(datfile)
    pdir = ("plots" + d[len("data"):]) if (d == "data" or d.startswith("data/")) else "plots"
    os.makedirs(pdir, exist_ok=True)
    png = png or os.path.join(pdir, "fit.png")
    gp = os.path.join(pdir, "fit.gp")
    plot(qr, Ginv_r, cnt, qmag, Ginv, eta_w, cross, p8, qmin, qmax, png, gp, datfile,
         Ginv_err=Ginv_err)
    return dict(N=int(header.get("N", 0)), p8=p8, png=png,
                eta_plateau=eta_p, eta_plateau_spread=spread_p,
                eta_window=eta_w, eta_window_err=err_w,
                eta_crossover=cross["eta"] if cross else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datfile", nargs="?", help="a .dat file; omit with --all")
    ap.add_argument("--qmin", type=float, default=None)
    ap.add_argument("--qmax", type=float, default=None)
    ap.add_argument("--nbins", type=int, default=60, help="radial shells")
    ap.add_argument("--quv", type=float, default=1.0,
                    help="upper q for crossover fit (avoid lattice UV)")
    ap.add_argument("--png", default=None)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--all", metavar="GLOB", nargs="?", const="data/N*/p*/data.dat",
                    help="batch-analyze every matching file (default data/N*/p*/data.dat)")
    ap.add_argument("--combined", action="store_true",
                    help="pool all N at each p8 into one weighted log-log fit; "
                         "writes plots/combined/{p<p8>.png, eta_vs_p8.png, combined_eta.csv}")
    ap.add_argument("--combined-p8", type=float, default=None,
                    help="combined fit for a single p8 only (implies --combined)")
    ap.add_argument("--legacy", metavar="GLOB", nargs="?",
                    const="example_data/N=*.dat",
                    help="combined multi-N fit on LEGACY files "
                         "(default example_data/N=*.dat); set coupling with --p8")
    ap.add_argument("--p8", type=float, default=0.3,
                    help="physical coupling of the legacy runs (fit ceiling); default 0.3")
    args = ap.parse_args()

    if args.legacy:
        combined_legacy(args.legacy, p8=args.p8, nbins=args.nbins)
        return

    if args.combined or args.combined_p8 is not None:
        if args.combined_p8 is not None:
            q, gi, ge, Ns = collect_pooled(args.combined_p8, nbins=args.nbins)
            if len(q) < 3:
                sys.exit(f"only {len(q)} pooled points at p8={args.combined_p8}")
            eta, err, npts = fit_pooled(q, gi, ge)
            nN = len(set(Ns.tolist()))
            print(f"p8={args.combined_p8:.2f}: eta={eta:.3f} +/- {err:.3f}  "
                  f"({nN} sizes, {npts} pts)")
            os.makedirs("plots/combined", exist_ok=True)
            plot_combined(q, gi, ge, Ns, eta, err, args.combined_p8,
                          f"plots/combined/p{args.combined_p8:.2f}.png")
        else:
            combined_all(nbins=args.nbins)
        return

    if args.all:
        import glob
        files = sorted(glob.glob(args.all))
        if not files:
            sys.exit(f"no files match {args.all}")
        for i, f in enumerate(files, 1):
            try:
                r = analyze_file(f, args.qmin, args.qmax, args.nbins, args.quv,
                                 quiet=True)
                print(f"[{i}/{len(files)}] {r['png']}  "
                      f"p8={r['p8']:.2f} N={r['N']} etaW={r['eta_window']}")
            except Exception as e:
                print(f"[{i}/{len(files)}] skip {f}: {e}")
        return

    if not args.datfile:
        sys.exit("give a .dat file, or use --all")
    analyze_file(args.datfile, args.qmin, args.qmax, args.nbins, args.quv,
                 args.png, args.quiet)


if __name__ == "__main__":
    main()
