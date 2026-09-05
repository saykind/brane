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


def running_eta(q, Ginv, nwin=40, width=1.6, min_pts=15):
    """eta_eff(q_r) WITHOUT angular averaging: at each q_r scale, fit the local
    log-log slope of ALL modes whose |q| falls in a fixed multiplicative window
    [c/width, c*width]. This is a local *fit* across every angle at that scale
    (not an average of G), so anisotropy is never collapsed into one value.
    Returns (q_centers, eta_eff)."""
    m = (q > 0) & (Ginv > 0) & np.isfinite(Ginv)
    q, Ginv = q[m], Ginv[m]
    if len(q) < min_pts:
        return np.array([]), np.array([])
    lq, lg = np.log(q), np.log(Ginv)
    centers = np.logspace(np.log10(q.min() * 1.02), np.log10(q.max() * 0.98), nwin)
    hw = np.log(width)
    cen, eta = [], []
    for c in centers:
        lc = np.log(c)
        sel = (lq >= lc - hw) & (lq <= lc + hw)
        if int(sel.sum()) < min_pts:
            continue
        slope = np.polyfit(lq[sel], lg[sel], 1)[0]
        cen.append(c); eta.append(4.0 - slope)
    return np.array(cen), np.array(eta)


def plateau_lowq(qc, ee, tol=0.05, minpts=5):
    """Extract the low-q plateau of the running exponent eta_eff(q_r): the true
    q->0 anomalous exponent. Finds the WIDEST contiguous window of eta_eff whose
    scatter (std) is below `tol`, preferring lower q. The erratic super-low-q
    points (finite-size breakdown) and the high-q crossover rise both fail the
    flatness test and are excluded automatically. Returns
    (eta, err, (q_lo, q_hi), npts) or None."""
    m = np.isfinite(ee) & (qc > 0)
    qc, ee = qc[m], ee[m]
    o = np.argsort(qc); qc, ee = qc[o], ee[o]
    n = len(qc)
    if n < minpts:
        return None
    best = None                       # ((npts, -q_lo), i, j)
    for i in range(n):
        for j in range(i + minpts, n + 1):
            seg = ee[i:j]
            if seg.std() <= tol:
                cand = (j - i, -qc[i])
                if best is None or cand > best[0]:
                    best = (cand, i, j)
    if best is None:                  # nothing flat enough: take flattest fixed window
        W = max(minpts, n // 3)
        _, i = min((ee[k:k + W].std(), k) for k in range(0, n - W + 1)); j = i + W
    else:
        _, i, j = best
    seg = ee[i:j]
    return (float(seg.mean()), float(seg.std() / np.sqrt(len(seg))),
            (float(qc[i]), float(qc[j - 1])), int(len(seg)))


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
def collect_pooled(p8_target, pattern="data/N*/p*/*/*.dat", nbins=60,
                   tol=1e-3):
    """Pool EVERY mode (q_r, G^-1, err, N) -- no angular averaging -- from every
    cell whose p8 matches p8_target, each restricted to its own window
    [3 a_N, p8]. (nbins kept for signature compatibility; unused.)"""
    import glob
    q, gi, ge, Ns = [], [], [], []
    for f in sorted(glob.glob(pattern)):
        qmag, G, Gerr, Ginv, header = load(f)
        p8 = float(header["p8"])
        if abs(p8 - p8_target) > tol:
            continue
        N = int(header["N"]); L = float(header["L"]); a = 2 * np.pi / L
        gv = (G > 0) & np.isfinite(G) & (qmag > 0)
        uq = np.unique(qmag[gv])
        qlo = uq[2] if len(uq) >= 3 else (uq[0] if len(uq) else 3 * a)  # 3rd smallest |q|
        m = gv & (qmag >= qlo) & (qmag <= p8)
        q.extend(qmag[m]); gi.extend(1.0 / G[m])
        ge.extend(np.where(Gerr[m] > 0, Gerr[m] / G[m] ** 2, 0.0))
        Ns.extend([N] * int(m.sum()))
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


def plot_combined(q, gi, ge, Ns, eta, err, png, wlo, whi, coupling=None):
    """Two panels, no angular averaging:
      left  - every mode (q_r, Ginv) colored by N + windowed power-law fit
      right - running exponent eta_eff(q_r) of the pooled cloud, with the
              auto-detected LOW-Q PLATEAU (the true q->0 exponent) shaded.
    Returns the plateau (eta, err, (q_lo,q_hi), npts) or None."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib missing; skipping combined plot.")
        return None
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.5))
    uN = sorted(set(Ns.tolist()))
    cmap = plt.get_cmap("viridis")
    for i, n in enumerate(uN):
        s = Ns == n
        c = cmap(i / max(1, len(uN) - 1))
        axL.scatter(q[s], gi[s], s=4, alpha=0.35, color=c, linewidths=0,
                    rasterized=True, label=f"N={n}")
    qq = np.logspace(np.log10(q.min()), np.log10(q.max()), 300)
    inpts = (q >= wlo) & (q <= whi)
    q0 = float(np.exp(np.mean(np.log(q[inpts]))))
    g0 = float(np.exp(np.mean(np.log(gi[inpts]))))
    axL.loglog(qq, g0 * (qq / q0) ** 4, "--", color="C3", lw=1.2,
               label=r"$q^4$ (harmonic)")
    if eta is not None:
        inw = (qq >= wlo) & (qq <= whi)
        gline = g0 * (qq / q0) ** (4 - eta)
        axL.loglog(qq[~inw], gline[~inw], ":", color="k", lw=1.2, alpha=0.6)
        axL.loglog(qq[inw], gline[inw], "-", color="k", lw=2,
                   label=rf"window fit $\eta={eta:.3f}$")
        axL.axvspan(wlo, whi, color="C1", alpha=0.12, label="fit window")
    axL.set_xscale("log"); axL.set_yscale("log")
    axL.set_xlabel(r"$q_r=|q|$"); axL.set_ylabel(r"$G^{-1}(q_r)$  (every mode)")
    ttl = (rf"Combined multi-$N$, $p_8={coupling if coupling is not None else whi:.2f}$"
           rf"  ({len(uN)} sizes, {len(q)} modes)")
    axL.set_title(ttl)
    axL.legend(frameon=False, fontsize=7, ncol=2); axL.grid(True, which="both", alpha=0.3)

    # right: running exponent + low-q plateau
    qc, ee = running_eta(q, gi)
    pl = plateau_lowq(qc, ee) if len(qc) else None
    if len(qc):
        axR.semilogx(qc, ee, "-", color="C0", lw=1.5, label=r"$\eta_{eff}(q_r)$")
    if eta is not None:
        axR.axhline(eta, ls="--", color="C1", alpha=0.7,
                    label=rf"window fit $\eta={eta:.3f}$")
    if pl is not None:
        etp, erp, (qlo, qhi), npl = pl
        axR.axvspan(qlo, qhi, color="C2", alpha=0.15)
        axR.axhline(etp, ls="-", color="C2", lw=2,
                    label=rf"low-$q$ plateau $\eta={etp:.3f}\pm{erp:.3f}$")
    axR.axhline(0.0, color="0.8", lw=0.8)
    axR.set_ylim(-0.5, 4.2)
    axR.set_xlabel(r"$q_r$"); axR.set_ylabel(r"$\eta_{eff}=4-d\ln G^{-1}/d\ln q$")
    axR.set_title("Running exponent + low-$q$ plateau")
    axR.legend(frameon=False, fontsize=8); axR.grid(True, which="both", alpha=0.3)

    fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    print(f"[plot] wrote {png}")
    return pl


def combined_all(pattern="data/N*/p*/*/*.dat", nbins=60, outdir="plots/combined"):
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
        pl = plot_combined(q, gi, ge, Ns, eta, err,
                           os.path.join(outdir, f"p{p8:.2f}.png"),
                           float(q.min()), p8)
        etp = pl[0] if pl else None
        erp = pl[1] if pl else None
        rows.append(dict(p8=p8, eta=eta, err=err, npts=npts, nN=nN,
                         eta_plateau=etp, eta_plateau_err=erp))
        estr = f"{eta:.3f}+/-{err:.3f}" if eta is not None else "n/a"
        pstr = f"{etp:.3f}+/-{erp:.3f}" if etp is not None else "n/a"
        print(f"  p8={p8:.2f}: low-q plateau eta={pstr} [PRIMARY]  "
              f"window eta={estr}  ({nN} sizes, {npts} pts)")
    # summary eta(p8)
    csv = os.path.join(outdir, "combined_eta.csv")
    with open(csv, "w") as f:
        f.write("p8,eta_plateau,eta_plateau_err,eta_window,eta_window_err,nN,npts\n")
        for r in rows:
            ep = r.get("eta_plateau"); epe = r.get("eta_plateau_err")
            f.write(f"{r['p8']:.2f},{ep if ep is None else round(ep,4)},"
                    f"{epe if epe is None else round(epe,4)},"
                    f"{r['eta']:.4f},{r['err']:.4f},{r['nN']},{r['npts']}\n")
    print(f"[csv] wrote {csv}")
    _plot_eta_of_p8(rows, os.path.join(outdir, "eta_vs_p8.png"))
    return rows


def combined_legacy(pattern="example_data/N=*.dat", p8=0.3,
                    lo=0.055, hi=0.11, nbins=40, outdir="plots/combined_legacy"):
    """Combined multi-N pooled fit on the LEGACY large-N data (all files share
    one physical coupling, default p8=0.3, so q8~p8=0.3).

    The fit window is the exact legacy one: q in [0.055, 0.11], hardcoded in
    legacy/plot.gp (fit2=0.11, fit1=.5*fit2) as an ABSOLUTE q-range. The model
    there, x^4*(a*p8/x)^eta, is a pure power law q^(4-eta), so this is just a
    log-log slope over that window -- well below the crossover q8~p8=0.3, where
    the anomalous plateau lives. Each N is radially averaged; the in-window
    points are pooled into one count-weighted slope. Returns (eta, err, npts,
    nN)."""
    import glob
    files = sorted(glob.glob(pattern),
                   key=lambda p: int(re.search(r"N=(\d+)", p).group(1)))
    if not files:
        sys.exit(f"no legacy files match {pattern}")
    os.makedirs(outdir, exist_ok=True)
    q, gi, ge, Ns = [], [], [], []          # ALL points (colored by N in plot)
    per_n = []                               # (N, eta_w, err_w) per file
    for f in files:
        qmag, G, N, L, a = load_legacy(f)
        qr, Gr, Ginv_r, cnt, Ginv_err = radial_average(qmag, G, nbins)
        q.extend(qr); gi.extend(Ginv_r); ge.extend(Ginv_err)
        Ns.extend([N] * len(qr))
        inw = int(((qr >= lo) & (qr <= hi)).sum())
        e1, er1, _ = fit_eta_window(qr, Ginv_r, cnt, lo, hi, Ginv_err)
        if e1 is not None:
            per_n.append((N, e1, er1))
        print(f"  N={N:3d}  a={a:.4f}  window=[{lo:.3f},{hi:.3f}]  in-window pts={inw}"
              f"  eta_N={'%.3f'%e1 if e1 is not None else 'n/a'}")
        # per-file fit.png (mirrors the modern analyze fit plot)
        analyze_legacy_file(f, p8, lo, hi, qmag, 1.0 / G, qr, Gr, cnt, Ginv_err)
    q, gi, ge, Ns = map(np.array, (q, gi, ge, Ns))
    # fit uses only in-window points
    inw = (q >= lo) & (q <= hi)
    eta, err, npts = fit_pooled(q[inw], gi[inw], ge[inw])
    nN = len(set(Ns.tolist()))
    print(f"\nCOMBINED (legacy, q8~p8={p8}, legacy window=[{lo:.3f},{hi:.3f}]): "
          f"eta = {eta:.3f} +/- {err:.3f}   "
          f"({nN} sizes N={Ns.min()}-{Ns.max()}, {npts} in-window pts)")
    plot_combined(q, gi, ge, Ns, eta, err,
                  os.path.join(outdir, "combined_legacy.png"),
                  lo, hi, coupling=p8)
    _plot_eta_of_N(per_n, eta, err, lo, hi, p8,
                   os.path.join(outdir, "eta_vs_N.png"))
    return eta, err, npts, nN


def _plot_eta_of_N(per_n, eta_pooled, err_pooled, lo, hi, p8, png):
    """eta measured per single N (windowed fit over [lo,hi]) vs N -- tests
    whether the exponent is N-independent (universal) at these large sizes."""
    if not per_n:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    N = np.array([r[0] for r in per_n], float)
    e = np.array([r[1] for r in per_n], float)
    er = np.array([r[2] for r in per_n], float)
    # trend test: weighted linear regression eta vs N
    w = 1.0 / np.where(er > 0, er, np.nan) ** 2
    ok = np.isfinite(w)
    slope = serr = None
    if ok.sum() >= 3:
        W = np.diag(w[ok]); X = np.vstack([N[ok], np.ones(ok.sum())]).T
        cov = np.linalg.inv(X.T @ W @ X)
        beta = cov @ X.T @ W @ e[ok]
        slope, serr = float(beta[0]), float(np.sqrt(cov[0, 0]))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.errorbar(N, e, yerr=er, fmt="o", color="C0", capsize=3,
                label=r"per-$N$ windowed $\eta$")
    ax.axhspan(eta_pooled - err_pooled, eta_pooled + err_pooled,
               color="C1", alpha=0.20)
    ax.axhline(eta_pooled, color="C1", lw=1.5,
               label=rf"pooled $\eta={eta_pooled:.3f}\pm{err_pooled:.3f}$")
    ax.axhline(0.78, ls=":", color="0.4", label=r"thesis $\eta=0.78$")
    ttl = rf"Legacy $\eta$ vs $N$  (window $[{lo:.3f},{hi:.3f}]$, $q_8\!\sim\!p_8={p8}$)"
    if slope is not None:
        sig = abs(slope) / serr if serr else 0.0
        ax.set_title(ttl + "\n" +
                     rf"trend $d\eta/dN={slope:+.4f}\pm{serr:.4f}$ "
                     rf"({sig:.1f}$\sigma$)")
    else:
        ax.set_title(ttl)
    ax.set_xlabel(r"$N$  ($L=2N+1$)"); ax.set_ylabel(r"windowed $\eta$")
    ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    print(f"[plot] wrote {png}"
          + (f"   (d eta/dN = {slope:+.4f} +/- {serr:.4f})" if slope is not None else ""))


def analyze_legacy_file(datfile, p8, lo, hi, qmag, Ginv_all, qr, Gr, cnt,
                        Ginv_err, outdir="plots/legacy"):
    """Per-file fit.png for a legacy example_data file: inverse Green + effective
    exponent, fit over the sub-crossover window [lo, hi]. Writes
    plots/legacy/N<N>/fit.png."""
    N = int(re.search(r"N=(\d+)", datfile).group(1))
    Ginv_r = 1.0 / Gr
    eta_w, err_w, nsh = fit_eta_window(qr, Ginv_r, cnt, lo, hi, Ginv_err)
    cross = fit_eta_crossover(qr, Ginv_r, cnt, lo, 1.0, eta0=eta_w or 0.78,
                              q8_0=p8)
    pdir = os.path.join(outdir, f"N{N}")
    os.makedirs(pdir, exist_ok=True)
    png = os.path.join(pdir, "fit.png")
    gp = os.path.join(pdir, "fit.gp")
    plot(qr, Ginv_r, cnt, qmag, Ginv_all, eta_w, cross, p8, lo, hi, png, gp,
         datfile, Ginv_err=Ginv_err)


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
def plot(q, Ginv, theta, qcen, eta_eff, eta_w, fitline, p8, qmin, qmax,
         out_png, out_gp, datfile, plateau=None):
    """No angular averaging: every mode is a point (q_r, Ginv), colored by its
    polar angle so direction-dependence (anisotropy) is visible. Left panel is
    the point cloud + windowed power-law fit; right panel is the running local
    exponent eta_eff(q_r)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        eta = eta_w if eta_w else 0.78
        with open(out_gp, "w") as f:
            f.write(f"""set terminal pngcairo size 900,650
set output '{out_png}'
set logscale xy
set xlabel 'q_r'; set ylabel 'G^{{-1}}'
plot '{datfile}' using 5:8 pt 7 ps 0.3 lc rgb '#cccccc' t 'all modes', \\
     x**4 lw 2 t 'q^4', x**(4-{eta}) lw 2 t 'q^{{4-eta}}'
""")
        print(f"[plot] matplotlib not found; wrote gnuplot script {out_gp}")
        return

    # fold angle into [0,45] deg by lattice symmetry: 0/90 = axis, 45 = diagonal
    tdeg = np.degrees(np.arctan2(np.abs(np.sin(theta)), np.abs(np.cos(theta))))
    tdeg = np.minimum(tdeg, 90 - tdeg)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5))

    sc = axL.scatter(q, Ginv, c=tdeg, s=4, alpha=0.45, cmap="viridis",
                     vmin=0, vmax=45, rasterized=True, linewidths=0)
    cb = fig.colorbar(sc, ax=axL, pad=0.01)
    cb.set_label("mode angle to nearest axis (deg): 0=axis, 45=diagonal", fontsize=8)
    qq = np.logspace(np.log10(q[q > 0].min()), np.log10(q.max()), 200)
    axL.plot(qq, qq ** 4, "--", color="C3", lw=1.3, label=r"$q^4$ (harmonic)")
    if fitline is not None:
        axL.plot(fitline[0], fitline[1], "-", color="k", lw=2,
                 label=rf"windowed fit $\eta={eta_w:.2f}$")
    axL.axvspan(qmin, qmax, color="C1", alpha=0.12, label="fit window")
    axL.set_xscale("log"); axL.set_yscale("log")
    axL.set_xlabel(r"$q_r=|q|$"); axL.set_ylabel(r"$G^{-1}(q)$  (every mode)")
    axL.set_title(rf"Inverse Green, all modes ($p_8={p8}$)")
    axL.legend(frameon=False, fontsize=8); axL.grid(True, which="both", alpha=0.3)

    if len(qcen):
        axR.semilogx(qcen, eta_eff, "-", color="C0", lw=1.5,
                     label=r"$\eta_{eff}(q_r)$ (sliding local fit)")
    if eta_w is not None:
        axR.axhline(eta_w, ls="--", color="C1", label=rf"windowed $\eta={eta_w:.2f}$")
    if plateau is not None:
        etp, erp, (qlo, qhi), _ = plateau
        axR.axvspan(qlo, qhi, color="C2", alpha=0.15)
        axR.axhline(etp, ls="-", color="C2", lw=2,
                    label=rf"low-$q$ plateau $\eta={etp:.3f}\pm{erp:.3f}$")
    axR.axhline(0.0, color="0.8", lw=0.8)
    axR.axvspan(qmin, qmax, color="C1", alpha=0.12)
    axR.set_ylim(-0.5, 4.2)
    axR.set_xlabel(r"$q_r$"); axR.set_ylabel(r"$\eta_{eff}=4-d\ln G^{-1}/d\ln q$")
    axR.set_title("Running exponent (no angular averaging)")
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

    # Default window: from a few grid steps (finite-size floor) up to the
    # crossover q8~p8. The running-exponent panel lets you judge/tighten it.
    qmin = qmin_arg if qmin_arg is not None else 3 * a
    qmax = qmax_arg if qmax_arg is not None else p8

    # --- all modes as points, NO angular averaging -------------------------
    data = np.loadtxt(datfile, comments="#")
    qx, qy, qm, Gc, Gec = (data[:, 2], data[:, 3], data[:, 4],
                           data[:, 5], data[:, 6])
    m = (qm > 0) & (Gc > 0) & np.isfinite(Gc)
    q, theta = qm[m], np.arctan2(qy[m], qx[m])
    Ginv_pt = 1.0 / Gc[m]
    Ginv_err_pt = np.where(Gec[m] > 0, Gec[m] / Gc[m] ** 2, 0.0)
    err_arg = Ginv_err_pt if np.all(Ginv_err_pt > 0) else None

    # primary: unbinned weighted log-log fit over the window (cnt=1 per mode)
    cnt1 = np.ones_like(q)
    eta_w, err_w, npts = fit_eta_window(q, Ginv_pt, cnt1, qmin, qmax, err_arg)

    # display fit line over the window
    fitline = None
    mw = (q >= qmin) & (q <= qmax) & (Ginv_pt > 0)
    if mw.sum() >= 3:
        sl, ic = np.polyfit(np.log(q[mw]), np.log(Ginv_pt[mw]), 1)
        qq = np.logspace(np.log10(q[mw].min()), np.log10(q[mw].max()), 60)
        fitline = (qq, np.exp(ic) * qq ** sl)

    # running local exponent (sliding fit over the raw cloud) + low-q plateau
    qcen, eta_eff = running_eta(q, Ginv_pt)
    plateau = plateau_lowq(qcen, eta_eff) if len(qcen) else None

    if not quiet:
        print(f"file           : {datfile}")
        print(f"N, L           : {header.get('N','?')}, {header.get('L','?')}")
        print(f"p8             : {p8}")
        print(f"window         : q_r in [{qmin:.3f}, {qmax:.3f}]  ({npts} modes, unbinned)")
        if plateau is not None:
            etp, erp, (qlo, qhi), npl = plateau
            print(f"low-q plateau  : {etp:.3f} +/- {erp:.3f}  eta_eff flat over "
                  f"q_r in [{qlo:.3f},{qhi:.3f}] ({npl} pts)  [PRIMARY: q->0 exponent]")
        if eta_w is not None:
            print(f"windowed eta   : {eta_w:.3f} +/- {err_w:.3f}  (all modes in [{qmin:.3f},{qmax:.3f}]; biased by crossover if qmax high)")
        print(f"Poisson ratio  : {header.get('nu','?')} +/- {header.get('nu_err','?')}")
        conv = header.get("converged", "?")
        print(f"run            : sweeps={header.get('sweeps','?')} "
              f"Delta2 rel.err={header.get('rel_err','?')} converged={conv}")

    # Mirror the input path under plots/, keyed by the data filename stem so
    # multiple configs (or example_data files) never collide:
    #   data/N60/p0.40/fixed800/therm300_nt16_it1_seed12345.dat
    #     -> plots/N60/p0.40/fixed800/therm300_nt16_it1_seed12345.png
    #   example_data/N=100.dat -> plots/example_data/N=100.png
    d = os.path.dirname(datfile)
    stem = os.path.splitext(os.path.basename(datfile))[0]
    if d == "data" or d.startswith("data/"):
        pdir = "plots" + d[len("data"):]
    elif d:
        pdir = os.path.join("plots", d)
    else:
        pdir = "plots"
    os.makedirs(pdir, exist_ok=True)
    png = png or os.path.join(pdir, stem + ".png")
    gp = os.path.join(pdir, stem + ".gp")
    plot(q, Ginv_pt, theta, qcen, eta_eff, eta_w, fitline, p8, qmin, qmax,
         png, gp, datfile, plateau=plateau)
    return dict(N=int(header.get("N", 0)), p8=p8, png=png,
                eta_window=eta_w, eta_window_err=err_w,
                eta_plateau=(plateau[0] if plateau else None))


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
    ap.add_argument("--all", metavar="GLOB", nargs="?", const="data/N*/p*/*/*.dat",
                    help="batch-analyze every matching file (default data/N*/p*/*/*.dat)")
    ap.add_argument("--combined", metavar="GLOB", nargs="?", const="data/N*/p*/*/*.dat",
                    help="pool all N at each p8 into one weighted log-log fit; "
                         "writes plots/combined/{p<p8>.png, eta_vs_p8.png, combined_eta.csv} "
                         "(default glob data/N*/p*/*/*.dat; pass e.g. 'example_data/N=*.dat')")
    ap.add_argument("--combined-p8", type=float, default=None,
                    help="combined fit for a single p8 only (implies --combined)")
    ap.add_argument("--legacy", metavar="GLOB", nargs="?",
                    const="legacy/example_data/N=*.dat",
                    help="combined multi-N fit on ORIGINAL legacy-format files "
                         "(default legacy/example_data/N=*.dat, the preserved backup); "
                         "reformatted example_data/ now works with --combined instead")
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
            plot_combined(q, gi, ge, Ns, eta, err,
                          f"plots/combined/p{args.combined_p8:.2f}.png",
                          float(q.min()), args.combined_p8)
        else:
            combined_all(pattern=args.combined, nbins=args.nbins)
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
