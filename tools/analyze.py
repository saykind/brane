#!/usr/bin/env python3
"""
analyze.py -- Extract the anomalous elasticity exponent eta and (echo) the
Poisson ratio from a brane Fourier-Monte-Carlo output file, and plot the
inverse Green function.

Physics
-------
The height-height correlator (Green function) is  G(q) = <|h_q|^2>.  In the
harmonic theory G(q) ~ q^-4, i.e. the inverse Green G^-1 ~ q^4.  Anomalous
elasticity renormalizes this in the scaling window to

        G^-1(q) ~ q^(4 - eta)  ,      q << p8

so a straight-line fit of log G^-1 vs log q has slope (4 - eta):

        eta = 4 - slope .

Rotational averaging
--------------------
In the continuum the theory is isotropic: G depends only on q_r = |q|, not on
the direction q_theta.  The lattice breaks this weakly at large q, but in the
small-q scaling window it is an excellent symmetry.  Rather than fitting only
the modes that happen to lie on the x/y axes and the diagonals (the original
thesis approach), we therefore map *every* mode (q_x, q_y) to its radius q_r,
average G over each thin annulus in q_r (the theta average), and fit the
resulting radially-averaged G(q_r).  This uses all L^2 modes and every
direction, so it is far less noisy.

Only numpy is required (polyfit).  matplotlib is used for plots if present;
otherwise a gnuplot script is emitted as a fallback.

Usage
-----
    python3 tools/analyze.py data/N=40.dat
    python3 tools/analyze.py data/N=40.dat --qmin 0.23 --qmax 0.40 --nbins 60
"""
import argparse
import sys
import numpy as np


def load(path):
    """Return (q1, q2, qmag, G, Ginv) and the header dict."""
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
            "  It looks like the legacy multi-line .dat format. The legacy\n"
            "  binary (legacy/a.out) overwrites data/N=<N>.dat in that format.\n"
            "  Regenerate with the new engine, e.g.:\n"
            f"      ./brane N=40 p8=0.4 out={path}\n"
            "  then re-run this analysis."
        )
    data = np.loadtxt(path, comments="#")
    q1, q2 = data[:, 0], data[:, 1]
    qmag = data[:, 4]
    G = data[:, 5]
    Ginv = data[:, 6]
    return q1, q2, qmag, G, Ginv, header


def radial_average(qmag, G, nbins):
    """Rotationally average G over log-spaced |q| shells.

    Returns (q_r, G_r, Ginv_r, count, sem) where G_r is the mean correlator in
    each annulus, Ginv_r = 1/G_r, and sem is the standard error of G_r.
    """
    mask = (qmag > 0) & (G > 0) & np.isfinite(G)
    q, g = qmag[mask], G[mask]
    edges = np.logspace(np.log10(q.min()), np.log10(q.max()), nbins + 1)
    idx = np.digitize(q, edges)
    qr, gr, cnt, sem = [], [], [], []
    for b in range(1, nbins + 1):
        sel = idx == b
        c = int(sel.sum())
        if c >= 1:
            gvals = g[sel]
            qr.append(q[sel].mean())
            gr.append(gvals.mean())
            cnt.append(c)
            sem.append(gvals.std(ddof=1) / np.sqrt(c) if c > 1 else 0.0)
    qr, gr, cnt, sem = map(np.array, (qr, gr, cnt, sem))
    return qr, gr, 1.0 / gr, cnt, sem


def fit_eta(qr, Ginv_r, cnt, qmin, qmax, weighted=True):
    """Fit log Ginv_r vs log q_r over [qmin, qmax] on the radial averages.

    Returns (eta, err, nshells). Each q_r shell contributes once; with
    weighted=True shells are weighted by sqrt(count) (reliability of the
    theta-average).
    """
    m = (qr >= qmin) & (qr <= qmax) & (Ginv_r > 0) & np.isfinite(Ginv_r)
    if m.sum() < 3:
        return None, None, int(m.sum())
    lx, ly = np.log(qr[m]), np.log(Ginv_r[m])
    w = np.sqrt(cnt[m]) if weighted else None
    coef, cov = np.polyfit(lx, ly, 1, w=w, cov=True)
    slope = coef[0]
    return 4.0 - slope, float(np.sqrt(cov[0, 0])), int(m.sum())


def plot(qr, Ginv_r, cnt, qmag, Ginv, eta, p8, qmin, qmax, out_png, out_gp, datfile):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        with open(out_gp, "w") as f:
            f.write(f"""set terminal pngcairo size 900,650
set output '{out_png}'
set logscale xy
set grid
set xlabel 'q_r'; set ylabel 'G^{{-1}}(q_r)'
set title 'Inverse Green function (eta={eta:.3f}, p8={p8})'
plot '{datfile}' using 5:7 pt 7 ps 0.3 lc rgb '#cccccc' title 'all modes', \\
     x**4 lw 2 title 'q^4', \\
     x**(4-{eta}) lw 2 title 'q^{{4-eta}}'
""")
        print(f"[plot] matplotlib not found; wrote gnuplot script {out_gp}")
        print(f"       run:  gnuplot {out_gp}   (or: uv sync)")
        return
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.loglog(qmag, Ginv, ".", ms=2, alpha=0.20, color="0.7",
              label="all modes (q_x, q_y)")
    # size radial points by how many modes were averaged in the shell
    sizes = 12 + 40 * (cnt / cnt.max())
    ax.scatter(qr, Ginv_r, s=sizes, color="C0", zorder=3,
               label=r"radial average $G^{-1}(q_r)$")
    qq = np.logspace(np.log10(qmag[qmag > 0].min()), np.log10(qmag.max()), 100)
    ax.loglog(qq, qq**4, "--", color="C3", lw=1.5, label=r"$q^4$ (harmonic)")
    ax.loglog(qq, qq**(4 - eta), "-", color="C2", lw=2,
              label=rf"$q^{{4-\eta}},\ \eta={eta:.2f}$")
    ax.axvspan(qmin, qmax, color="C1", alpha=0.12, label="fit window")
    ax.set_xlabel(r"$q_r = |q|$")
    ax.set_ylabel(r"$G^{-1}(q_r)$")
    ax.set_title(rf"Radially-averaged inverse Green ($p_8={p8}$, $\eta={eta:.3f}$)")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    print(f"[plot] wrote {out_png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datfile")
    ap.add_argument("--qmin", type=float, default=None)
    ap.add_argument("--qmax", type=float, default=None)
    ap.add_argument("--nbins", type=int, default=60, help="radial shells")
    ap.add_argument("--unweighted", action="store_true",
                    help="equal weight per shell (default weights by sqrt(count))")
    ap.add_argument("--png", default=None)
    args = ap.parse_args()

    q1, q2, qmag, G, Ginv, header = load(args.datfile)
    p8 = float(header.get("p8", 0.3))
    a = 2 * np.pi / float(header.get("L", 81))

    # default fit window: above finite-size (3rd mode) up to the crossover ~p8
    qmin = args.qmin if args.qmin is not None else 3 * a
    qmax = args.qmax if args.qmax is not None else max(p8, 5 * a)

    qr, Gr, Ginv_r, cnt, sem = radial_average(qmag, G, args.nbins)
    eta, err, nshells = fit_eta(qr, Ginv_r, cnt, qmin, qmax,
                                weighted=not args.unweighted)

    print(f"file           : {args.datfile}")
    print(f"N, L           : {header.get('N','?')}, {header.get('L','?')}")
    print(f"p8             : {p8}")
    print(f"samples        : {header.get('samples','?')}")
    print(f"radial shells  : {len(qr)} (all {len(qmag)} modes averaged by |q|)")
    print(f"fit window     : q_r in [{qmin:.4f}, {qmax:.4f}]  ({nshells} shells)")
    if eta is not None:
        wtag = "unweighted" if args.unweighted else "weighted by sqrt(count)"
        print(f"eta            : {eta:.3f} +/- {err:.3f}   ({wtag})")
    else:
        print("eta            : (too few shells in window; widen --qmin/--qmax)")
    print(f"Poisson ratio  : {header.get('nu','?')}  (from simulation)")

    png = args.png or (args.datfile.rsplit(".", 1)[0] + ".png")
    gp = args.datfile.rsplit(".", 1)[0] + ".gp"
    if eta is not None:
        plot(qr, Ginv_r, cnt, qmag, Ginv, eta, p8, qmin, qmax, png, gp, args.datfile)


if __name__ == "__main__":
    main()
