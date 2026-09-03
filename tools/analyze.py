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

Only numpy is required (polyfit).  matplotlib is used for plots if present;
otherwise a gnuplot script is emitted as a fallback.

Usage
-----
    python3 tools/analyze.py data/N=40.dat
    python3 tools/analyze.py data/N=40.dat --qmin 0.23 --qmax 0.35 --nbins 40
"""
import argparse
import sys
import numpy as np


def load(path):
    """Return (q1, q2, qmag, G, Ginv) and the header dict."""
    header = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                # parse "key=value" tokens from the first header line
                for tok in line[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        header[k] = v
                continue
            break
    data = np.loadtxt(path, comments="#")
    q1, q2 = data[:, 0], data[:, 1]
    qmag = data[:, 4]
    G = data[:, 5]
    Ginv = data[:, 6]
    return q1, q2, qmag, G, Ginv, header


def radial_bin(qmag, y, nbins):
    """Average y in log-spaced |q| bins; return (qcenter, ymean)."""
    mask = (qmag > 0) & (y > 0) & np.isfinite(y)
    q, y = qmag[mask], y[mask]
    edges = np.logspace(np.log10(q.min()), np.log10(q.max()), nbins + 1)
    idx = np.digitize(q, edges)
    qc, ym = [], []
    for b in range(1, nbins + 1):
        sel = idx == b
        if sel.sum() >= 1:
            qc.append(q[sel].mean())
            ym.append(y[sel].mean())
    return np.array(qc), np.array(ym)


def fit_eta(qmag, Ginv, qmin, qmax):
    """Linear fit of log Ginv vs log q in [qmin, qmax]; return eta, err, npts."""
    mask = (qmag >= qmin) & (qmag <= qmax) & (Ginv > 0) & np.isfinite(Ginv)
    if mask.sum() < 3:
        return None, None, int(mask.sum())
    lx, ly = np.log(qmag[mask]), np.log(Ginv[mask])
    coef, cov = np.polyfit(lx, ly, 1, cov=True)
    slope = coef[0]
    slope_err = float(np.sqrt(cov[0, 0]))
    return 4.0 - slope, slope_err, int(mask.sum())


def plot(qc, Ginvc, qmag, Ginv, eta, p8, qmin, qmax, out_png, out_gp, datfile):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        # gnuplot fallback
        with open(out_gp, "w") as f:
            f.write(f"""set terminal pngcairo size 900,650
set output '{out_png}'
set logscale xy
set grid
set xlabel 'q'; set ylabel 'G^{{-1}}(q)'
set title 'Inverse Green function (eta={eta:.3f}, p8={p8})'
plot '{datfile}' using 5:7 pt 7 ps 0.4 title 'MC', \\
     x**4 lw 2 title 'q^4', \\
     x**(4-{eta}) lw 2 title 'q^{{4-eta}}'
""")
        print(f"[plot] matplotlib not found; wrote gnuplot script {out_gp}")
        print(f"       run:  gnuplot {out_gp}   (or: pip3 install matplotlib)")
        return
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.loglog(qmag, Ginv, ".", ms=2, alpha=0.25, color="0.6", label="all modes")
    ax.loglog(qc, Ginvc, "o", ms=5, color="C0", label="radial mean")
    qq = np.logspace(np.log10(qmag[qmag > 0].min()),
                     np.log10(qmag.max()), 100)
    ax.loglog(qq, qq**4, "--", color="C3", lw=1.5, label=r"$q^4$ (harmonic)")
    ax.loglog(qq, qq**(4 - eta), "-", color="C2", lw=2,
              label=rf"$q^{{4-\eta}},\ \eta={eta:.2f}$")
    ax.axvspan(qmin, qmax, color="C1", alpha=0.12, label="fit window")
    ax.set_xlabel("q")
    ax.set_ylabel(r"$G^{-1}(q)$")
    ax.set_title(rf"Inverse Green function ($p_8={p8}$, $\eta={eta:.3f}$)")
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
    ap.add_argument("--nbins", type=int, default=40)
    ap.add_argument("--png", default=None)
    args = ap.parse_args()

    q1, q2, qmag, G, Ginv, header = load(args.datfile)
    p8 = float(header.get("p8", 0.3))
    a = 2 * np.pi / float(header.get("L", 81))

    # default fit window: above finite-size (3rd mode) up to the crossover ~p8
    qmin = args.qmin if args.qmin is not None else 3 * a
    qmax = args.qmax if args.qmax is not None else max(p8, 5 * a)

    eta, err, npts = fit_eta(qmag, Ginv, qmin, qmax)
    qc, Ginvc = radial_bin(qmag, Ginv, args.nbins)

    print(f"file           : {args.datfile}")
    print(f"N, L           : {header.get('N','?')}, {header.get('L','?')}")
    print(f"p8             : {p8}")
    print(f"samples        : {header.get('samples','?')}")
    print(f"fit window     : q in [{qmin:.4f}, {qmax:.4f}]  ({npts} modes)")
    if eta is not None:
        print(f"eta            : {eta:.3f} +/- {err:.3f}   (slope err)")
    else:
        print("eta            : (too few modes in window; widen --qmin/--qmax)")
    print(f"Poisson ratio  : {header.get('nu','?')}  (from simulation)")

    png = args.png or (args.datfile.rsplit(".", 1)[0] + ".png")
    gp = args.datfile.rsplit(".", 1)[0] + ".gp"
    if eta is not None:
        plot(qc, Ginvc, qmag, Ginv, eta, p8, qmin, qmax, png, gp, args.datfile)


if __name__ == "__main__":
    main()
