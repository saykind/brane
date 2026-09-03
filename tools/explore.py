#!/usr/bin/env python3
"""
explore.py -- Map how the *measured* exponent depends on lattice size N and
coupling p8, using the brane engine + the estimators in analyze.py.

Reminder (see docs/model.md, README): eta is a UNIVERSAL exponent and does not
depend on p8. What changes with N and p8 is how much of the anomalous scaling
window (3a << q << q8~p8) fits inside the lattice, hence how close the *measured
effective* exponent gets to the true eta ~ 0.78. This script visualizes exactly
that finite-size / window effect:

  * eta vs N at fixed p8  -> should rise toward ~0.78 as N grows (window widens
    downward as a=2pi/L shrinks).
  * eta vs p8 at fixed N  -> rises with p8 because the crossover q8~p8 moves up,
    exposing more anomalous regime -- a MEASUREMENT effect, not real p8 dependence.

We report the windowed slope (stable, with error bar) and flag whether
eta_eff(q) actually has a plateau (spread < 0.25); points without a plateau are
drawn hollow to signal "eta not really determined at this size".

Usage:
    uv run tools/explore.py                     # default grids, ~10-15 min
    uv run tools/explore.py --quick             # small/fast grids
"""
import argparse
import subprocess
import sys
import numpy as np

import analyze  # same directory


def run_sim(N, p8, nt, therm, sweeps, seed=12345):
    out = f"data/explore_N{N}_p{int(round(p8*100)):03d}.dat"
    subprocess.run(["./brane", f"N={N}", f"p8={p8}", f"nt={nt}",
                    f"therm={therm}", f"sweeps={sweeps}", f"seed={seed}",
                    f"out={out}"], check=True, capture_output=True, text=True)
    return out


def measure(datfile):
    qmag, G, Ginv, header = analyze.load(datfile)
    L = float(header["L"]); a = 2 * np.pi / L
    p8 = float(header["p8"])
    qr, Gr, Ginv_r, cnt = analyze.radial_average(qmag, G, 60)
    qmin, qmax = 3 * a, max(p8, 5 * a)
    eta, err, _ = analyze.fit_eta_window(qr, Ginv_r, cnt, qmin, qmax)
    _, spread, _ = analyze.plateau_eta(qr, Ginv_r, cnt, qmin, qmax)
    nu = float(header.get("nu", "nan"))
    return eta, err, spread, nu


def sweep(configs, nt, therm, sweeps):
    rows = []
    for N, p8 in configs:
        f = run_sim(N, p8, nt, therm, sweeps)
        eta, err, spread, nu = measure(f)
        rows.append(dict(N=N, p8=p8, eta=eta, err=err, spread=spread, nu=nu))
        print(f"  N={N:3d} p8={p8:.2f} -> eta={eta:.3f}+/-{err:.3f} "
              f"(spread {spread:.2f}{'' if spread<0.25 else ' NO plateau'})",
              flush=True)
    return rows


def plot(rows_N, rows_p8, p8_fixed, N_fixed, out_png):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib missing; skipping plot (uv sync).")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    def draw(ax, xs, rows, xlabel, title):
        xs = np.array(xs)
        eta = np.array([r["eta"] for r in rows])
        err = np.array([r["err"] for r in rows])
        reliable = np.array([r["err"] < 0.06 for r in rows])
        ax.axhline(0.78, ls="--", color="0.5", label=r"universal $\eta\approx0.78$")
        ax.axhline(0.0, color="0.85", lw=0.8)
        # filled = statistically precise fit, hollow = large error bar
        ax.errorbar(xs[reliable], eta[reliable], yerr=err[reliable], fmt="o",
                    color="C0", capsize=3, label="fit err < 0.06")
        ax.errorbar(xs[~reliable], eta[~reliable], yerr=err[~reliable], fmt="o",
                    mfc="white", color="C0", capsize=3, label="fit err > 0.06")
        ax.set_xlabel(xlabel); ax.set_ylabel(r"measured effective $\eta$")
        ax.set_title(title); ax.set_ylim(0, 1.2)
        ax.legend(frameon=False, fontsize=8); ax.grid(alpha=0.3)

    draw(ax1, [r["N"] for r in rows_N], rows_N, "N  (L = 2N+1)",
         rf"$\eta$ vs size at $p_8={p8_fixed}$")
    draw(ax2, [r["p8"] for r in rows_p8], rows_p8, r"$p_8$",
         rf"$\eta$ vs coupling at $N={N_fixed}$ (window effect)")
    fig.tight_layout(); fig.savefig(out_png, dpi=140)
    print(f"[plot] wrote {out_png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nt", type=int, default=12)
    ap.add_argument("--therm", type=int, default=60)
    ap.add_argument("--sweeps", type=int, default=80)
    ap.add_argument("--p8fix", type=float, default=0.4)
    ap.add_argument("--Nfix", type=int, default=40)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--png", default="data/explore.png")
    args = ap.parse_args()

    if args.quick:
        Ns = [16, 24, 32]
        p8s = [0.3, 0.5, 0.8]
        args.therm, args.sweeps = 40, 50
    else:
        Ns = [24, 30, 36, 42, 48, 54, 60]
        p8s = [0.2, 0.3, 0.4, 0.6, 0.8, 1.0]

    print(f"eta vs N at p8={args.p8fix} (therm={args.therm} sweeps={args.sweeps} nt={args.nt}):")
    rows_N = sweep([(N, args.p8fix) for N in Ns], args.nt, args.therm, args.sweeps)
    print(f"eta vs p8 at N={args.Nfix}:")
    rows_p8 = sweep([(args.Nfix, p8) for p8 in p8s], args.nt, args.therm, args.sweeps)

    # CSV
    with open("data/explore.csv", "w") as f:
        f.write("scan,N,p8,eta,err,spread,nu\n")
        for r in rows_N:
            f.write(f"N,{r['N']},{r['p8']},{r['eta']:.4f},{r['err']:.4f},{r['spread']:.4f},{r['nu']:.4f}\n")
        for r in rows_p8:
            f.write(f"p8,{r['N']},{r['p8']},{r['eta']:.4f},{r['err']:.4f},{r['spread']:.4f},{r['nu']:.4f}\n")
    print("[csv] wrote data/explore.csv")

    plot(rows_N, rows_p8, args.p8fix, args.Nfix, args.png)


if __name__ == "__main__":
    main()
