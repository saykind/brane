#!/usr/bin/env python3
"""
heatmap.py -- 2D map of the measured effective exponent eta over the (N, p8)
plane, to disentangle finite-size (N) from coupling/window (p8) effects.

Interpretation guide:
  * vertical contours (eta depends on p8, not N)  -> at these sizes the coupling
    / fit-window sets the measurement; growing N alone barely helps.
  * horizontal/diagonal contours (eta rises with N) -> genuine finite-size
    convergence toward the universal eta ~ 0.78.

Each (N, p8) cell is saved as data/N<N>/p<p8>/data.dat, so runs accumulate on
disk; --replot-all combines every cell present. eta is measured exactly as in
analyze.py: rotationally average G(q_r), then a weighted log-log slope over the
window [3a, p8].

Usage:
    uv run tools/heatmap.py                                  # default grid
    uv run tools/heatmap.py --Ns 16,24,32 --p8s 0.3,0.5,0.8 --sweeps 60
    uv run tools/heatmap.py --replot-all                     # combine all cells
    uv run tools/heatmap.py --replot-all --refine 2          # smoother map
"""
import argparse
import os
import subprocess
import numpy as np

import analyze  # same directory


def measure_file(path):
    """Return (N, p8, eta) measured from a single brane output file."""
    qmag, G, Gerr, Ginv, header = analyze.load(path)
    N = int(header["N"]); L = float(header["L"]); a = 2 * np.pi / L
    p8 = float(header["p8"])
    qr, Gr, Ginv_r, cnt, Ginv_err = analyze.radial_average(qmag, G, 60, Gerr)
    eta, err, _ = analyze.fit_eta_window(qr, Ginv_r, cnt, 3 * a, p8, Ginv_err)
    return N, p8, eta


def run_and_measure(N, p8, nt, therm, sweeps, eps):
    out = f"data/N{N}/p{p8:.2f}/data.dat"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    subprocess.run(["./brane", f"N={N}", f"p8={p8}", f"nt={nt}",
                    f"therm={therm}", f"sweeps={sweeps}", f"eps={eps}",
                    f"out={out}"],
                   check=True, capture_output=True, text=True)
    return measure_file(out)


def collect_all(pattern="data/N*/p*/*/*.dat"):
    """Gather (N, p8, eta) from every cell file on disk (combines all runs)."""
    import glob
    pts = []
    for f in sorted(glob.glob(pattern)):
        try:
            N, p8, eta = measure_file(f)
            if eta is not None and np.isfinite(eta):
                pts.append((N, p8, eta))
        except Exception as e:
            print(f"  skip {f}: {e}")
    return pts


def plot_eta(pts, png, refine=0):
    """The single heatmap style: filled contours over a Delaunay triangulation
    of the (N, p8, eta) points, white contour lines, and the computed points
    marked as dots. Works for any (regular or irregular) combined grid.

    refine>0 subdivides the triangulation (cubic interpolation) for a smoother
    map -- same style, just finer -- without any new simulations.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ImportError:
        print("matplotlib missing; cannot plot.")
        return
    P = np.array([p[1] for p in pts], float)
    NN = np.array([p[0] for p in pts], float)
    E = np.array([p[2] for p in pts], float)

    triang = mtri.Triangulation(P, NN)
    if refine and refine > 0:
        plot_tri, plot_E = mtri.UniformTriRefiner(triang).refine_field(
            E, subdiv=refine)
    else:
        plot_tri, plot_E = triang, E

    fig, ax = plt.subplots(figsize=(8.5, 6))
    tcf = ax.tricontourf(plot_tri, plot_E, levels=np.linspace(0.15, 0.9, 26),
                         cmap="viridis", extend="both")
    fig.colorbar(tcf, ax=ax, label=r"effective $\eta$")
    try:
        cs = ax.tricontour(plot_tri, plot_E, levels=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
                           colors="white", linewidths=1.0, linestyles="--")
        ax.clabel(cs, fmt="%.2f", fontsize=8)
    except Exception:
        pass
    ax.plot(P, NN, "o", ms=3, color="white", mec="black", mew=0.4)
    ax.set_xlabel(r"$p_8$  (coupling / crossover $q_8\sim p_8$)")
    ax.set_ylabel(r"$N$  (size, $L=2N+1$)")
    tag = f"  [x{refine} refined]" if refine else ""
    ax.set_title(rf"Effective $\eta(N, p_8)$ -- {len(pts)} runs combined{tag}")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print(f"[plot] wrote {png}  ({len(pts)} cells)")


def write_csv(pts, path="plots/heatmap_all.csv"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("N,p8,eta\n")
        for N, p8, e in sorted(pts):
            f.write(f"{N},{p8:.3f},{e:.4f}\n")
    print(f"[csv] wrote {path} ({len(pts)} cells)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ns", default="16,24,32,40,48")
    ap.add_argument("--p8s", default="0.25,0.4,0.55,0.7,0.85,1.0")
    ap.add_argument("--nt", type=int, default=12)
    ap.add_argument("--therm", type=int, default=40)
    ap.add_argument("--sweeps", type=int, default=60)
    ap.add_argument("--eps", type=float, default=0.01,
                    help="convergence target passed to brane (rel err on Delta2)")
    ap.add_argument("--png", default="plots/heatmap.png")
    ap.add_argument("--replot-all", action="store_true",
                    help="replot EVERY data/N*/p*/*/*.dat cell on disk (combine all runs)")
    ap.add_argument("--refine", type=int, default=0,
                    help="triangulation subdivisions for a smoother map (no new sims)")
    args = ap.parse_args()
    os.makedirs("plots", exist_ok=True)

    if not args.replot_all:
        Ns = [int(x) for x in args.Ns.split(",")]
        p8s = [float(x) for x in args.p8s.split(",")]
        for N in Ns:
            for p8 in p8s:
                _, _, e = run_and_measure(N, p8, args.nt, args.therm,
                                          args.sweeps, args.eps)
                estr = f"{e:.3f}" if e is not None else "None (window too narrow)"
                print(f"  N={N:3d} p8={p8:.2f} -> eta={estr}", flush=True)

    pts = collect_all()
    write_csv(pts)
    plot_eta(pts, args.png, refine=args.refine)


if __name__ == "__main__":
    main()
