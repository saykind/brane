#!/usr/bin/env python3
"""
heatmap.py -- 2D map of the measured effective exponent eta over the (N, p8)
plane, to disentangle finite-size (N) from coupling/window (p8) effects.

Interpretation guide:
  * vertical contours (eta depends on p8, not N)  -> at these sizes the coupling
    / fit-window sets the measurement; growing N alone barely helps.
  * horizontal/diagonal contours (eta rises with N) -> genuine finite-size
    convergence toward the universal eta ~ 0.78.

eta is measured exactly as in analyze.py: rotationally average G(q_r), then a
weighted log-log slope over the window [3a, max(p8, 5a)].

Usage:
    uv run tools/heatmap.py                       # default grid (~12 min)
    uv run tools/heatmap.py --Ns 16,24,32 --p8s 0.3,0.5,0.8 --sweeps 60
"""
import argparse
import subprocess
import numpy as np

import analyze  # same directory


def refine_grid(Ns, p8s, eta, k):
    """Bilinearly upsample eta onto a k-times finer uniform grid (no new sims;
    inserts (k-1) interpolated points between each pair of real points)."""
    Ns = np.asarray(Ns, float); p8s = np.asarray(p8s, float)
    fp = np.linspace(p8s[0], p8s[-1], (len(p8s) - 1) * k + 1)
    fN = np.linspace(Ns[0], Ns[-1], (len(Ns) - 1) * k + 1)
    tmp = np.vstack([np.interp(fp, p8s, eta[i, :]) for i in range(len(Ns))])
    fine = np.vstack([np.interp(fN, Ns, tmp[:, j]) for j in range(len(fp))]).T
    return fN, fp, fine


def load_csv(path):
    """Read a heatmap.csv written by this tool -> (Ns, p8s, eta 2d)."""
    with open(path) as f:
        head = f.readline().strip().split(",")
        p8s = [float(h.split("=")[1]) for h in head[1:]]
        Ns, rows = [], []
        for line in f:
            parts = line.strip().split(",")
            Ns.append(int(parts[0]))
            rows.append([float(x) for x in parts[1:]])
    return Ns, p8s, np.array(rows)


def plot_map(Ns, p8s, eta, png, refine=1):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib missing; CSV only.")
        return
    Ns_a, p8s_a = np.array(Ns, float), np.array(p8s, float)

    if refine > 1:
        fN, fp, fine = refine_grid(Ns_a, p8s_a, eta, refine)
    else:
        fN, fp, fine = Ns_a, p8s_a, eta

    fig, ax = plt.subplots(figsize=(8, 6))
    ext = [fp[0], fp[-1], fN[0], fN[-1]]
    im = ax.imshow(fine, origin="lower", extent=ext, aspect="auto",
                   cmap="viridis", vmin=0.2, vmax=0.9,
                   interpolation="bilinear")
    fig.colorbar(im, ax=ax, label=r"effective $\eta$ (bilinear-interpolated)"
                 if refine > 1 else r"measured effective $\eta$")
    P, NN = np.meshgrid(fp, fN)
    try:
        cs = ax.contour(P, NN, fine, levels=[0.5, 0.6, 0.7, 0.78], colors="white",
                        linewidths=[1, 1, 1, 2], linestyles=["--", "--", "--", "-"])
        ax.clabel(cs, fmt="%.2f", fontsize=8)
    except Exception:
        pass
    # mark the REAL computed points and annotate them
    for i, N in enumerate(Ns):
        for j, p8 in enumerate(p8s):
            if np.isfinite(eta[i, j]):
                ax.plot(p8, N, "o", ms=3, color="white", mec="black", mew=0.5)
                ax.text(p8, N + (Ns_a[-1] - Ns_a[0]) * 0.012, f"{eta[i,j]:.2f}",
                        ha="center", va="bottom", fontsize=6, color="white")
    ax.set_xlabel(r"$p_8$  (coupling / crossover $q_8\sim p_8$)")
    ax.set_ylabel(r"$N$  (size, $L=2N+1$)")
    tag = f"  [bilinear x{refine}]" if refine > 1 else ""
    ax.set_title(rf"Effective $\eta(N, p_8)${tag}  (dots = computed; white line = 0.78)")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print(f"[plot] wrote {png}")


def run_and_measure(N, p8, nt, therm, sweeps):
    out = f"data/hm_N{N}_p{int(round(p8*100)):03d}.dat"
    subprocess.run(["./brane", f"N={N}", f"p8={p8}", f"nt={nt}",
                    f"therm={therm}", f"sweeps={sweeps}", f"out={out}"],
                   check=True, capture_output=True, text=True)
    return measure_file(out)


def measure_file(path):
    """Return (N, p8, eta) measured from a single brane output file."""
    qmag, G, Ginv, header = analyze.load(path)
    N = int(header["N"]); L = float(header["L"]); a = 2 * np.pi / L
    p8 = float(header["p8"])
    qr, Gr, Ginv_r, cnt = analyze.radial_average(qmag, G, 60)
    eta, err, _ = analyze.fit_eta_window(qr, Ginv_r, cnt, 3 * a, max(p8, 5 * a))
    return N, p8, eta


def collect_all(pattern="data/hm_*.dat"):
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


def scattered_to_grid(pts):
    """Reshape scattered (N,p8,eta) points onto a regular grid (NaN if missing)."""
    Ns = sorted(set(p[0] for p in pts))
    p8s = sorted(set(round(p[1], 3) for p in pts))
    eta = np.full((len(Ns), len(p8s)), np.nan)
    iN = {n: i for i, n in enumerate(Ns)}
    iP = {p: j for j, p in enumerate(p8s)}
    for N, p8, e in pts:
        eta[iN[N], iP[round(p8, 3)]] = e
    return Ns, p8s, eta


def plot_scattered(pts, png):
    """Colormap from an arbitrary (possibly irregular) set of (N,p8,eta) points
    via Delaunay triangulation -- lets us combine grids from multiple runs."""
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
    fig, ax = plt.subplots(figsize=(8.5, 6))
    tcf = ax.tricontourf(triang, E, levels=np.linspace(0.35, 0.9, 23),
                         cmap="viridis", extend="both")
    fig.colorbar(tcf, ax=ax, label=r"effective $\eta$")
    try:
        cs = ax.tricontour(triang, E, levels=[0.5, 0.6, 0.7, 0.78],
                           colors="white", linewidths=[1, 1, 1, 2],
                           linestyles=["--", "--", "--", "-"])
        ax.clabel(cs, fmt="%.2f", fontsize=8)
    except Exception:
        pass
    ax.plot(P, NN, "o", ms=3, color="white", mec="black", mew=0.4)
    ax.set_xlabel(r"$p_8$  (coupling / crossover $q_8\sim p_8$)")
    ax.set_ylabel(r"$N$  (size, $L=2N+1$)")
    ax.set_title(rf"Effective $\eta(N, p_8)$ -- {len(pts)} runs combined "
                 r"(white line = 0.78)")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print(f"[plot] wrote {png}  ({len(pts)} cells)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ns", default="16,24,32,40,48")
    ap.add_argument("--p8s", default="0.25,0.4,0.55,0.7,0.85,1.0")
    ap.add_argument("--nt", type=int, default=12)
    ap.add_argument("--therm", type=int, default=40)
    ap.add_argument("--sweeps", type=int, default=60)
    ap.add_argument("--png", default="data/heatmap.png")
    ap.add_argument("--from-csv", default=None,
                    help="replot an existing heatmap.csv instead of simulating")
    ap.add_argument("--replot-all", action="store_true",
                    help="replot EVERY data/hm_*.dat cell on disk (combine all runs)")
    ap.add_argument("--refine", type=int, default=1,
                    help="bilinear upsample factor for a smoother map (no new sims)")
    args = ap.parse_args()

    if args.replot_all:
        pts = collect_all()
        with open("data/heatmap_all.csv", "w") as f:
            f.write("N,p8,eta\n")
            for N, p8, e in sorted(pts):
                f.write(f"{N},{p8:.3f},{e:.4f}\n")
        print(f"[csv] wrote data/heatmap_all.csv ({len(pts)} cells)")
        if args.refine > 1:
            Ns, p8s, eta = scattered_to_grid(pts)
            if np.isnan(eta).any():
                print("  grid incomplete -> triangulation plot (no bilinear refine)")
                plot_scattered(pts, args.png)
            else:
                plot_map(Ns, p8s, eta, args.png, refine=args.refine)
        else:
            plot_scattered(pts, args.png)
        return

    if args.from_csv:
        Ns, p8s, eta = load_csv(args.from_csv)
        plot_map(Ns, p8s, eta, args.png, refine=args.refine)
        return

    Ns = [int(x) for x in args.Ns.split(",")]
    p8s = [float(x) for x in args.p8s.split(",")]
    eta = np.full((len(Ns), len(p8s)), np.nan)

    for i, N in enumerate(Ns):
        for j, p8 in enumerate(p8s):
            _, _, e = run_and_measure(N, p8, args.nt, args.therm, args.sweeps)
            eta[i, j] = e if e is not None else np.nan
            print(f"  N={N:3d} p8={p8:.2f} -> eta={e:.3f}", flush=True)

    with open("data/heatmap.csv", "w") as f:
        f.write("N," + ",".join(f"p8={p:.2f}" for p in p8s) + "\n")
        for i, N in enumerate(Ns):
            f.write(f"{N}," + ",".join(f"{eta[i,j]:.4f}" for j in range(len(p8s))) + "\n")
    print("[csv] wrote data/heatmap.csv")

    plot_map(Ns, p8s, eta, args.png, refine=args.refine)


if __name__ == "__main__":
    main()
