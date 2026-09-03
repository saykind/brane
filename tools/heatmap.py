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


def run_and_measure(N, p8, nt, therm, sweeps):
    out = f"data/hm_N{N}_p{int(round(p8*100)):03d}.dat"
    subprocess.run(["./brane", f"N={N}", f"p8={p8}", f"nt={nt}",
                    f"therm={therm}", f"sweeps={sweeps}", f"out={out}"],
                   check=True, capture_output=True, text=True)
    qmag, G, Ginv, header = analyze.load(out)
    L = float(header["L"]); a = 2 * np.pi / L
    qr, Gr, Ginv_r, cnt = analyze.radial_average(qmag, G, 60)
    eta, err, _ = analyze.fit_eta_window(qr, Ginv_r, cnt, 3 * a, max(p8, 5 * a))
    return eta, err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ns", default="16,24,32,40,48")
    ap.add_argument("--p8s", default="0.25,0.4,0.55,0.7,0.85,1.0")
    ap.add_argument("--nt", type=int, default=12)
    ap.add_argument("--therm", type=int, default=40)
    ap.add_argument("--sweeps", type=int, default=60)
    ap.add_argument("--png", default="data/heatmap.png")
    args = ap.parse_args()

    Ns = [int(x) for x in args.Ns.split(",")]
    p8s = [float(x) for x in args.p8s.split(",")]
    eta = np.full((len(Ns), len(p8s)), np.nan)
    err = np.full_like(eta, np.nan)

    for i, N in enumerate(Ns):
        for j, p8 in enumerate(p8s):
            e, de = run_and_measure(N, p8, args.nt, args.therm, args.sweeps)
            eta[i, j] = e if e is not None else np.nan
            err[i, j] = de if de is not None else np.nan
            print(f"  N={N:3d} p8={p8:.2f} -> eta={e:.3f}", flush=True)

    # CSV
    with open("data/heatmap.csv", "w") as f:
        f.write("N," + ",".join(f"p8={p:.2f}" for p in p8s) + "\n")
        for i, N in enumerate(Ns):
            f.write(f"{N}," + ",".join(f"{eta[i,j]:.4f}" for j in range(len(p8s))) + "\n")
    print("[csv] wrote data/heatmap.csv")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib missing; CSV only.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    # cell-centered pcolormesh
    p8s_a, Ns_a = np.array(p8s), np.array(Ns)
    def edges(v):
        v = np.asarray(v, float)
        e = np.concatenate([[v[0] - (v[1]-v[0])/2],
                            (v[:-1]+v[1:])/2,
                            [v[-1] + (v[-1]-v[-2])/2]])
        return e
    pm = ax.pcolormesh(edges(p8s_a), edges(Ns_a), eta, cmap="viridis",
                       vmin=0.2, vmax=0.9, shading="flat")
    cb = fig.colorbar(pm, ax=ax, label=r"measured effective $\eta$")
    # contour at the universal value
    P, Nn = np.meshgrid(p8s_a, Ns_a)
    try:
        cs = ax.contour(P, Nn, eta, levels=[0.6, 0.7, 0.78], colors="white",
                        linewidths=[1, 1, 2], linestyles=["--", "--", "-"])
        ax.clabel(cs, fmt="%.2f", fontsize=8)
    except Exception:
        pass
    # annotate cells
    for i, N in enumerate(Ns):
        for j, p8 in enumerate(p8s):
            if np.isfinite(eta[i, j]):
                ax.text(p8, N, f"{eta[i,j]:.2f}", ha="center", va="center",
                        fontsize=7, color="w")
    ax.set_xlabel(r"$p_8$  (coupling / crossover $q_8\sim p_8$)")
    ax.set_ylabel(r"$N$  (size, $L=2N+1$)")
    ax.set_title(r"Effective $\eta(N, p_8)$  (white solid = universal 0.78)")
    fig.tight_layout(); fig.savefig(args.png, dpi=140)
    print(f"[plot] wrote {args.png}")


if __name__ == "__main__":
    main()
