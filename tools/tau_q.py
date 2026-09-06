#!/usr/bin/env python3
"""tau_q.py -- per-mode autocorrelation time tau(q) from qseries files.

The engine's `qseries=<path>` option records the per-sweep |h_q|^2 for the
qx-axis ray q=(j,0), j=1..N (replica 0). This tool computes the integrated
autocorrelation time tau for each mode (each column) via Sokal windowing (shared
with autocorr.py) and plots tau vs |q|.

If several qseries files are given (independent seeds), tau is computed per file
per mode and averaged, with the spread shown as an error band -- this is the
clean way to get tau(q), since a single chain's tau estimate is noisy.

The physical question: is tau(q) roughly FLAT across |q| (good -- no critical
slowing down, the Troster OFMC goal), or does it blow up as q->0 (bad)?

Usage:
    uv run tools/tau_q.py <qseries...> [--burn 100] [--out tau_q.png]
"""
import argparse
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from autocorr import integrated_tau


def read_qmag(path):
    with open(path) as f:
        for line in f:
            if line.startswith("# qmag:"):
                return np.array([float(x) for x in line.split(":", 1)[1].split()])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("qseries", nargs="+")
    ap.add_argument("--burn", type=int, default=100, help="drop first N sweeps")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    qmag = read_qmag(a.qseries[0])
    if qmag is None:
        print("no '# qmag:' header found"); sys.exit(1)
    nmodes = len(qmag)

    # tau[file, mode]
    taus = np.full((len(a.qseries), nmodes), np.nan)
    for fi, path in enumerate(a.qseries):
        d = np.loadtxt(path, comments="#")
        cols = d[:, 1:]                       # drop the sweep column
        if a.burn:
            cols = cols[a.burn:]
        for m in range(min(nmodes, cols.shape[1])):
            t, _, _ = integrated_tau(cols[:, m])
            taus[fi, m] = t

    tau_mean = np.nanmean(taus, axis=0)
    tau_std = np.nanstd(taus, axis=0, ddof=1) if len(a.qseries) > 1 else np.zeros(nmodes)
    nseed = len(a.qseries)
    tau_sem = tau_std / np.sqrt(nseed) if nseed > 1 else np.zeros(nmodes)

    # crossover q_c ~ p8 (read from header if present)
    p8 = None
    with open(a.qseries[0]) as f:
        m = re.search(r"p8=([\d.]+)", f.read(400))
        if m:
            p8 = float(m.group(1))

    print(f"tau(q) from {nseed} seed(s), burn={a.burn}")
    print(f"{'|q|':>8} {'tau':>8} {'+/-':>7}")
    for m in range(nmodes):
        print(f"{qmag[m]:>8.4f} {tau_mean[m]:>8.2f} {tau_sem[m]:>7.2f}")

    out = a.out or (os.path.splitext(a.qseries[0])[0] + ".tau_q.png")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.errorbar(qmag, tau_mean, yerr=tau_sem, fmt="-o", ms=4, capsize=2,
                color="tab:blue", label=f"τ(q), {nseed} seed(s)")
    if p8:
        ax.axvline(p8, color="tab:purple", ls=":", lw=1.3,
                   label=f"crossover q$_c\\approx$p8={p8:g}")
    ax.set_xscale("log")
    ax.set_xlabel("|q|")
    ax.set_ylabel("integrated autocorrelation time τ (sweeps)")
    ax.set_title("per-mode autocorrelation time τ(q)")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
