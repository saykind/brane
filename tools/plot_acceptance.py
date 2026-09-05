#!/usr/bin/env python3
"""plot_acceptance.py -- visualize Metropolis acceptance diagnostics.

Two panels from a single run's diagnostic files:
  (1) acceptance vs sweep    -- from <base>.trace  (does it drift as the chain
      thermalizes?)
  (2) acceptance vs |q|      -- from <base>.accept (is the momentum-dependent
      step tuned to ~uniform acceptance across wave vectors, a la Troster 2013
      OFMC, or does it fall off at low q?)

Usage:
    uv run tools/plot_acceptance.py <base.dat> [out.png]
where <base.dat>.trace and <base.dat>.accept must exist (the engine writes them
when a run completes). Default output is <base.dat>.acceptance.png.
"""
import sys
import re
import numpy as np


def read_p8(path):
    try:
        with open(path) as f:
            head = f.read(500)
        m = re.search(r"p8=([\d.]+)", head)
        return float(m.group(1)) if m else None
    except OSError:
        return None


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    base = sys.argv[1]
    tracef, accf = base + ".trace", base + ".accept"
    out = sys.argv[2] if len(sys.argv) > 2 else base + ".acceptance.png"

    tr = np.loadtxt(tracef, comments="#")          # sweeps Delta2 rel_err accept wall_s
    ac = np.loadtxt(accf, comments="#")            # q1 q2 qmag proposed accepted rate
    p8 = read_p8(accf)

    sw, acc_sweep = tr[:, 0], tr[:, 3]
    q, prop, rate = ac[:, 2], ac[:, 3], ac[:, 5]

    # bin acceptance-vs-|q| into log-spaced |q| bins: proposal-weighted mean +/- spread
    edges = np.geomspace(q.min(), q.max(), 26)
    ctr, mean, lo, hi = [], [], [], []
    for i in range(len(edges) - 1):
        m = (q >= edges[i]) & (q < edges[i + 1])
        if m.sum() < 2:
            continue
        w = prop[m]
        mu = np.average(rate[m], weights=w)
        ctr.append(np.sqrt(edges[i] * edges[i + 1]))
        mean.append(mu)
        lo.append(np.percentile(rate[m], 10))
        hi.append(np.percentile(rate[m], 90))
    ctr, mean, lo, hi = map(np.array, (ctr, mean, lo, hi))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

    # (1) acceptance vs sweep
    axL.plot(sw, acc_sweep, "-o", ms=3, color="tab:blue")
    axL.axhspan(0.35, 0.55, color="tab:green", alpha=0.12,
                label="ideal 35-55%")
    axL.axhline(0.5, color="0.5", ls="--", lw=0.8)
    axL.set_xlabel("sweep")
    axL.set_ylabel("Metropolis acceptance (per block)")
    axL.set_title("acceptance vs sweep")
    axL.set_ylim(0, 1)
    axL.legend(loc="best", fontsize=9)
    axL.grid(alpha=0.3)

    # (2) acceptance vs |q|
    axR.scatter(q, rate, s=4, color="0.7", alpha=0.35, label="per mode")
    axR.plot(ctr, mean, "-o", ms=4, color="tab:red",
             label="proposal-weighted mean")
    axR.fill_between(ctr, lo, hi, color="tab:red", alpha=0.15,
                     label="10-90 pct")
    axR.axhspan(0.35, 0.55, color="tab:green", alpha=0.12)
    axR.axhline(0.5, color="0.5", ls="--", lw=0.8)
    if p8:
        axR.axvline(p8, color="tab:purple", ls=":", lw=1.3,
                    label=f"crossover q$_c\\approx$p8={p8:g}")
    axR.set_xscale("log")
    axR.set_xlabel("|q|")
    axR.set_ylabel("Metropolis acceptance (per mode)")
    axR.set_title("acceptance vs |q|")
    axR.set_ylim(0, 1)
    axR.legend(loc="best", fontsize=9)
    axR.grid(alpha=0.3, which="both")

    fig.suptitle(f"Metropolis acceptance diagnostics: {base}")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
