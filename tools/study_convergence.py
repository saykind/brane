#!/usr/bin/env python3
"""study_convergence.py -- parse a verbose `brane` run log and plot how the
statistical error and the observable evolve with the number of sweeps.

The engine's verbose lines look like:
    sweeps=20  Delta2=13753.059473  Delta2 rel.err=0.0589  (target 0.0000)

We extract (sweeps, Delta2, rel_err) and produce a two-panel figure:
  left  : Delta2 vs sweeps            -> thermalization (running mean plateaus)
  right : rel_err vs sweeps (log-log) -> error decay, with a 1/sqrt(n) reference

Usage:
    uv run tools/study_convergence.py study/convergence_N100_p0.40.log \
        --png plots/study/convergence_N100_p0.40.png
"""
import argparse
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LINE = re.compile(r"sweeps=(\d+)\s+Delta2=([-\d.eE+]+)\s+Delta2 rel\.err=([-\d.eE+]+)")


def parse(path):
    """Read (sweeps, Delta2, rel_err) from either a <out>.dat.trace file
    (tab columns: sweeps Delta2 rel_err wall_s) or a verbose stdout log."""
    sw, d2, re_ = [], [], []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        m = LINE.search(line)
        if m:                                   # verbose stdout log format
            sw.append(int(m.group(1)))
            d2.append(float(m.group(2)))
            re_.append(float(m.group(3)))
            continue
        parts = line.split()                    # trace-file columns
        if len(parts) >= 3:
            try:
                sw.append(int(float(parts[0])))
                d2.append(float(parts[1]))
                re_.append(float(parts[2]))
            except ValueError:
                pass
    return np.array(sw), np.array(d2), np.array(re_)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--png", default="plots/study/convergence.png")
    ap.add_argument("--therm-tol", type=float, default=0.01,
                    help="Delta2 within this rel. band of its final value = thermalized")
    a = ap.parse_args()

    sw, d2, rel = parse(a.log)
    if len(sw) < 3:
        sys.exit(f"not enough data points parsed from {a.log}")

    # --- thermalization: first sweep where running-mean Delta2 is within tol of final
    final = d2[-1]
    within = np.abs(d2 - final) / abs(final) < a.therm_tol
    therm_sw = int(sw[np.argmax(within)]) if within.any() else None

    # --- error decay: fit rel_err = C * sweeps^p on the valid (>0) tail
    good = rel > 0
    p = C = None
    if good.sum() >= 3:
        lp = np.polyfit(np.log(sw[good]), np.log(rel[good]), 1)
        p, C = lp[0], np.exp(lp[1])

    import os
    os.makedirs(os.path.dirname(a.png) or ".", exist_ok=True)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5))

    axL.plot(sw, d2, "o-", ms=4)
    if therm_sw is not None:
        axL.axvline(therm_sw, color="tab:green", ls=":",
                    label=f"within {a.therm_tol:.0%} of final @ {therm_sw} sweeps")
        axL.legend()
    axL.set_xlabel("measurement sweeps"); axL.set_ylabel(r"$\Delta_2=\sum_q\langle|h_q|^2\rangle$ (running mean)")
    axL.set_title("Thermalization")

    axR.loglog(sw[good], rel[good], "o-", ms=4, label="measured rel.err")
    if p is not None:
        axR.loglog(sw[good], C * sw[good] ** p, "--",
                   label=f"fit ∝ n^{p:.2f}")
        axR.loglog(sw[good], rel[good][0] * (sw[good] / sw[good][0]) ** -0.5, ":",
                   color="gray", label=r"ideal $n^{-1/2}$")
    axR.set_xlabel("measurement sweeps"); axR.set_ylabel(r"rel. error on $\Delta_2$")
    axR.set_title("Statistical error decay"); axR.legend()

    fig.tight_layout(); fig.savefig(a.png, dpi=140); plt.close(fig)
    print(f"[plot] wrote {a.png}")

    # --- report ---
    print(f"points            : {len(sw)}  (sweeps {sw[0]}..{sw[-1]})")
    print(f"Delta2 final      : {final:.4f}")
    if therm_sw is not None:
        print(f"thermalized by    : ~{therm_sw} sweeps (running mean within {a.therm_tol:.0%})")
    if p is not None:
        print(f"error decay slope : {p:.3f}  (ideal -0.5 = 1/sqrt(n))")
        for eps in (0.01, 0.005, 0.003, 0.002):
            n = (eps / C) ** (1.0 / p)
            print(f"  rel.err < {eps:<6}: ~{n:.0f} sweeps")


if __name__ == "__main__":
    main()
