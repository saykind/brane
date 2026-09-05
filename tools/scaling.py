#!/usr/bin/env python3
"""
scaling.py -- Measure how the new replica-parallel engine scales with the
number of cores (threads/replicas), and emit a table + plot.

Throughput is single-mode updates/sec = nt * sweeps * (l*l) / walltime.
Because sweeps and l*l are fixed, speedup(nt) = nt * t(1) / t(nt).

Usage:
    uv run tools/scaling.py                 # N=40, default thread list
    uv run tools/scaling.py --N 36 --sweeps 20 --reps 2
"""
import argparse
import re
import subprocess
import sys


def run_one(binary, N, p8, nt, sweeps):
    out = subprocess.run(
        [binary, f"N={N}", f"p8={p8}", f"nt={nt}",
         "therm=0", f"sweeps={sweeps}", "out=data/scaling.dat"],
        capture_output=True, text=True, check=True).stdout
    m = re.search(r"time = ([\d.]+) s", out)
    if not m:
        sys.exit(f"could not parse time from:\n{out}")
    return float(m.group(1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=40)
    ap.add_argument("--p8", type=float, default=0.4)
    ap.add_argument("--sweeps", type=int, default=12)
    ap.add_argument("--reps", type=int, default=2, help="repeats; take best time")
    ap.add_argument("--binary", default="./brane")
    ap.add_argument("--threads", default="1,2,4,6,8,10,12,14,16")
    ap.add_argument("--png", default="plots/scaling.png")
    args = ap.parse_args()

    threads = [int(x) for x in args.threads.split(",")]
    L = 2 * args.N + 1
    upd_per_sweep = L * L

    print(f"scaling: N={args.N} (L={L}) p8={args.p8} sweeps={args.sweeps} "
          f"reps={args.reps}\n")
    rows = []
    for nt in threads:
        t = min(run_one(args.binary, args.N, args.p8, nt, args.sweeps)
                for _ in range(args.reps))
        thru = nt * args.sweeps * upd_per_sweep / t
        rows.append((nt, t, thru))

    t1 = rows[0][1]
    hdr = f"{'cores':>5} {'time[s]':>9} {'Mupd/s':>9} {'speedup':>9} {'effic.':>8}"
    print(hdr)
    print("-" * len(hdr))
    table = []
    for nt, t, thru in rows:
        speedup = nt * t1 / t
        effic = speedup / nt
        table.append((nt, t, thru / 1e6, speedup, effic))
        print(f"{nt:5d} {t:9.2f} {thru/1e6:9.3f} {speedup:8.2f}x {100*effic:7.0f}%")

    # ---- plot -----------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not available; table only. `uv run tools/scaling.py`)")
        return

    nts = [r[0] for r in table]
    speed = [r[3] for r in table]
    eff = [r[4] for r in table]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.plot(nts, nts, "--", color="0.6", label="ideal (linear)")
    ax1.plot(nts, speed, "o-", color="C0", lw=2, label="measured")
    ax1.axvspan(12, max(nts), color="C1", alpha=0.10,
                label="E-cores (12P+4E)")
    ax1.set_xlabel("threads / replicas")
    ax1.set_ylabel("speedup vs 1 core")
    ax1.set_title(f"Replica-parallel speedup (N={args.N})")
    ax1.legend(frameon=False)
    ax1.grid(alpha=0.3)

    ax2.plot(nts, [100 * e for e in eff], "s-", color="C2", lw=2)
    ax2.axhline(100, ls="--", color="0.6")
    ax2.axvspan(12, max(nts), color="C1", alpha=0.10)
    ax2.set_xlabel("threads / replicas")
    ax2.set_ylabel("parallel efficiency  [%]")
    ax2.set_ylim(0, 110)
    ax2.set_title("4-12: memory bandwidth; 12-16: E-cores")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    import os
    os.makedirs(os.path.dirname(args.png) or ".", exist_ok=True)
    fig.savefig(args.png, dpi=140)
    print(f"\n[plot] wrote {args.png}")


if __name__ == "__main__":
    main()
