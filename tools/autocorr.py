#!/usr/bin/env python3
"""autocorr.py -- integrated autocorrelation time tau from a per-sweep Delta2
series (written by `brane ... series=<path>`).

tau_int tells you how many sweeps the chain needs to produce one *independent*
sample: the number of effective samples is N_eff = N / (2 tau_int), and the
statistical error scales as 1/sqrt(N_eff). Reducing tau ("decorrelation") is the
lever that gets more information per sweep -- over-relaxation, parallel
tempering, smarter moves all aim to lower it.

Uses the standard normalized ACF with Sokal automatic windowing (window W is the
smallest with W >= c*tau(W), c=5).

Usage:
    uv run tools/autocorr.py <series-file> [--c 5] [--burn 0]
"""
import argparse
import numpy as np


def integrated_tau(x, c=5.0):
    x = np.asarray(x, float)
    n = len(x)
    x = x - x.mean()
    var = np.dot(x, x) / n
    if var == 0:
        return 0.5, 0, np.array([1.0])
    # normalized autocorrelation via FFT
    f = np.fft.rfft(x, n=2 * n)
    acf = np.fft.irfft(f * np.conjugate(f))[:n].real
    acf /= acf[0]
    tau = 0.5
    W = n - 1
    for t in range(1, n):
        tau += acf[t]
        if t >= c * (2 * tau) and tau > 0:   # Sokal window: W >= c*tau_int
            W = t
            break
    return tau, W, acf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("series")
    ap.add_argument("--c", type=float, default=5.0, help="Sokal window factor")
    ap.add_argument("--burn", type=int, default=0, help="drop first N sweeps")
    a = ap.parse_args()

    data = np.loadtxt(a.series, comments="#")
    x = data[:, 1] if data.ndim == 2 else data
    if a.burn:
        x = x[a.burn:]
    n = len(x)
    tau, W, acf = integrated_tau(x, a.c)
    neff = n / (2 * tau)
    print(f"series          : {a.series}")
    print(f"samples (sweeps): {n}  (burn={a.burn})")
    print(f"tau_int         : {tau:.2f} sweeps   (window W={W})")
    print(f"N_eff           : {neff:.0f} independent samples  (= N / 2 tau)")
    print(f"decorrelation   : ~{2*tau:.1f} sweeps per independent sample")
    # a couple of ACF values for context
    for t in (1, 2, 5, 10, 20, 50):
        if t < len(acf):
            print(f"  ACF({t:3d}) = {acf[t]:+.3f}")


if __name__ == "__main__":
    main()
