#!/usr/bin/env python3
"""reformat_legacy.py -- convert a LEGACY brane dump (example_data/N=<N>.dat)
into the modern .dat format (rich key=value header + q1 q2 qx qy qmag G Gerr
Ginv columns), so the same analysis tools work on both.

Legacy layout (legacy/storage.c dump): L*L modes in row-major (q1=0..L-1,
q2=0..L-1) order, three lines each -- "c0 c1", "Re Im", "g" -- then a trailing
"C px0 px1". G(q) = g / c1. The legacy engine is a SINGLE chain with no replicas,
so there is NO per-mode error estimate: Gerr is written as 0 and flagged NA in
the header. Coupling p8 and N are not stored; p8 defaults to 0.3 (legacy default)
and N comes from the filename. Unknown run params are written as NA.

Usage:
    uv run tools/reformat_legacy.py example_data/N=100.dat            # in place
    uv run tools/reformat_legacy.py example_data/N=100.dat --p8 0.3 --out /tmp/x.dat
"""
import argparse
import math
import re
import sys
import numpy as np

PI = math.pi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--p8", type=float, default=0.3, help="legacy coupling (default 0.3)")
    ap.add_argument("--out", default=None, help="output path (default: overwrite infile)")
    ap.add_argument("--d0", type=float, default=2.6, help="legacy base step (default 2.6)")
    a = ap.parse_args()

    m = re.search(r"N=(\d+)", a.infile)
    if not m:
        sys.exit(f"cannot find N=<int> in filename {a.infile}")
    N = int(m.group(1)); L = 2 * N + 1; sp = 2.0 * PI / L
    p8 = a.p8; Y = (2.0 * PI / 3.0) * p8 * p8; N8 = int(p8 / PI * N)

    toks = np.fromstring(open(a.infile).read().replace("\t", " "), sep=" ")
    need = L * L * 5
    if toks.size < need:
        sys.exit(f"{a.infile}: expected >= {need} numbers for N={N}, got {toks.size} "
                 "(is this already reformatted?)")
    body = toks[:need].reshape(L * L, 5)          # [c0, c1, re, im, g]
    c0, c1, g = body[:, 0], body[:, 1], body[:, 4]
    trailer = toks[need:need + 3]
    C = trailer[0] if trailer.size >= 1 else float(np.max(c1))
    px0 = trailer[1] if trailer.size >= 2 else float("nan")
    px1 = trailer[2] if trailer.size >= 3 else float("nan")

    # G(q) = g / c1 over the row-major (unsigned) storage grid.
    with np.errstate(divide="ignore", invalid="ignore"):
        Gstore = np.where(c1 > 0, g / c1, 0.0)     # length L*L, index = qi1*L+qi2

    # SN[unsigned] = sin(sp * signed_q)
    uidx = np.arange(L)
    sgn = np.where(uidx <= N, uidx, uidx - L)
    SN = np.sin(sp * sgn)

    # Poisson ratio via the legacy calcPR formula (window q in [-N8, N8)).
    nu = float("nan")
    if N8 > 0 and C > 0 and np.isfinite(px0) and np.isfinite(px1):
        Kx = Ky = 0.0
        for q1 in range(-N8, N8):
            for q2 in range(-N8, N8):
                gg = Gstore[((q1 + L) % L) * L + ((q2 + L) % L)]
                Kx += SN[(q1 + L) % L] ** 2 * gg
                Ky += SN[(q2 + L) % L] ** 2 * gg
        denom = px0 / C - Kx * Kx
        if denom != 0:
            nu = -(px1 / C - Kx * Ky) / denom

    out = a.out or a.infile
    samples = int(C) if np.isfinite(C) else "NA"
    with open(out + ".tmp", "w") as f:
        f.write("# Fourier MC membrane (legacy example_data, reformatted)\n")
        f.write(f"# N={N} L={L} n={N} p8={p8:.4f} N8={N8} Y={Y:.6f} d0={a.d0:.4f} seed=NA\n")
        f.write("# nt=1 it=NA cores=NA\n")   # legacy = single chain (intra-chain threads unknown)
        f.write(f"# therm=NA sweeps={samples} sweeps_cap=NA min_sweeps=NA block=NA "
                "meas_every=1 steps_per_sweep=NA\n")
        f.write("# eps=NA rel_err=NA converged=NA\n")
        nu_s = f"{nu:.6f}" if np.isfinite(nu) else "NA"
        f.write(f"# samples={samples} accept_rate=NA wall_s=NA nu={nu_s} nu_err=NA\n")
        f.write("# engine_sha=legacy source=example_data Gerr=NA(single-chain,no-error-estimate)\n")
        f.write("# q1 q2 qx qy qmag G Gerr Ginv\n")
        for q1 in range(-N, N + 1):
            i1 = (q1 + L) % L
            for q2 in range(-N, N + 1):
                if q1 == 0 and q2 == 0:
                    continue
                i2 = (q2 + L) % L
                G = Gstore[i1 * L + i2]
                qx = sp * q1; qy = sp * q2
                qm = math.sqrt(qx * qx + qy * qy)
                ginv = 1.0 / G if G > 0 else 0.0
                f.write(f"{q1}\t{q2}\t{qx:.8f}\t{qy:.8f}\t{qm:.8f}\t"
                        f"{G:.10e}\t{0.0:.10e}\t{ginv:.10e}\n")
    import os
    os.replace(out + ".tmp", out)
    print(f"[reformat] {a.infile} -> {out}  (N={N} L={L} p8={p8} samples={samples} nu={nu_s})")


if __name__ == "__main__":
    main()
