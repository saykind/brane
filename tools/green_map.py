#!/usr/bin/env python3
"""
green_map.py -- 2D heatmaps of the Green function G(qx, qy) and its inverse
G^-1(qx, qy) over the Brillouin zone, from a brane output file.

G(q) = <|h_q|^2> spans a huge dynamic range (~q^-4, diverging at q->0), so the
color scale is logarithmic. The zero mode (0,0) is left blank. These maps make
the lattice anisotropy visible: near q=0 the contours are circular (isotropic,
G ~ q^-4), and they distort toward the square Brillouin-zone symmetry at large q.

Usage:
    uv run tools/green_map.py data/N=40.dat
    uv run tools/green_map.py data/N=40.dat --png data/green_N40.png
"""
import argparse
import sys
import numpy as np


def load_grid(path):
    """Return qx1d, qy1d, G2d, Ginv2d, header (2d arrays indexed [iy, ix])."""
    header = {}
    saw = False
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "qmag" in line and "Ginv" in line:
                    saw = True
                for tok in line[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1); header[k] = v
                continue
            break
    if not saw:
        sys.exit(f"error: '{path}' is not a brane new-format file "
                 "(regenerate with ./brane ... out=...).")
    d = np.loadtxt(path, comments="#")
    q1 = d[:, 0].astype(int); q2 = d[:, 1].astype(int)
    qx = d[:, 2]; qy = d[:, 3]; G = d[:, 5]; Ginv = d[:, 6]
    N = int(header.get("N", (q1.max())))
    L = 2 * N + 1
    a = 2 * np.pi / L
    G2 = np.full((L, L), np.nan)
    Gi2 = np.full((L, L), np.nan)
    # index: ix = q1+N (qx axis), iy = q2+N (qy axis)
    G2[q2 + N, q1 + N] = np.where(G > 0, G, np.nan)
    Gi2[q2 + N, q1 + N] = np.where(Ginv > 0, Ginv, np.nan)
    qax = a * np.arange(-N, N + 1)
    return qax, G2, Gi2, header, a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datfile")
    ap.add_argument("--png", default=None)
    args = ap.parse_args()

    qax, G2, Gi2, header, a = load_grid(args.datfile)
    p8 = header.get("p8", "?"); N = header.get("N", "?")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
    except ImportError:
        sys.exit("matplotlib required for green_map (uv sync).")

    ext = [qax[0] - a / 2, qax[-1] + a / 2, qax[0] - a / 2, qax[-1] + a / 2]
    fig, (axG, axI) = plt.subplots(1, 2, figsize=(12, 5.2))

    def finite(a2):
        v = a2[np.isfinite(a2)]
        return v.min(), v.max()

    gmin, gmax = finite(G2)
    imG = axG.imshow(G2, origin="lower", extent=ext, cmap="magma",
                     norm=LogNorm(vmin=gmin, vmax=gmax), aspect="equal")
    fig.colorbar(imG, ax=axG, label=r"$G(q)=\langle|h_q|^2\rangle$")
    axG.set_title(rf"Green function $G(q_x,q_y)$  ($N={N},\ p_8={p8}$)")

    imin, imax = finite(Gi2)
    imI = axI.imshow(Gi2, origin="lower", extent=ext, cmap="viridis",
                     norm=LogNorm(vmin=imin, vmax=imax), aspect="equal")
    fig.colorbar(imI, ax=axI, label=r"$G^{-1}(q)$")
    axI.set_title(r"Inverse Green $G^{-1}(q_x,q_y)$")

    for ax in (axG, axI):
        ax.set_xlabel(r"$q_x$"); ax.set_ylabel(r"$q_y$")

    fig.tight_layout()
    png = args.png or (args.datfile.rsplit(".", 1)[0] + "_green.png")
    fig.savefig(png, dpi=140)
    print(f"[plot] wrote {png}")


if __name__ == "__main__":
    main()
