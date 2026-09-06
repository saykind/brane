#!/usr/bin/env python3
"""
braneio.py -- Shared readers for brane .dat output (header + columns).

Consolidates the parsing that analyze / reformat_legacy / tau_q /
plot_acceptance each used to re-implement. Two formats exist:

* Modern brane files carry a leading key=value '#' header followed by columns
      q1 q2 qx qy qmag G Gerr Ginv
  (older files predate Gerr and have 7 columns; Ginv is then the last column).
* Legacy example_data dumps use a different multi-line layout, handled by the
  read_legacy_raw() helper (used by tools/reformat_legacy.py to convert the raw
  thesis dumps into the modern format).

Only numpy is required.
"""
import os
import re
import sys
import numpy as np


# ---- modern format ---------------------------------------------------------
def read_header(path):
    """Parse the leading '#' key=value header of a brane .dat/.trace/.accept
    file. Returns (header_dict, saw_col_header) where saw_col_header is True iff
    the 'q1 q2 ... qmag ... Ginv' column line is present (i.e. a modern brane
    output file, not the legacy multi-line dump)."""
    header, saw_col_header = {}, False
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if "qmag" in line and "Ginv" in line:
                saw_col_header = True
            for tok in line[1:].split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    header[k] = v
    return header, saw_col_header


def _require_modern(path, saw_col_header):
    if not saw_col_header:
        sys.exit(
            f"error: '{path}' is not a brane (new-format) output file.\n"
            "  It looks like the legacy multi-line .dat format. Regenerate with\n"
            f"      ./brane N=40 p8=0.4 out={path}\n"
            "  or convert it with tools/reformat_legacy.py.\n"
        )


def load(path):
    """Return (qmag, G, Gerr, Ginv, header) for a modern brane .dat file."""
    header, saw = read_header(path)
    _require_modern(path, saw)
    data = np.loadtxt(path, comments="#")
    # columns: q1 q2 qx qy qmag G Gerr Ginv   (Gerr added; older files lack it)
    if data.shape[1] >= 8:
        return data[:, 4], data[:, 5], data[:, 6], data[:, 7], header  # qmag,G,Gerr,Ginv
    # backward compat: old 7-column files (no Gerr)
    return data[:, 4], data[:, 5], np.zeros(len(data)), data[:, 6], header


def p8_from_header(path, default=None, nbytes=500):
    """Read p8 from the header of a brane file (.dat/.trace/.accept/qseries).
    Returns a float, or `default` if not found / unreadable."""
    try:
        with open(path) as f:
            head = f.read(nbytes)
    except OSError:
        return default
    m = re.search(r"p8=([\d.]+)", head)
    return float(m.group(1)) if m else default


# ---- legacy example_data format --------------------------------------------
def read_legacy_raw(path):
    """Tokenize a LEGACY brane dump (example_data/N=<N>.dat) into its raw parts.

    Legacy layout (legacy/storage.c dump): L*L modes in row-major
    (q1=0..L-1, q2=0..L-1) order, three lines each -- "c0 c1", "Re Im", "g" --
    followed by a trailing "C px0 px1". Here c1 is the measurement count and
    g = sum_measurements |h_q|^2, so G(q) = g / c1.

    N is taken from the filename. Returns (body, trailer, N, L, a) where
    body is [L*L, 5] = [c0, c1, re, im, g], trailer is up to 3 numbers
    [C, px0, px1], and a = 2*pi/L (continuum convention).
    """
    m = re.search(r"N=(\d+)", path)
    if not m:
        sys.exit(f"cannot find N=<int> in filename {path}")
    N = int(m.group(1)); L = 2 * N + 1; a = 2 * np.pi / L
    toks = np.fromstring(open(path).read().replace("\t", " "), sep=" ")
    need = L * L * 5
    if toks.size < need:
        sys.exit(f"{path}: expected >= {need} numbers for N={N}, got {toks.size} "
                 "(is this already reformatted?)")
    body = toks[:need].reshape(L * L, 5)      # [c0, c1, re, im, g]
    trailer = toks[need:need + 3]
    return body, trailer, N, L, a


# ---- rotational (radial) averaging -----------------------------------------
def radial_average(qmag, G, nbins, Gerr=None):
    """Rotationally average G over log-spaced |q| shells.

    Returns (qr, Gr, Ginv_r, cnt, Ginv_err). Ginv_err is the propagated
    statistical error of 1/<G> in each shell: if per-mode errors Gerr (from the
    replica spread) are given, the shell error of the mean is
    sqrt(sum Gerr_i^2)/n, else it falls back to the in-shell std / sqrt(n).
    """
    mask = (qmag > 0) & (G > 0) & np.isfinite(G)
    q, g = qmag[mask], G[mask]
    ge = (Gerr[mask] if Gerr is not None else np.zeros_like(g))
    edges = np.logspace(np.log10(q.min()), np.log10(q.max()), nbins + 1)
    idx = np.digitize(q, edges)
    qr, gr, cnt, gerr = [], [], [], []
    for b in range(1, nbins + 1):
        sel = idx == b
        n = int(sel.sum())
        if n < 1:
            continue
        gm = g[sel].mean()
        # statistical error of the shell mean
        if Gerr is not None and np.any(ge[sel] > 0):
            em = np.sqrt(np.sum(ge[sel] ** 2)) / n
        else:
            em = (g[sel].std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        qr.append(q[sel].mean()); gr.append(gm); cnt.append(n); gerr.append(em)
    qr, gr, cnt, gerr = map(np.array, (qr, gr, cnt, gerr))
    Ginv_r = 1.0 / gr
    Ginv_err = gerr / gr ** 2          # error of 1/<G>
    return qr, gr, Ginv_r, cnt, Ginv_err
