#!/usr/bin/env python3
"""Rigid vs non-rigid models of the ITS IB chips from the normal equations accumulated by milleSplit.cxx.

  milleSplitModels.py labels mp2tree.txt iblabels.txt          # chip TX label -> sensor ID, input of milleSplit
  milleSplitModels.py fit sens.txt ms_B0.bin [ms_more.bin ...] # models, from one or more summed accumulations

sens.txt is the per-sensor table of sensDump.C (columns: id layer hb stave chipInStave ...).
Models: rigid (6 DOFs per chip, shared by its y and z rows), rigid + a y-only TY term beta (derivative
tan(phi_inc) on the y rows only) grouped per sensor / layer x chip index / layer / layer x (end, middle),
and the full y/z split (12 per chip). The chi2 gains are exact (pseudo-inverse); the quoted parameter values
use a 50 um / 5 mrad prior, needed for the weak modes of a fit in which every chip is free.
The mean-vertex calibration (X,Y,Z,SlopeX,SlopeY of 8 time slots, accumulated after the chip parameters) is free in
every model: the records measure it relative to a reference at (0,0), ~1.8 mm from the beam, and fixing it fakes a
large rphi-only (y-only) term."""
import re, sys
import numpy as np

NS, NPS = 432, 12
NPIB = NS * NPS
NMV = 8 * 5  # mean-vertex calibration, always free
NP = NPIB + NMV


def writeLabels(treef, outf):
    nst = {0: 6, 1: 8, 2: 10}  # IB staves per half barrel
    off = {0: 0, 1: 108, 2: 252}  # first sensor ID of the layer
    n = 0
    with open(outf, "w") as out:
        for l in open(treef):
            m = re.search(r'ITS/ITSULayer([0-2])/ITSUHalfBarrel(\d)/ITSUStave(\d+)/ITSUChip(\d+) \(sens\) RB\[ TX\((\d+)\)', l)
            if m:
                L, h, s, c, lab = map(int, m.groups())
                out.write(f"{lab} {off[L] + (h * nst[L] + s) * 9 + c}\n")
                n += 1
    print(f"{n} IB chip labels written to {outf}")


def read(fn):
    with open(fn, 'rb') as f:
        n = np.frombuffer(f.read(4), np.int32)[0]
        if n != NP:
            sys.exit(f"{fn}: {n} parameters, {NP} expected; re-accumulate with the current milleSplit.cxx (mean vertex free)")
        N = np.frombuffer(f.read(8 * n * n), np.float64).reshape(n, n).copy()
        B = np.frombuffer(f.read(8 * n), np.float64).copy()
    return N, B


def fit(sensf, files):
    N, B = read(files[0])
    for fn in files[1:]:
        n, b = read(fn)
        N += n
        B += b
    sd = np.loadtxt(sensf)
    chip = sd[:, 4].astype(int)
    lay = sd[:, 1].astype(int)
    reg = 1 / 0.005**2

    def solve(Pc):
        P = np.zeros((NP, Pc.shape[1] + NMV))
        P[:NPIB, :Pc.shape[1]] = Pc
        P[NPIB:, Pc.shape[1]:] = np.eye(NMV)
        Nu = P.T @ N @ P
        Br = P.T @ B
        R = np.diag(np.r_[np.full(Pc.shape[1], reg), np.ones(NMV)])  # no real prior on the mean vertex (~1.8 mm off)
        x = np.linalg.solve(Nu + R, Br)
        err = np.sqrt(np.diag(np.linalg.inv(Nu + R)))
        w, V = np.linalg.eigh(Nu)
        k = w > w.max() * 1e-12
        gain = ((V[:, k].T @ Br)**2 / w[k]).sum()  # exact chi2 decrease
        return x[:Pc.shape[1]], err[:Pc.shape[1]], gain

    P0 = np.zeros((NPIB, NS * 6))  # rigid: y and z copies tied
    for s in range(NS):
        for k in range(6):
            P0[s * NPS + k, s * 6 + k] = P0[s * NPS + 6 + k, s * 6 + k] = 1
    _, _, g0 = solve(P0)
    print(f"rigid (6/chip): gain {g0:.0f}")

    def extra(groups, ng):  # y-only TY, i.e. the TY copy of the y rows
        E = np.zeros((NPIB, ng))
        for s in range(NS):
            E[s * NPS + 1, groups[s]] = 1
        return E

    models = {
        "beta per sensor": (np.arange(NS), NS),
        "beta per layer x chipidx": (lay * 9 + chip, 27),
        "beta per layer": (lay, 3),
        "beta per layer x (end/mid)": (lay * 2 + ((chip == 0) | (chip == 8)), 6),
    }
    res = {}
    for nm, (g, ng) in models.items():
        x, e, gn = solve(np.hstack([P0, extra(g, ng)]))
        res[nm] = (x[NS * 6:], e[NS * 6:])
        print(f"rigid + {nm:28s}: extra gain {gn - g0:8.0f} for {ng} params")
    _, _, gs = solve(np.eye(NPIB))
    print(f"full split (12/chip): extra gain {gs - g0:.0f} for {NPIB - NS * 6} params")
    b, e = res["beta per layer x chipidx"]
    for L in range(3):
        print(f" Lr{L} beta (y-only TY, um) per chip:  " + "  ".join(f"c{c}:{1e4 * b[L * 9 + c]:+5.1f}+-{1e4 * e[L * 9 + c]:.1f}" for c in range(9)))
    b, e = res["beta per sensor"]
    for L in range(3):
        for c in [0, 4, 8]:
            ss = np.where((lay == L) & (chip == c))[0]
            print(f"   Lr{L} c{c} per-stave beta: mean {1e4 * b[ss].mean():+5.1f}, rms {1e4 * b[ss].std():4.1f}, <err> {1e4 * e[ss].mean():4.1f}")


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "labels":
        writeLabels(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 4 and sys.argv[1] == "fit":
        fit(sys.argv[2], sys.argv[3:])
    else:
        sys.exit(__doc__)
