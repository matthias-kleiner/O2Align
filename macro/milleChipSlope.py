#!/usr/bin/env python3
"""Fits, from the Mille records, the IB chip RY (and optional incidence-angle terms) separately for the dy rows.

  milleChipSlope.py mp2data_0.bin mp2tree.txt NRECORDS

The local track parameters are profiled out per row; a term tx*z_loc (incidence tangent times the local z) that replaces
the per-chip RY is the signature of a mismatch between the plane x used for the track and the cluster x.
tz, tx are the incidence tangents in the sensor frame, taken from the TY derivatives of the z and y rows
(single-row approximation of the profiling)."""
import re, struct, sys, math, collections

binf, treef, NREC = sys.argv[1], sys.argv[2], int(sys.argv[3])
RE = re.compile(r'ITSULayer([0-2])/ITSUHalfBarrel(\d)/ITSUStave(\d+)/ITSUChip(\d+) \(sens\) RB\[(.*)\]')
lab = {}
for l in open(treef):
    m = RE.search(l)
    if m:
        for d, v in re.findall(r'(\w+)\((\d+)\)', m.group(5)):
            lab[int(v)] = (int(m.group(1)), int(m.group(4)), d)

def chol(A, n):
    Lm = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = A[i][j] - sum(Lm[i][k] * Lm[j][k] for k in range(j))
            if i == j:
                Lm[i][j] = math.sqrt(max(s, 1e-14))
            else:
                Lm[i][j] = s / Lm[j][j]
    return Lm

def cholsolve(Lm, rhs, n):
    y = [0.0] * n
    for i in range(n):
        y[i] = (rhs[i] - sum(Lm[i][k] * y[k] for k in range(i))) / Lm[i][i]
    x = [0.0] * n
    for i in reversed(range(n)):
        x[i] = (y[i] - sum(Lm[k][i] * x[k] for k in range(i + 1, n))) / Lm[i][i]
    return x

# per layer: parameter list: a_0..a_8, b_0..b_8, kz, kx  (index = 0..17, 18, 19)
NP = 21
N = {L: [[0.0] * NP for _ in range(NP)] for L in range(3)}
B = {L: [0.0] * NP for L in range(3)}
chi_sum = {L: 0.0 for L in range(3)}
cnt = {L: 0 for L in range(3)}

def parse(f):
    h = f.read(4)
    if len(h) < 4:
        return None
    n2 = struct.unpack('i', h)[0]
    dbl = n2 < 0
    n = abs(n2) // 2
    vals = struct.unpack(('%dd' if dbl else '%df') % n, f.read(n * (8 if dbl else 4)))
    ints = struct.unpack('%di' % n, f.read(4 * n))
    rows, state = [], 0
    for i in range(1, n):
        v, k = vals[i], ints[i]
        if k == 0:
            if state in (0, 2):
                rows.append([v, 0.0, {}, {}]); state = 1
            else:
                rows[-1][1] = v; state = 2
        elif state == 1:
            rows[-1][2][k] = v
        else:
            rows[-1][3][k] = v
    return [r for r in rows if r[1] > 0]

nrec = used = 0
with open(binf, 'rb') as f:
    while nrec < NREC:
        rows = parse(f)
        if rows is None:
            break
        nrec += 1
        locs = sorted({k for r in rows for k in r[2]})
        if not locs or len(locs) > 40:
            continue
        li = {k: i for i, k in enumerate(locs)}
        nl = len(locs)
        A = [[0.0] * nl for _ in range(nl)]
        bres = [0.0] * nl
        for r in rows:
            w = 1.0 / (r[1] * r[1])
            it = [(li[k], d) for k, d in r[2].items()]
            for i, di in it:
                bres[i] += w * di * r[0]
                for j, dj in it:
                    A[i][j] += w * di * dj
        Lm = chol(A, nl)
        xr = cholsolve(Lm, bres, nl)
        for ri, r in enumerate(rows):
            base = None
            for k in r[3]:
                t = lab.get(k)
                if t and t[2] == 'RY':
                    base = k
                    L, c = t[0], t[1]
            if base is None:
                continue
            txd, tzd = r[3].get(base - 4, 0.0), r[3].get(base - 2, 0.0)
            if abs(txd) <= abs(tzd) or ri + 1 >= len(rows):
                continue                              # not a y row
            z = rows[ri + 1]
            if (base - 4) not in z[3] and (base - 2) not in z[3]:
                continue
            if abs(z[3].get(base - 2, 0.0)) <= abs(z[3].get(base - 4, 0.0)):
                continue                              # next row is not the z row of the same chip
            tx = r[3].get(base - 3, 0.0)              # TY derivative of the y row: incidence in the bending plane
            tz = z[3].get(base - 3, 0.0)              # TY derivative of the z row: incidence along z
            dry = r[3][base]
            w = 1.0 / (r[1] * r[1])
            g = {c: txd, 9 + c: dry, 18: tz, 19: tx, 20: tx * dry}
            u = {p: [w * d * r[2].get(k, 0.0) for k in locs] for p, d in g.items()}
            x = {p: cholsolve(Lm, u[p], nl) for p in g}
            for p, dp in g.items():
                for q, dq in g.items():
                    N[L][p][q] += w * dp * dq - sum(u[p][i] * x[q][i] for i in range(nl))
                B[L][p] += w * dp * r[0] - sum(u[p][i] * xr[i] for i in range(nl))
            cnt[L] += 1
            used += 1
print('records', nrec, 'IB y rows', used, cnt)

def fit(L, active):
    idx = [i for i in range(NP) if i in active]
    n = len(idx)
    M = [[N[L][i][j] for j in idx] for i in idx]
    # regularise the empty classes
    for k in range(n):
        if M[k][k] <= 0:
            M[k][k] = 1.0
    v = [B[L][i] for i in idx]
    Lm = chol(M, n)
    x = cholsolve(Lm, v, n)
    err = []
    for k in range(n):
        e = [1.0 if j == k else 0.0 for j in range(n)]
        s = cholsolve(Lm, e, n)
        err.append(math.sqrt(max(s[k], 0.0)))
    gain = sum(x[k] * v[k] for k in range(n))      # chi2 decrease
    return dict(zip(idx, x)), dict(zip(idx, err)), gain

base_set = set(range(18))
for L in range(3):
    s0, e0, g0 = fit(L, base_set)
    s1, e1, g1 = fit(L, base_set | {18, 19})
    s2, e2, g2 = fit(L, {*range(9), 18, 19, 20})
    s3, e3, g3 = fit(L, base_set | {20})
    s4, e4, g4 = fit(L, {*range(9), 20})
    print('=== layer', L, ' chi2 gain: (TX,RY) %.0f | (TX,RY,kz,kx) %.0f | (TX,kz,kx) %.0f' % (g0, g1, g2))
    print('   chi2 gain (TX,RY,tx*zloc) %.0f | (TX,tx*zloc only) %.0f | (TX,kz,kx,tx*zloc) %.0f ;  coef tx*zloc: %.0f +- %.0f urad/unit tx (with RY), %.0f +- %.0f (no RY)' % (g3, g4, g2, s3[20]*1e6, e3[20]*1e6, s4[20]*1e6, e4[20]*1e6))
    print('   chip: RY (TX,RY)        RY (TX,RY,tx*zloc)')
    for c in range(9):
        print('   c%d  %7.0f +-%4.0f      %7.0f +-%4.0f' % (c, s0[9 + c] * 1e6, e0[9 + c] * 1e6, s3[9 + c] * 1e6, e3[9 + c] * 1e6))
