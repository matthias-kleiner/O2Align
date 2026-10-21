#!/usr/bin/env python3
"""Per-layer summary of the fitted IB chip RY (rotation about the sensor normal) in result.json.

  chipRYSummary.py result.json [injectChipRY.txt]

Without the injection list: mean/rms of the chip RY per IB layer, in urad (0 is expected for a control
run started from the already fitted geometry).
With it (closure): regression of the fitted RY on the injected one. A consistent chain gives slope -1
(the fit undoes the injection), corr ~ 1; slope +1 means a sign error between the derivatives and the
geometry; slope ~ 0 means the RY is not seen by the fit.
"""
import json
import re
import statistics as st
import sys

RE = re.compile(r'ITS/ITSULayer([0-2])/ITSUHalfBarrel\d+/ITSUStave\d+/ITSUChip\d+$')

res = {e['symName']: e['rigidBody'][4] for e in json.load(open(sys.argv[1])) if RE.match(e['symName'])}
inj = {}
if len(sys.argv) > 2:
    for line in open(sys.argv[2]):
        n, v = line.split()
        inj[n] = float(v)

for L in '012':
    names = [n for n in res if RE.match(n).group(1) == L]
    fit = [res[n] * 1e6 for n in names]
    if not fit:
        print(f'Lr{L}: no chip in the result')
        continue
    s = f'Lr{L}: {len(fit)} chips, fitted RY mean {st.mean(fit):8.1f} rms {st.pstdev(fit):7.1f} urad'
    ns = [n for n in names if n in inj]
    if ns:
        x = [inj[n] * 1e6 for n in ns]
        y = [res[n] * 1e6 for n in ns]
        mx, my = st.mean(x), st.mean(y)
        sxx = sum((a - mx) ** 2 for a in x)
        sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
        syy = sum((b - my) ** 2 for b in y)
        slope = sxy / sxx if sxx > 0 else float('nan')
        corr = sxy / (sxx * syy) ** 0.5 if sxx > 0 and syy > 0 else float('nan')
        s += f' | injected mean {mx:8.1f} | fit = {slope:+.2f}*inj {my - slope * mx:+.1f} (corr {corr:+.2f})'
    print(s)
