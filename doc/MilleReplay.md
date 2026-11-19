# Validating a pede solution on the Mille data (`o2-dev-mille-replay`)

`o2-dev-mille-replay` ([src/milleReplay.cxx](../src/milleReplay.cxx), built and installed with the package) replays the Mille binaries written by `--output MilleData`, with a pede solution applied, and reports how
the data take it, per input file. It links against nothing but the standard library (no O2, ROOT or pede at run time) and needs no new data pass: all the TRD/TOF collision and cosmic
binaries of a fit (about 9 M records) take about 3 minutes on 28 threads. Regenerating the `gbltrack` debug trees for the
same check takes hours.

## What it computes

For every record, `o2-dev-mille-replay`:
1. subtracts the solution from the residuals, r → r − Σ g·p, using the global derivatives stored in the record;
2. refits the local (track) parameters as pede does in its final iteration: pass 1 without down-weighting, passes 2–3
   Huber, then Cauchy, each from the residuals of the previous pass (`-dw`, default 4 passes as `outlierdownweighting 4`);
3. applies pede's rejections: the down-weighted χ² of the last pass above `-cut` × the 3σ quantile of χ²(ndf), and a
   down-weight fraction above `-dwcut` (default 0.15);
4. for the accepted records, accumulates per global parameter j (free or fixed) the gradient term b = Σ w·g·r and the
   curvature with the local parameters profiled out, Hp = Σ w·g² − uᵀA⁻¹u.

The resulting rejection fractions reproduce pede's per-file ones within half a percent.

b/Hp is the shift the data still ask for in that parameter alone, b²/Hp the χ² it would gain. Since they ignore the
correlations between parameters, they are a diagnostic, not a new solution:
- a free parameter of a converged fit has b ≈ 0, unless it enters a constraint: the common mode of children constrained
  to a free parent is balanced by the Lagrange multiplier, so judge such modes by their group (below);
- a fixed parameter shows what freeing it would bring;
- running each dataset separately (B0, B+, B−, cosmics) separates geometry, which every dataset accepts, from effects
  that depend on the field or the sample, where the datasets ask for opposite shifts.

Groups add the derivatives of all the parameters whose volume matches a regex and whose DOF has a given name, and treat
the sum as one parameter: the common radial shift of the chambers of a TRD layer, the common rφ shift of all sectors
(a global rotation), etc.

## Running

```bash
source install/init.sh   # puts install/bin on PATH
o2-dev-mille-replay -tree <recording>/mp2tree.txt -res <fit>/millepede.res -groups groups.txt \
    -sel 'TRD.*|TOF.*|PVT.*|ITS|ITS/ITSULayer[0-6]/ITSUHalfBarrel[0-1]' -j 28 -o rB0 B0/mp2data_*.bin B0Bis/mp2data_*.bin
```

| option | meaning |
|---|---|
| `-tree` | `mp2tree.txt` of the recording run: it lists the labels of every DOF present in the records (the tree of a later, more restricted fit omits the fixed ones) |
| `-res` | solution to apply; without it the records are replayed at the recorded reference |
| `-con` | add the `Measurement` priors of an `mp2con.txt` to b and Hp |
| `-groups` | lines `name DOF regex`, e.g. `TRD_L5_radial TX TRD/sm[0-9]+/st[0-4]/pl5`; `#` starts a comment |
| `-sel` | accumulate per-parameter results only for the volumes matching the regex (default: all); groups are always accumulated |
| `-dw`, `-dwcut`, `-cut` | pede's `outlierdownweighting`, `dwfractioncut` and the final χ² cut factor |
| `-maxloc`, `-maxrec` | skip records with more local parameters; process at most N records per file |
| `-j` | threads, one input file each at a time |

Records containing a `PVT/meanVertex` label (the multi-track vertex records, or tracks with a vertex constraint) are
counted separately from the others.

## Output

All the sums are additive over input files.

| file | content |
|---|---|
| `<out>_files.txt` | per input file: records, skipped; for records without and with vertex labels: accepted, rejected by χ², rejected by down-weight fraction, ⟨χ²/ndf⟩ of the accepted |
| `<out>_pars.txt` | per label: volume, DOF, time slot, solution value, records, rows, b, Hu (w/o profiling), Hp, b/Hp, b²/Hp (cm, rad) |
| `<out>_groups.txt` | the same per group |

Each parameter is in the frame of its volume's DOF set. For TRD chambers that is the sensor frame (`TX` radial, `TY` rφ,
`TZ` along the beam), while for TRD supermodules it is the TGeo frame, see [DOFConfig_rules.md](DOFConfig_rules.md#trd-alignable-match-strings).
