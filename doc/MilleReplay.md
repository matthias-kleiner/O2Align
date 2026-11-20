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
   curvature with the local parameters profiled out, Hp = Σ w·g² − uᵀA⁻¹u (defined below).

The resulting rejection fractions reproduce pede's per-file ones within half a percent.

### Definitions

Everything is per record (a track, or a composed multi-track trajectory). The sums run over the measurement rows i of
the record and are then added over the records. Row i of a Mille record holds:
- r_i, the residual (measurement − prediction);
- σ_i, its error;
- the local derivatives d_ik = ∂pred_i/∂q_k with respect to the record's own track parameters q_k (GBL offsets and
  kinks);
- the global derivatives g_ij = ∂pred_i/∂p_j with respect to the alignment parameters p_j.

With the solution applied, r_i → r_i − Σ_j g_ij p_j.

| symbol | definition | meaning |
|---|---|---|
| w_i | ω_i / σ_i² | weight of the row in the final local fit; ω_i ∈ [0,1] is pede's outlier down-weight of the last pass (Huber/Cauchy on the normalized residual of the previous pass), 1 for a well-behaved point |
| A_kl | Σ_i w_i d_ik d_il | n_loc × n_loc normal matrix of the local (track) fit |
| u_k | Σ_i w_i g_ij d_ik | mixed global–local block: how much the derivative of p_j looks like a change of the track parameters |
| Hu | Σ_i w_i g_ij² | curvature of χ² in p_j with the track parameters fixed |
| Hp | Hu − uᵀ A⁻¹ u | curvature with the track parameters refitted (Schur complement), i.e. what the tracks cannot absorb |
| r_i^fit | r_i − Σ_k d_ik q̂_k | residual after the final local fit |
| b | Σ_i w_i g_ij r_i^fit | gradient of −χ²/2 in p_j |

The local χ² of a record, as a function of a correction Δq of the track parameters and a shift Δp_j of one global
parameter, is χ² = Σ_i w_i (r_i^fit − Σ_k d_ik Δq_k − g_ij Δp_j)². Eliminating Δq, i.e. letting the track refit
itself for every Δp_j, gives χ²(Δp_j) = χ²₀ − 2 b Δp_j + Hp Δp_j². b needs no profiling term because the local fit
leaves Σ_i w_i d_ik r_i^fit = 0. Hence:
- Δp_j = b/Hp, the shift the data still ask for;
- Δχ² = b²/Hp, the χ² it would gain;
- σ(Δp_j) = 1/√Hp.

Limiting cases:
- u = 0: the derivative is orthogonal to every track degree of freedom, Hp = Hu, and the data constrain the parameter
  fully.
- g_j a linear combination of the d_k: the shift is just a different track, Hp = 0, and the data cannot see it.

For a group, g_ij is the sum of the derivatives of its members, and the group is treated as one parameter.

b/Hp is the shift the data still ask for in that parameter alone, b²/Hp the χ² it would gain. Since they ignore the
correlations between parameters, they are a diagnostic, not a new solution:
- a free parameter of a converged fit has b ≈ 0 when summed over all the data of the fit and the prior term (`-con`), not
  per dataset. It is also not ≈ 0 if the parameter enters a constraint: the common mode of children constrained to a
  free parent is balanced by the Lagrange multiplier, so judge such modes by their group (below);
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
| `-pv no\|only` | process only the records without / with vertex labels |
| `-pt min,max` | process only the records all of whose tracks have an estimated pT in [min,max) GeV (below) |
| `-j` | threads, one input file each at a time |

The records store no track parameters. For `-pt` the pT of each track of a record is estimated from its multiple-scattering
kink rows (zero measurement, no global derivative): their errors are θ0·√((1−snp²)(1+tgl²)) and θ0·(1+tgl²), with
θ0 = 0.0136/p·√(x/X0)·(1+0.038 ln(x/X0)). tgl is taken from the pair with the smallest error ratio, and p from the median
over the first (ITS) kinks with nominal x/X0 (0.36 % IB, 0.9 % OB, times the path factor). The resolution is about 20 % per
track, the scale good to about 10 %. Tracks with fewer than 4 ITS kinks get no estimate, and records without any are
dropped. At B = 0 the kinks use the fixed `meanPtB0`, so the selection is meaningless there.

Records containing a `PVT/meanVertex` label (the multi-track vertex records, or tracks with a vertex constraint) are
counted separately from the others.

## Output

All the sums are additive over input files.

| file | content |
|---|---|
| `<out>_files.txt` | per input file: records, skipped; for records without and with vertex labels: accepted, rejected by χ², rejected by down-weight fraction, ⟨χ²/ndf⟩ of the accepted; records dropped by `-pt` |
| `<out>_pars.txt` | per label: volume, DOF, time slot, solution value, records, rows, b, Hu (w/o profiling), Hp, b/Hp, b²/Hp (cm, rad) |
| `<out>_groups.txt` | the same per group |
| `<out>_groupH.txt` | `# i name` lines, the vector b of the groups (line `b`), then the nG × nG profiled curvature matrix between groups (lines `H<i>`), H_ab = Σ w·g_a·g_b − u_aᵀA⁻¹u_b, whose diagonal is the Hp of `_groups.txt`; priors are not included |

Each parameter is in the frame of its volume's DOF set. For TRD chambers that is the sensor frame (`TX` radial, `TY` rφ,
`TZ` along the beam), while for TRD supermodules it is the TGeo frame, see [DOFConfig_rules.md](DOFConfig_rules.md#trd-alignable-match-strings).

## Reading the tables

### Columns of `_pars.txt` and `_groups.txt`

| column | meaning |
|---|---|
| `nRec`, `nRow` | accepted records and measurement rows with a non-zero derivative on the parameter (for a group, on any member) |
| `b` | Σ w·g·r over those rows: g is the derivative of the prediction (the group's derivatives are summed), r the residual after the solution and the refit of the local parameters, and w = (down-weight)/σ² |
| `Hu` | Σ w·g², the curvature with the track parameters kept fixed |
| `Hp` | Hu − uᵀA⁻¹u, the curvature that is left after the track parameters absorb what they can |
| `shift` | b/Hp, the change to *add* to the solution value of this parameter (or to every member of the group), in cm or rad, if it were the only thing to move |
| `dchi2` | b²/Hp, the χ² decrease that shift would bring, which is also (shift/σ)² with σ = 1/√Hp |

A row of zeros means no accepted record has a derivative on the parameter. Usually the label is missing from the
`-tree` file: give the tree of the *recording* run, since a fit's own tree omits its fixed DOFs.

### How to read them

1. **Significance.** σ = 1/√Hp is the statistical error of the shift. dchi2 is the significance squared, so on a
   converged fit and a statistically clean sample it is χ²(1)-distributed: values of a few are noise. b, Hu and Hp
   grow linearly with the number of records. A replay of k of the N files of a dataset therefore underestimates the
   full-sample dchi2 by about N/k, and overestimates σ by about √(N/k).
2. **Degeneracy with the tracks.** Hp/Hu is the fraction of the parameter that the track parameters cannot absorb.
   - Hp/Hu ≈ 1: a local motion that the tracks cannot mimic (in-plane chamber rotations, single-layer shifts).
   - Hp/Hu ≪ 1: a motion that is mostly a change of the track parameters, so the residuals constrain it only weakly
     (common rφ shifts of a whole stack or sector, radial shifts).
   - When Hp/Hu is small, a large dchi2 is still a real request from the data, but one that a small track-model
     error (energy loss, multiple scattering, field) can produce just as well.
3. **A free parameter of the fit itself** sits at the minimum of the fit's total χ², which includes every dataset and
   the priors. Its b for a *subset* of the data is therefore not zero. Instead, the b of all datasets, plus the
   prior term, add up to ≈ 0.
   - Datasets asking for the **same** shift: the fit is not converged, or a constraint or prior holds the parameter
     back. Give `-con` to include the priors; the mean-of-children constraints are not included, see 5.
   - Datasets asking for **opposite** shifts: tension. The parameter is a compromise between datasets that do not
     agree on the geometry, which points to a non-geometric effect correlated with the field or the sample
     (calibration, track model). Run all of B0, B+, B− and cosmics. An effect that changes sign with the polarity
     (B+ vs B−) and one that depends only on |B| (B0 vs B±) point to different causes.
4. **A fixed parameter** shows what freeing it would bring, on top of the solution. Freeing is worth it if the
   datasets agree on the shift and dchi2 is large.
5. **Joint fits from `_groupH.txt`.** Any mode built as a linear combination of groups, x = J·α (columns of J =
   coefficients of the groups), has b_α = Jᵀb and H_α = JᵀHJ. Its joint shift is H_α⁻¹b_α, with covariance H_α⁻¹, and
   Δχ² = b_αᵀH_α⁻¹b_α. The correlation of two group estimates is −H_ab/√(H_aa·H_bb). Examples:
   - a global rotation δφ of all TRD supermodules is TX = R·δφ, RY = δφ for every sector, with R the radius of the
     sector origin (329.6 cm for TRD, 385 cm for TOF);
   - an ITS rφ mode ∝ r² is layer RZ_k ∝ r_k.
   Fitting correlated modes jointly avoids reading one mode's request as another's.
6. **Overlapping groups and parameters are not independent.** Each shift is computed as if nothing else moved, so do
   not add the dchi2 of groups that share a mode. For example, the sum over the stacks of a supermodule, or over the
   layers, of the chamber rφ shifts is the supermodule rφ shift. Moreover, the mean-of-children constraint moves the
   all-chamber common mode into the parent. The informative quantities are therefore the *differences* between the
   chamber groups, e.g. layer to layer, not their common part. Compare that common part with the parent group
   (`TRD/sm..` `TX`).
7. **Rejections bias b.** Only the accepted records contribute. When a dataset loses a large fraction of its records
   (see `_files.txt`), its b describes the surviving, better-modelled tracks.
