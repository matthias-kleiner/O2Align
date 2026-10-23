# IB end-chip dy-vs-Z pattern: origin and proposed fix

Written for: whoever continues the ITS alignment tuning (knows O2Align, Millepede-II and checkResid).
State: 2026-10-24. Supersedes the conclusion "rphi-only, non-rigid term" of the 2026-10-22/23 version.

## 1. Summary

- The Lr0 end-chip excess in dy vs Z (page 3 of the checkResid PDF) comes from the Lr0 hits themselves.
  A bias of the Lr6→Lr1 track that is extrapolated to Lr0 does not explain it.
- The shape is that of a rigid **chip RX**: a tilt about the chip's rφ axis, so that the normal
  displacement grows along z. v8ns, v10 and v11 keep RX fixed for the IB chips (only TX, TY, TZ, RY are
  free), so nothing in those fits can absorb it.
- In the Mille data, freeing chip RX on top of the v8ns DOF set gains Δχ² ≈ 33k / 70k / 72k
  (B0 / B+ / B−) for 432 parameters. The fitted pattern agrees with the checkResid estimate.
- The pooled RY "endChipTilt" of `handout_endchip_dy_vs_z.md` rotates about the wrong axis: in the
  local frame RY is the in-plane rotation. The proposed configuration is `doc/ITS_RB_v13.json` (section 5).
- A separate, larger, non-rigid term remains open (section 6). It is not specific to the end chips.

## 2. checkResid with cluster size and snp (v8ns trees)

Trees: `/data/devAlg/solA2v8ns/checkResid_LHC25ab_562862_AlgRS_PV1.root` (pp, B0) and
`checkResid_LHC25ae_AlgRS_PV1.root` (OO, B+), with `points.clSize` and `points.snp` added.

- **No cluster-size dependence.** On Lr0 chip0 the dy-vs-local-z profile is the same for every
  `clSize` from 1 to more than 7. The two datasets agree.
- **snp is the hit position.** Within a chip, snp and the rφ hit position across the chip (v) are
  degenerate for primary tracks. In B0, snp minus the position-based angle y/x is a constant per sensor
  (±0.08, sign set by the half barrel) with only 0.005 rms. So a dependence on snp is a dependence on
  where the track crosses the chip, and primaries cannot separate an incidence effect from a position
  effect.
- **Shape on Lr0 chip0.** The residual behaves like dy ≈ A(u) + B(u)·v, where u is the local z pointing
  towards the stave end. B grows towards the outer end: dy runs from −2 to +25 µm across v at the outer
  end, but only from −10 to −2 µm at the inner end. Read as a normal displacement δX(u)·tanφ_inc, this
  gives δX ≈ 13 → 28 → 62 µm across chip0 (B+). That is mostly linear in u, i.e. a tilt.
- **z side.** The dz-based estimate of the same δX(u) rises the same way on chips 0–1, but at about
  half the size. Tracks reaching the outer end come from |zv| ≈ 10–13 cm with small |tanλ|, so the dz
  lever there is weak.

## 3. Hit displacement or track bias?

The Lr0 residual is an extrapolation of the Lr6→Lr1 fit, since the outward prediction holds only the
vertex, which gives no direction. Three tests decide whether the cause is in the Lr0 hits or in that
track.

1. **No dependence on pT or charge (B+).** On the chip0 outer end, dy0 is +11.6, +11.5, +10.1 and
   +11.6 µm for the four combinations of charge and pT (< 0.8 / > 0.8 GeV/c). An offset of the Lr1 or
   Lr2 hits would reach Lr0 with a weight that roughly halves between these pT bins: the toy below gives
   ∂dy0/∂δ1 = −1.45 at 0.6 GeV/c and −0.66 at 3 GeV/c. A bias from curvature or from the field model
   would flip with charge.
2. **The OB has almost no weight.** A toy of the checkResid leave-one-out Kalman combination
   (straight lines, or a helix with a loose curvature prior for B≠0, with multiple scattering and
   approximate material) gives the OB almost no weight in the IB residuals at these momenta. In practice
   "the Lr6→Lr1 track" means the Lr1/Lr2 hits.
3. **The vertex anchors the attribution.** The mean residuals of the same tracks (dyv, dy0, dy1, dy2)
   were solved for the hit offsets (δ0, δ1, δ2) with the toy response matrix.
   - Only offsets relative to the vertex can be determined, because a common shift of all points is
     invisible.
   - With a fixed vertex, the Lr0 end-chip hits come out displaced by +15…+40 µm at the outer end of
     chip0, and δ1, δ2 are smaller.
   - Forcing δ0 = 0 needs a "vertex offset" of −3 to −32 µm that depends on u and v of the Lr0 hit,
     even after averaging over all 12 staves. A real PV translation would vary like sin φ, average to
     zero over the azimuth, and could not know where the track crossed an Lr0 chip.

The vertex "vortex" for |zv| > 11 cm tracks (a φ-independent dyv of −20…−97 µm) is the same effect
seen from the vertex. Displacing the Lr0 hit drags the 7-point extrapolation to r = 0 by about −3 to −4
times dy0.

## 4. Mille fit

Data: `/data/devAlg/mille{B0,BPlus,BMinus}/mp2data_0.bin`, recorded with all six chip DOFs. The extended
accumulator adds per-chip columns beyond the rigid y/z copies of `macro/milleSplit.cxx`. Local
parameters are profiled out per record.

Δχ² over the rigid model with 6 DOFs per chip; each row has 432 parameters, which is also the gain
expected by chance:

| term | B0 | B+ | B− |
|---|---|---|---|
| y-only normal β (∝ tanφ_inc, y rows only) | 79.7k | 94.6k | 88.9k |
| in-plane rφ scale ∝ v (identical to β) | 79.7k | 94.5k | 88.9k |
| out-of-plane bow ∝ w² (y and z rows) | 10.4k | 19.0k | 18.4k |
| twist ∝ v·w (y and z rows) | 3.3k | 5.8k | 5.7k |
| in-plane bow ∝ w² | 1.3k | 1.6k | 1.6k |
| y-only normal × w | 1.3k | 2.4k | 2.3k |

The global fit therefore needs no non-rigid term that varies along z on the end chips. The z dependence
seen in checkResid is rigid.

The same data refitted on top of the v8ns set (TX, TY, TZ, RY per chip):

| added | parameters | B0 | B+ | B− |
|---|---|---|---|---|
| RX, all chips | 432 | 32.9k | 69.9k | 72.5k |
| RX, end chips only | 96 | 3.8k | 7.6k | 7.5k |
| RZ, all chips | 432 | 94.4k | 68.0k | 68.0k |
| RX + RZ | 864 | 125k | 136k | 138k |
| RX + RZ + β | 1296 | 205k | 230k | 227k |

Fitted RX, as the stave mean of the normal displacement at the chip's outer end. The pattern is the same
in all three fields; the stave rms is 10–15 µm on the end chips, and the error on the mean is about 1 µm:

| | chip0 | chip1 | chip2 | chips 3–7 | chip8 |
|---|---|---|---|---|---|
| Lr0, Mille | ~31 µm | ~15 µm | ~8 µm | 1–5 µm | ~0 |
| Lr0, checkResid δX(u) slope | ~25 µm | ~12 µm | ~6 µm | small | small |
| Lr1 / Lr2, Mille | ~11 / ~7 µm | | | | ~5 / ~3 µm |

The tilt is strongest at the z>0 end and grows chip by chip towards it, like a stave bending there. The
direction (lift outwards or inwards) has not been checked against the sign conventions.

## 5. Proposed configuration: `doc/ITS_RB_v13.json` (from v8 / v8ns)

| level | change | reason |
|---|---|---|
| IB chips | free RX and RZ (all 6 DOFs); priors 3 mrad on each, the v8 TX/TY/TZ priors kept | RX is the end-chip tilt; RZ is also strongly preferred |
| IB staves | free RX (prior 0.5 mrad); `pinChildrenMeanNoCross: ["RZ"]` | the chip tilts are coherent along a stave, so the mean RX per stave is not zero; the automatic constraint passes that mean to the stave. Stave RZ stays fixed, so the chips' mean roll is pinned |
| IB half barrels | `pinChildrenMeanNoCross: ["RX"]` | pins the collective stave tilt (the cone mode, radius linear in z) |
| removed | the stave `weightedMean` rule | v8ns used equal-weight constraints |

- **Running it:** the Mille data above carry all derivatives, so regenerate the steering files only
  (`--output MilleSteer`) and rerun pede.
- **Checks:** the constraint rank in `millepede.log`, the radial polar weak modes in `weakModes.json`,
  and the checkResid comparison with v8ns.
- **Leaner variant:** if the stave RX is weakly determined, free only chip RX/RZ and pin RX and RZ on
  the staves. The cost is about 5 µm distributed onto every chip.

## 6. Open item: the y-only normal term β

- **Size:** it is the largest single gain: 118–151k over the v8ns set, 80–95k over a fully rigid set.
- **Degeneracy:** it is exactly degenerate with an in-plane rφ scale across each chip, because v is
  equivalent to tanφ_inc for primaries.
- **Values:** the per-chip values scatter by 35–120 µm with a mean near zero. The rms is the same for
  every chip index of a layer and grows with the layer radius (Lr0 : Lr1 : Lr2 ≈ 77 : 100 : 123 µm in
  B0, about half of that in B±).
- **Origin:** not known. It could be an incidence-dependent cluster response, or a weak mode of the
  per-chip fit. No rigid DOF can absorb it.

## 7. Tools

- **In the repository:** `macro/milleSplit.cxx`, `macro/milleSplitModels.py` and `macro/sensDump.C`,
  documented in `doc/milleSplit.md`.
- **Session scratchpad only:**
  - the extended accumulator with bow/twist/shear columns (`milleSplit2.cxx`) and its fits
    (`fit3.py`, `fit4.py`);
  - the checkResid flattener with clSize/snp (`flatIB2.C`);
  - the Kalman response toys (`toy.py`, `toyB.py`);
  - the vertex-anchored inversions (`c1.py` – `c4.py`).

  Copy them into `macro/` if they are needed again.
