# IB end-chip dy-vs-Z pattern: origin and proposed fix

Written for: whoever continues the ITS alignment tuning (knows O2Align, Millepede-II and checkResid).
State: 2026-10-24. Supersedes the conclusion "rphi-only, non-rigid term" of the 2026-10-22/23 version:
that term was an artefact of the offline Mille fit (section 5).

## 1. Summary

- The Lr0 end-chip excess in dy vs Z (page 3 of the checkResid PDF) is reproduced by rigid motions of the
  IB chips. The IB Mille data are consistent with rigid chips: there is no rφ-only term.
- **Two ingredients are missing from v8ns/v10/v11:**
  - chip **RX**, a tilt about the chip's rφ axis, so that the normal displacement grows along z. Those
    configurations keep it fixed.
  - large **end-chip shifts**. On Lr0 chip0, the data want TY ≈ −37 µm (stave rms 24 µm) and TX ≈ +11 µm.
    The chip priors inherited from v8 (TY 5 µm, TX 10 µm) forbid this.
- **v13** freed RX and RZ but kept those priors. Pede then compensated on the end chips with an RX of the
  wrong sign, and the dy-vs-Z pattern did not improve.
- **v13b** (`doc/ITS_RB_v13b.json`) relaxes the end-chip priors. It is not tested yet.
- The pooled RY "endChipTilt" of `handout_endchip_dy_vs_z.md` rotates about the wrong axis: in the local
  frame RY is the in-plane rotation.

## 2. checkResid with cluster size and snp (v8ns trees)

Trees: `/data/devAlg/solA2v8ns/checkResid_LHC25ab_562862_AlgRS_PV1.root` (pp, B0) and
`checkResid_LHC25ae_AlgRS_PV1.root` (OO, B+), with `points.clSize` and `points.snp` added.

- **No cluster-size dependence.** On Lr0 chip0 the dy-vs-local-z profile is the same for every `clSize`
  from 1 to more than 7. The two datasets agree.
- **snp is the hit position.** Within a chip, snp and the rφ hit position across the chip (v) are
  degenerate for primary tracks. In B0, snp minus the position-based angle y/x is a constant per sensor
  (±0.08, sign set by the half barrel) with only 0.005 rms. That offset is the beam position, at
  Y ≈ −1.8 mm. So a dependence on snp is a dependence on where the track crosses the chip.
- **Shape on Lr0 chip0.** The residual behaves like dy ≈ A(u) + B(u)·v, where u is the local z pointing
  towards the stave end. B grows towards the outer end: dy runs from −2 to +25 µm across v at the outer
  end, but only from −10 to −2 µm at the inner end. Read as a normal displacement δX(u)·tanφ_inc, this
  gives δX ≈ 13 → 28 → 62 µm across chip0 (B+): a shift plus a tilt.
- **Fine binning (0.3 cm).** On the z>0 side, Lr0 shows a per-chip sawtooth: dy rises within each chip and
  drops at every chip boundary. Chip0 goes −6.5 → +13 µm, chip1 −2 → +2 and chip2 −1 → +1 µm. The z<0
  side has a weaker, mirrored version.

## 3. Hit displacement or track bias?

The Lr0 residual is an extrapolation of the Lr6→Lr1 fit, since the outward prediction holds only the
vertex, which gives no direction.

1. **pT and charge (B+).** On the chip0 outer end, dy0 is +11.6, +11.5, +10.1 and +11.6 µm for the four
   combinations of charge and pT (< 0.8 / > 0.8 GeV/c).
   - This excludes a bias from curvature or from the field model, and an offset that differs between Lr1
     and Lr2: the weight of such an offset at Lr0 changes about 2× between the pT bins.
   - It does **not** exclude a common shift of Lr1 and Lr2 relative to Lr0. That reaches Lr0 with a weight
     of about −1 at any momentum (−1.03 at 0.6 GeV/c, −1.07 at 3 GeV/c in the toy).
2. **The OB has almost no weight.** A toy of the checkResid leave-one-out Kalman combination (straight
   lines, or a helix with a loose curvature prior for B≠0, with multiple scattering and approximate
   material) gives the OB almost no weight in the IB residuals at these momenta.
3. **The vertex anchors the attribution.** The mean residuals of the same tracks (dyv, dy0, dy1, dy2) were
   solved for the hit offsets (δ0, δ1, δ2) with the toy response matrix.
   - Only offsets relative to the vertex can be determined.
   - With a fixed vertex, the Lr0 end-chip hits come out displaced by +15…+40 µm at the outer end of
     chip0, and δ1, δ2 are smaller.
   - Forcing δ0 = 0 needs a "vertex offset" that depends on where the track crosses the Lr0 chip, even
     after averaging over all 12 staves. No PV translation can do that.

The vertex "vortex" for |zv| > 11 cm tracks (a φ-independent dyv of −20…−97 µm) is the same effect seen
from the vertex. The Mille fit of section 4, which uses the vertex, attributes the effect to the end chips
themselves.

## 4. Mille fit (mean vertex free)

Data: `/data/devAlg/mille{B0,BPlus,BMinus}/mp2data_0.bin`, recorded with all six chip DOFs; the same files
as the v13 pede run.
- **Parameters:** the extended accumulator adds per-chip columns beyond the rigid y/z copies. The
  mean-vertex calibration of every time slot is free; it comes out at (+0.32, −1.82) mm, as in pede.
- **Other globals:** the other non-IB parameters stay at their recorded reference values.
- **Locals:** the local parameters are profiled out per record.

Δχ² (B0 / B+ / B−). The gain expected by chance equals the number of parameters.

| over the rigid model (6 DOFs per chip) | parameters | B0 | B+ | B− |
|---|---|---|---|---|
| full y/z split | 5184 | 4.6k | 10k | 9.6k |
| y-only normal β (∝ tanφ_inc, y rows only) | 432 | 2.3k | 6.3k | 6.3k |
| out-of-plane bow ∝ w² (y and z rows) | 432 | 10.2k | 18.7k | 18.2k |
| y-only normal × w | 432 | 0.6k | 1.0k | 0.9k |

| over the v8ns set (TX, TY, TZ, RY per chip) | parameters | B0 | B+ | B− |
|---|---|---|---|---|
| RX, all chips | 432 | 31.3k | 68.2k | 70.0k |
| RX, end chips only | 96 | 3.6k | 7.7k | 7.4k |
| RZ, all chips | 432 | 19.3k | 36.1k | 36.6k |
| RX + RZ | 864 | 48.2k | 101k | 103k |

- **Rigid chips are enough.** Only a modest out-of-plane bow remains beyond them.
- **The RZ pattern.** RZ shows a coherent twist along Lr0, from +1.5 mrad on chip0 to −1.3 mrad on chip8.
  Pede finds it too (section 6).

Lr0 chip0, stave mean (rms), all six chip DOFs free, B+:

| TX | TY | TZ | RX (normal shift at w = 1.5 cm) | RY (at 1.5 cm) | RZ |
|---|---|---|---|---|---|
| +11 (13) µm | −37 (24) µm | +1 (5) µm | −28 (15) µm | −6.5 (5) µm | +1.5 (1.1) mrad |

RX, as the stave mean of the normal displacement at w = 1.5 cm. The pattern is the same in all three
fields:

| | chip0 | chip1 | chip2 | chips 3–7 | chip8 |
|---|---|---|---|---|---|
| Lr0 | −30 µm | −14.5 µm | −7 µm | −1…+5 µm | ~0 to +1.5 µm |
| Lr1 / Lr2 | −10 / −7 µm | | | | −6 / −3 µm |

Chip8 has a large stave-to-stave scatter instead of a coherent shift (TY rms 34 µm, TZ rms 24 µm). The
direction (towards or away from the beam) has not been checked against the sign conventions.

## 5. The artefact: y-only normal term with the mean vertex fixed

A first version of the accumulator kept every non-IB global fixed, including the mean vertex. The vertex
records measure the mean vertex relative to a reference at (0, 0), so each vertex was pulled by about
1.8 mm (~40σ), purely in the transverse plane. That bends the tracks in rφ only. The fit could absorb it
only with a y-only chip term, and produced:
- β: Δχ² 80–150k;
- full y/z split: 88–101k;
- a per-stave first harmonic in the stave azimuth, ∝ R·cos α, twice as large in B0 as in B±.

All of this disappears with the mean vertex free (section 4). `macro/milleSplit.cxx` now frees the mean
vertex, and the rule is noted in `CLAUDE.md`.

## 6. Configurations tested and proposed

### v13 (`doc/ITS_RB_v13.json`, result in `/data/devAlg/solA2v13`)

Changes from v8 / v8ns:

| level | change |
|---|---|
| IB chips | free RX and RZ (all 6 DOFs), priors 3 mrad; the v8 priors TX 10, TY 5, TZ 10 µm kept |
| IB staves | free RX (prior 0.5 mrad); `pinChildrenMeanNoCross: ["RZ"]` |
| IB half barrels | `pinChildrenMeanNoCross: ["RX"]` |
| removed | the stave `weightedMean` rule (v8ns used equal-weight constraints) |

**Result.**
- dz vs Z of Lr1/Lr2 and the vertex dy at Z > 6 cm improved.
- The Lr0 dy vs Z did not improve: per-bin medians at the z>0 end went +13 → +10 µm, at the z<0 end
  +6 → +10 µm (worse), and the sawtooth remains.

**Why.** For Lr0 chips 1–7, the v13 values agree with the free Mille fit in every DOF to within a few µm,
including the RZ twist. On chip0, the priors block the shift the data want, and pede compensates with an
RX of the wrong sign:

| Lr0 chip0, stave mean | TX | TY | RX (at 1.5 cm) | RY (at 1.5 cm) |
|---|---|---|---|---|
| Mille fit, mean vertex free, no priors | +11 µm | −37 µm | −28 µm | −6.5 µm |
| pede v13 | −2 µm | +0.6 µm | +13 µm | +3 µm |

### v13b (`doc/ITS_RB_v13b.json`, not tested)

It is v13 plus one rule after the chip prior rule. The rule relaxes the priors of IB chips 0 and 8 to
TX 30, TY 50 and TZ 30 µm, keeping RX and RZ at 3 mrad.

Checks when it is run:
- the Lr0 dy vs Z at both ends and the per-chip sawtooth;
- whether Lr0 chip0 TY comes out near −37 µm and RX near −28 µm (at 1.5 cm);
- the constraint rank in `millepede.log` and the weak modes, since the end-chip shifts are now loosely
  constrained.

The Mille data above carry all derivatives: regenerate the steering files with `--output MilleSteer` and
rerun pede.

## 7. Open points

- **Lr0 sawtooth.** The per-chip sawtooth seen in checkResid should be removed by chip RX/RY. If it
  survives v13b, compare checkResid with a residual computed from the full vertex + Lr1..Lr6 fit, which is
  what pede effectively sees.
- **Out-of-plane bow.** The bow ∝ w² (Δχ² 10–19k for 432 parameters) is the only non-rigid shape with a
  significant gain. It cannot be fitted with the current DOF sets.

## 8. Tools

- **In the repository:** `macro/milleSplit.cxx`, `macro/milleSplitModels.py` and `macro/sensDump.C`,
  documented in `doc/milleSplit.md`. The mean vertex is free since 2026-10-24.
- **Session scratchpad only:**
  - the extended accumulator with bow/twist/shear columns and a free mean vertex (`milleSplit3.cxx`) and
    its fits (`fit6.py` – `fit8.py`);
  - the checkResid flattener with clSize/snp (`flatIB2.C`);
  - the Kalman response toys (`toy.py`, `toyB.py`);
  - the vertex-anchored inversions (`c1.py` – `c4.py`).

  Copy them into `macro/` if they are needed again.
