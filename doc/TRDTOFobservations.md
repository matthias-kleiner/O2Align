# TRD/TOF alignment with the ITS fixed: observations

Date: 2026-11-20

The main conclusions of the TRD/TOF iterations on `Mille_IBEdgeS4R` (recorded with `ITSnochipTRDTOF_RB.json` on the IBEdge
geometry, `applyTRDshift=false`, `useMultiTrackPVConstraint=2`), and of the fits `solB0A1`, `solB1A1`, `solB1A2`. The
chronology and the details are in [history.md](history.md). The tools are `o2-dev-mille-replay` ([MilleReplay.md](MilleReplay.md))
and the `gbltrack` debug trees ([DebugInfo.md](DebugInfo.md)).

Data: one run per collision dataset (B0 562862, B+ 562837, B− 564644, plus the `Bis` files), and 123 cosmic runs
(LHC25aa, LHC26aa–aj). The tracks are ITS + TRD and/or TOF, without TPC points. Shifts are given in mm at the sector
origin (TRD R = 329.6 cm, TOF R = 385 cm); TRD/TOF sector frames: TX rφ, TY beam, TZ radial; TRD chamber (sensor) frame:
TX radial, TY rφ, TZ beam.

## 1. TRD pad shift: seen by all datasets, absorbed by the sector rφ

- Without the old pad-to-rφ bug (`applyTRDshift=false`), the IBEdge TRD alignment is off by about one pad. Since the pad
  width grows with the radius, this is close to a global rotation of TRD (sector TX = R·δφ, RY = δφ).
- At the IBEdge reference, a joint fit of the TRD and TOF global rotations gives the TRD rotation +4.9 / +5.1 / +5.3 mm
  for B0 / B+ / B− (δφ ≈ 1.5–1.6 mrad) and +6.2 mm for cosmics. All datasets see it.
- After `solB1A2` the TRD rotation is accepted by all datasets (+0.03…+0.28 mm), and run by run in cosmics
  (χ²/ndf 82/83).
- TRD sector TX and RY are well separated (correlation 0.2–0.3, layer lever arm ±32 cm).
- Caveat for single-group reading: the TRD and TOF global rotations are correlated by 0.5–0.8, so a single-group TRD rφ
  shift carries the TOF tension (section 2). Use joint fits (`_groupH.txt`).

## 2. TOF rφ: B on/off tension from a non-uniform in-pad residual distribution

### The tension

- At the IBEdge reference and at `solB1A2`, the TOF global rotation relative to TRD is requested with opposite signs by
  the B = 0 and B ≠ 0 data:

  | | B0 | B+ | B− | cosmics |
  |---|---|---|---|---|
  | IBEdge | −0.74 | +1.87 | +1.62 | +0.52 |
  | `solB1A2` | −0.83 | +1.45 | +1.14 | −0.26 |

- It is polarity-even, survives removing the composed PV records, and is already present with the same size (2.0–2.6 mm)
  at IBEdge, i.e. it was not created by the TRD/TOF fit.

### Not an ITS curvature weak mode

- An ITS rφ mode ∝ r² (a charge-independent curvature at B ≠ 0) moves the B ≠ 0 requests of TRD and TOF in the ratio
  (r_TRD/r_TOF)² ≈ 0.74, and the B0 ones by much less.
- Closing the gap would need d₆ (its rφ at ITS L6) ≈ +3…+5 µm for TRD and −28…−40 µm for TOF. No single value does both.
- The data themselves ask for d₆ = −11 µm (B0) and −16 µm (cosmics).

### Cause: the TOF residual shape at B ≠ 0

- TOF measurement: the centre of the main contributing pad (σY 7.5 mm, 25 mm pads). The residual of ITS+TRD+TOF tracks is
  therefore a 25 mm box smeared by the prediction error (3.7 mm at B0, 4.6 mm at B+).
- At B0 the box is flat (−4.9…+21 mm, centre +8.2 at IBEdge). At B+ and B− the right edge is at the same place within
  ≈ 0.5 mm, but the box carries an extra narrow peak at 0…+5 mm:

  | ITS+TRD+TOF, TOF dyUnb | median | centre of 5–95 % | skew (84−50)−(50−16) |
  |---|---|---|---|
  | B0 | +8.4 | +8.2 | +0.0 |
  | B+ | +5.1 | +7.4 | +3.2 |
  | B− | +5.5 | +7.95 | +3.4 |
  | cosmics | +7.2 | +7.2 | +0.5 |

- The skew does not depend on the incidence angle: it is +2…+4 mm already at |snp| < 0.05, where B0 has +0.1. It is the
  same for both polarities, absent in cosmics (≤ 2 mm run by run), and the B0 and B+ runs are only 25 runs apart. All
  TOF clusters have the same error, so cluster types cannot be told apart in the Mille data.
- With pad-centre measurements and a 4 mm prediction error, a uniformly illuminated pad can only produce the flat B0 box.
  The B ≠ 0 peak means that the predicted crossings are not uniform within the matched pads, or that the main
  contributing channel is not always the pad crossed by the track.
- The pad positions (box edges) agree between B on and off. A least-squares fit with σ ≈ 7.5 mm turns the shape change
  into a 2–3 mm rφ shift. The tension is therefore not a geometry difference.
- At B = 0, ITS+TOF tracks compare TOF directly with ITS: B0 +0.2…+0.4 mm, cosmics −0.7 mm at IBEdge. The TOF global
  rotation of IBEdge is good within ≈ 0.7 mm.

**Open:** the origin of the B ≠ 0 peak, in TOF clusterization or TOF matching. To check on the reconstruction side, e.g.
in the TOF matching residuals (dx) of runs 562837 vs 562862.

### Track composition at B ≠ 0

- 57 % of the TOF tracks (B+) are ITS+TOF only. At B ≠ 0 their curvature is fixed by the TOF point itself, so a TOF rφ
  shift is nearly degenerate with it; their TOF rφ information comes only from the ITS sagitta.
- Their unbiased TOF residual binned in the KF q/pT is self-correlated (a ±6 mm wiggle around q/pT = 0): q/pT is not a
  usable binning variable for TOF.
- `dyPre` is not an ITS-only residual. Its seed is the inward KF refit over all the points of the track, so at TOF it is
  biased towards zero by the point itself and 2–3 times wider than `dyUnb`.

## 3. TRD radial

### Common radial: drift time-zero, not geometry

- Requests of the common TRD radial shift: cosmics +4.6 mm in every period (LHC25aa +4.9, LHC26aa +4.4), with about 1 mm
  of run-to-run scatter beyond the errors (χ²/ndf 336/83). Collisions: B0 +1.3, B+ +0.1, B− −2.3 mm, three different runs.
- A drift-time offset δt moves every tracklet by vd·δt along the drift (radial) direction, which a common radial shift
  reproduces exactly: 4.6 mm ≈ 0.3 µs at vd ≈ 1.5 cm/µs. The common radial is read as the time-zero of each run (cosmic
  trigger vs collisions).
- The "polarity-odd" TRD radial seen earlier is explained by B0, B+ and B− being three different runs.
- Up-weighting the cosmics would write their time-zero into the geometry.

### Chamber-level radial: geometry

- After removing the sector and layer means, the per-chamber radial requests correlate between all datasets. The common
  part is 0.5–0.7 mm rms between collisions and cosmics, and cosmics 2025 vs 2026 (independent samples) correlate by 0.73
  (1.0 mm rms).
- Per track, B± carry ≈ 6.5× the radial information of B0 or cosmics (⟨tan²ψ⟩ 0.065 vs 0.010). With all the files, the
  chamber radial error is 0.22–0.28 mm at B± and 0.43 mm at B0.

### Layer-common radial: real but small; B± low-pT tracks distort it

- Layer deviations from the common mean: cosmics +0.4 (L0) → −0.4 mm (L5), the same in 2025 and 2026.
- B± replay: L0 −1.0/−1.5 → L5 +3.0/+3.6 mm (χ² ≈ 13–14k). This pattern depends on pT. In the debug trees, B± tracks with
  |q/pT| < 0.7 give +0.5…+0.8 → −0.5…−0.9 mm, like the cosmics; tracks with |q/pT| 1.0–1.5 give −1.4…−2.0 → +3.3…+4.1 mm.
  The low-pT tracks carry most of the bending lever arm and dominate. Geometry cannot depend on pT, so this is a
  momentum-dependent track-model effect inside TRD (multiple scattering, energy loss or curvature). Its cause is open.
- The stack radial pattern is not consistent between datasets (st2 −1.1…−4.4 mm, the others vary).

## 4. TOF radial and TOF strips

- TOF common radial in the joint fit: B0 −0.2 ± 0.2, B+ +0.1 ± 0.2, B− +0.0 ± 0.2, cosmics +0.2 ± 0.04 mm. All datasets
  agree; nothing to fix.
- TOF strips: only a ≈ 0.7 mm z pattern; kept fixed.

## 5. Track model and rejections

- Freeing the mean vertex reduced the pede rejections from 28.4 % (`solB0A1`) to 22.2 % (`solB1A1`). The rest are mostly
  B0 composed records with a track softer than `meanPtB0` = 0.7 GeV, whose multiple scattering is underestimated. At
  `solB1A2`, 56 % of the B0 PV records and 22 % of the others are χ²-rejected, against 16 % and 5 % at B+.
- The B0 composed records amplify the B0 requests (TRD rφ +0.63 → +0.38 mm and TRD radial +1.3 → +0.7 mm without them) but
  do not create them.

## 6. Configuration for the next recording (`extraClsErrYTOF=0.5`)

[ITSEnvTRDTOFch3_RB.json](ITSEnvTRDTOFch3_RB.json):

| volume | DOFs | priors / constraints | reason |
|---|---|---|---|
| ITS | fixed | | |
| mean vertex | free | | (section 5) |
| TRD sectors | TX, TY, RX, RY, RZ free; TZ fixed | `pinChildrenMean` `TZ` | rφ absorbs the pad shift; common radial = time-zero (3) |
| TRD chambers | TX, TY, RX | 0 ± 1 mm, 1 mm, 1 mrad | chamber radial pattern is geometry (3) |
| TRD layers L0–L4 | | `groupModes` `radialScaling` per layer (L5 follows from the sector pins) | B± low-pT layer pattern (3) |
| TOF sectors | all | `curl` (no common rotation about the beam); TX 0 ± 1 mm, TY/TZ 0 ± 2 mm, rotations 1 mrad | B ≠ 0 TOF residual shape (2); 18 × 1 mm priors would allow only 0.24 mm on the common rotation, while the data measure it to ≈ 25 µm |
| TOF strips | fixed | | (4) |

- The per-layer equations need `AlignParams.groupModesRelToRoots=false` when the steering files are written: with the
  default, they are relative to the chambers' own free motion and are skipped as empty.
- Alternative: record the B± collisions with `AlignParams.minPt` ≈ 1.5 GeV (default 0.4). The high-pT B± tracks see the
  same layer pattern as the cosmics and still have 2–4× the B0 radial information per track. Then drop the per-layer
  `groupModes` and keep only the sector pin. 1.0 GeV is not enough: the |q/pT| 0.7–1.0 bin already shows the low-pT
  pattern.
- Expected in `mp2con.txt`: 18 TRD sector `TZ` pins, 5 `Group mode … over TRD-L*` equations and one TOF `curl` equation.

## 7. `solB1A3` (ch3 config, recording with minPt 1.2 and LUT corrections)

- **Rejections.** B0 is fixed (6.5 / 8.7 % for single / PV records, was 22 / 56 %). B± PV records went up to 24 %, and
  cosmics to 12 % with `meanPtB0Cosm=6` (4 % with 2).
- **Chamber radial excursions are a constraint artefact.** Up to −38 mm, 6.8 % of the chambers beyond 3 mm, while all
  datasets ask to undo them. The equal-weight sector pins and per-layer equations put the same Lagrange force on every
  chamber, and the weakly measured chambers restore the mean.
- **The radial requests of B = 0 and B ≠ 0 have opposite layer patterns,** also with `minPt=1.2`: B0 and cosmics want
  +1 → −2 mm (L0 → L5), B± −2.3 → +3.3 mm. Per chamber, B± anti-correlate with the cosmics, while the two cosmic years
  agree. TRD radial information at B ≠ 0 is therefore affected by a field-dependent effect, like TOF rφ.
- **TRD sector RY and the layer-linear chamber TY form a nearly unconstrained pair** (−3 mrad and −0.9 → +1.2 mm).
  B = 0 data reject it.
- **Accepted by all datasets:** the TRD rotation (pad shift), the TRD sector z offsets, the chamber in-plane RX. The TOF
  `curl` holds the TOF rotation at its B = 0 value; the B ≠ 0 data still ask for about +2 mm.
- **ITS.** All datasets ask for an ITS rφ ∝ r² mode of d₆ ≈ −9…−17 µm at L6.

## 8. Open questions

1. The origin of the narrow peak in the B ≠ 0 TOF residual (TOF clusterization or matching).
2. The momentum-dependent TRD layer pattern at B ≠ 0 (multiple scattering, energy loss or curvature in TRD).
3. A per-run TRD time-zero (x-offset) calibration DOF, time-sliced like the mean vertex. The common radial would stay
   degenerate with the average time-zero and would have to be anchored by collisions with a calibrated time-zero.
4. The cosmic track model (`meanPtB0Cosm` = 2 GeV) and the B0 soft-track multiple scattering, which drive the rejections.
