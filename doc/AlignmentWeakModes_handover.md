# ITS Rigid-Body DOF Config: Weak-Mode Investigation (handover)

Investigation of why incremental Millepede alignment rounds on ITS (tests in
`/data/devAlg/tests/{ITSBarFix,StavesOnly,savStavesOnlyPinHB}`, data from runs LHC25ae
(B=0.5T) and LHC25ab/562862 (B=0)) showed a large pede chi2/ndf improvement (factor ~5)
while `checkResid` residuals either developed new systematic patterns or got worse
instead of better. All three tests derive from `doc/ITS_RB_from_hstaves_IB.json`, the
config used to produce the Mille data, by changing which DOFs are free at the stave/
half-stave/chip levels and whether `pinChildrenMean` is used.

## Summary of findings

### 1. `ITSBarFix` — chip-level freedom under a fully-free stave (pages 3-5 of
   `checkResid_hman_ITSBarFix_fine.pdf`: strong pattern vs Z, correlated with chip, on
   IB staves only)

Config: IB stave free in all 6 DOF (`TX,TY,TZ,RX,RY,RZ`); IB chip free in `TX,TY,TZ`
only (no rotation), directly below it.

- `RX`/`RY` (tilts about the axes transverse to the stave) transport to a child at
  offset *z* along the stave as a translation proportional to *z*. A chip with
  translation but no rotation freedom is the only thing available to absorb that
  lever-arm term, so a genuine stave tilt and a *linear-in-chip-index* pattern across
  the 9 IB chip `TY` values are (numerically) nearly indistinguishable.
- `Volume::writeChildrenMeanConstraints` only auto-constrains the mean (0th moment) of
  free children to the parent — it does not resolve this position-weighted
  (1st-moment/dipole) ambiguity.
- Confirmed directly in `millepede_labeled.res`: chip `TY` ramps smoothly and almost
  linearly across chip index on every IB stave checked (e.g. Stave0/Lr0: +114 → −68 µm
  chip0→chip8), with slope of the same order as that stave's own fitted `RX`. `TX`/`TZ`
  stay small and non-systematic. Only IB shows this (only IB chips got this DOF rule),
  matching the "chip-correlated, Z-dependent, IB-only" symptom exactly.

### 2. `StavesOnly` — half-barrel fixed, stave free, no `pinChildrenMean`
   (OB residuals get *worse* than the initial geometry)

Config: only the stave level is free (IB chip and OB half-stave now fixed); everything
above the stave (half-barrel, envelope) is fixed.

- A fixed parent with free descendants and **no** `pinChildrenMean` writes **no
  constraint at all** by design (`mp2con.txt` was literally 0 bytes). The common
  (mean) rigid-body mode of all staves under a given half-barrel is left fully
  unconstrained — a true near-null direction of the global matrix.
- Confirmed in `millepede_labeled.res`: the mean stave `TX` (≈ radial) and its
  stave-to-stave RMS both grow monotonically with layer radius, from ≈ −21 µm (Lr0) to
  **−537 µm (Lr6)**, with RMS up to 2.5 mm at Lr6 — the classic signature of an
  unconstrained/weakly-observable global mode absorbing noise, worse at larger R
  because the same angular ambiguity maps to a bigger absolute number there.

### 3. `savStavesOnlyPinHB` — same as above + `pinChildrenMean:"all"` at the
   half-barrel (little improvement in `checkResid`, chi2 still drops ×5)

- The pin works: `mp2con.txt` is no longer empty, the RZ mean is now exactly 0 per
  layer, and the pathological growth-with-radius bias from test 2 is gone (Lr6 TX mean
  −537 → −81 µm; RMS 2.5 mm → 64 µm). The fit is now well-posed.
- But `checkResid` barely improves. Reason: a single rigid stave transform is the
  *only* freedom available, and it can only capture the stave-average component of
  whatever is really misaligned. Test 1 already showed the real misalignment has
  genuine sub-stave (chip-level) structure that a per-stave rigid body cannot
  represent — it's a basis mismatch, not a conditioning problem.
- **Why residuals get worse than the (already decent) initial geometry, specifically
  at OB:** pede's own chi2/ndf is an *in-sample* metric — computed from the very same
  ~42M equations (`NEGB` in `millepede.log`) used to fit the 1182 free parameters
  (`NVGB`), so adding free parameters mechanically reduces it whether or not they
  reflect real geometry. The stave-only basis is still poorly observable at OB (fewer,
  larger staves ⇒ larger lever arms, larger correlated errors), so each stave's "best
  fit" correction is dominated by noise of the same order as (or larger than) any real
  residual misalignment already present. Applying that noisy correction on top of an
  already-good OB geometry displaces previously-correct staves by a comparable random
  amount — pede's chi2 improves because it is explaining away in-sample fluctuations,
  while the independent, leave-one-out `checkResid` estimator exposes that the actual
  geometry got worse.

## Matching engine: `fnmatch` replaced by `std::regex`

`Volume::applyDOFConfig` (`src/Volume.cxx`) matched `"match"` against each symbolic
name with POSIX `fnmatch` (glob syntax: `?`, `*`, `[...]`). This has been replaced by
`std::regex` (ECMAScript grammar): `"match"` is now a real regex, matched against the
*whole* symbolic name with `std::regex_match` (not a substring search), with each
pattern also tried with an implicit leading `.*` — the same "also works as a bare
suffix" convenience `fnmatch`'s implicit leading `*` gave, just expressed in regex.
Patterns are compiled once per rule before the volume traversal, not once per
`(rule, volume)` pair.

This is a breaking change to every existing DOF config JSON: glob's `?` (exactly one
char) becomes `.` (regex: any one char), and glob's `*` (any run of chars) becomes
`.*`. Bracket expressions (`[0-9]`, `[3-6]`, …) are unchanged — they mean the same
thing in both.

One real improvement falls out of this, beyond just syntax: a numeric ID field should
be written `[0-9]+` (one or more digits) rather than `.*`. Unlike `.*`, `[0-9]+`
*cannot* match `/`, so it is immune to the trailing-wildcard hazard described below —
it can never accidentally swallow a child level — and it matches any number of digits
in one rule, which removes the whole one-/two-digit OB stave-ID duplication (and its
bug, see next section) in every rule that has one, not just the half-stave rule.

## Proposal: revised `doc/ITS_RB_from_hstaves_IB.json`

The file that produced the Mille data frees **four nested levels** (envelope,
half-barrel, stave, chip) all in 6 DOF, anchored only by a soft 100 µm prior on the
envelope, with no `pinChildrenMean` anywhere, plus two bugs: the half-stave rule only
matched two-digit OB stave IDs (staves 0–9 got no half-stave freedom at all), and the
chip rule's layer range (`[0-3]`) was a harmless but misleading no-op for layer 3.

Both bugs disappear once every numeric ID field is written `[0-9]+`: there is no
longer a one-digit/two-digit rule pair to go out of sync, for the stave level *or* the
half-stave level, and the stray layer-3 chip pattern is simply dropped (its ID range
corrected to `[0-2]`, where it unambiguously means IB only).

Proposed replacement, combining the lessons above (fewest weak modes, still enough
free DOFs to capture the chip-scale structure seen in test 1), now in regex:

```json
[
  { "match": "ITS", "rigidBody": "fixed", "pinChildrenMean": ["TX", "TY", "TZ", "RZ"] },

  { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-1]", "rigidBody": ["TX", "TY", "TZ", "RZ"] },
  { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-1]", "measurement": { "TX": [0.0, 0.003], "TY": [0.0, 0.003], "TZ": [0.0, 0.003] } },

  { "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]", "rigidBody": "fixed", "pinChildrenMean": "all" },

  { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+", "rigidBody": ["TX", "TY", "TZ", "RZ"] },

  { "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+", "rigidBody": "all" },

  { "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+/ITSUHalfStave[0-9]+", "rigidBody": "all" },

  { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+/ITSUChip[0-9]+", "rigidBody": ["TX", "TY", "TZ"] }
]
```

Eight rules now (two of them the IB half-barrel's DOF definition and its measurement,
which — per the applyDOFConfig rule — must be separate rules, the measurement always
following the one that frees the DOFs).

Changes and rationale:

- **`ITS` fixed, with `pinChildrenMean` moved up here** from the half-barrel level.
  Overall ITS position/orientation belongs to global (GLO) alignment, so the envelope
  itself stays fixed, but since the IB half-barrels below it are now free (next
  bullet), the envelope needs its own pin to avoid reopening the `StavesOnly`-style
  unconstrained common mode, this time one level higher. `TX,TY,TZ,RZ` is enough:
  `RX,RY` have no free contributor from the IB side (the half-barrel itself doesn't
  free them) to be degenerate with, and leaving them out of the pin keeps the OB
  staves' own `RX,RY` (looked through the fixed-and-already-self-pinned OB
  half-barrel) from being pulled into one large, awkward cross-half-barrel
  constraint.
- **IB half-barrel freed in `TX,TY,TZ,RZ`, anchored by a 30 µm soft `measurement`
  prior instead of fixed + hard `pinChildrenMean`.** This is a direct response to a
  real physical suspicion, not a weak-mode fix: the two IB half-shells meet at a
  mechanical gap, and fixing them (as the original draft did) *structurally forbids*
  the fit from ever recovering a genuine common displacement of one half-shell
  relative to the other — `pinChildrenMean` at the half-barrel would have forced the
  mean of its own staves to exactly zero relative to the (assumed-perfect) fixed
  half-barrel pose, zeroing out precisely the mode a gap-closure error would produce.
  A soft prior (value 0, sigma 0.003 cm = 30 µm, the quoted uncertainty on the gap)
  keeps the same regularizing role — it still prevents this mode from blowing up the
  way the unconstrained `StavesOnly` common mode did — but *lets the data move it*
  if the evidence supports an offset beyond that tolerance, rather than hard-zeroing
  it. `RZ` is freed alongside the translations (a relative half-shell twist is as
  plausible as a pure offset) but left without its own measurement — cm-scale sigma
  doesn't translate to a rotation, and it's already regularized indirectly by the
  envelope-level pin and by the automatic stave/half-barrel mean constraint (both
  being free DOFs of matching name). Before spending effort on this, check
  `checkResid`'s dy/dz vs phi specifically **near phi ≈ 0/2π and ≈ π**, where the two
  half-barrels physically meet — a localized step there (rather than a smooth
  stave-to-stave pattern) is the expected fingerprint of this defect and would
  confirm it's worth carrying.
- Note there's no intermediate "layer" volume to pin a half-barrel *pair* against:
  per `CLAUDE.md`, the layer number is a path component, not a real alignable, so
  `ITS` (the envelope) is the *direct* parent of all 14 half-barrels across all 7
  layers at once. The envelope-level pin above therefore ties the aggregate mean
  across all of them (mostly the 6 IB ones, since the OB ones already net to ~0
  through their own pin), not a clean per-layer pairing — good enough to remove the
  degeneracy, but it's the per-half-barrel measurement, not this pin, that's doing
  the real job of letting an individual IB half-barrel's relative offset through.
- **OB half-barrel unchanged (fixed + `pinChildrenMean:"all"`)**: the gap suspicion is
  specifically about the IB half-shells; no comparable evidence or mechanical reason
  has come up for OB, so it keeps the simpler, more conservative treatment.
- **IB stave restricted to `TX,TY,TZ,RZ`** (drops `RX,RY`). Directly removes the
  degenerate channel identified in test 1 — the stave tilts whose lever-arm transport
  aliases with the per-chip translation ramp — without growing the parameter count
  (the alternative, giving the chip matching rotational freedom, would add 3 more
  parameters per chip that are barely observable on a ~mm-scale sensor at near-normal
  incidence, i.e. trade one weak mode for several worse ones). Cost: a genuine
  whole-stave tilt is no longer modeled here; if needed later, add it back via an
  external/survey prior rather than fitting it degenerately against chip positions.
- **IB chip kept at `TX,TY,TZ`** — this granularity is what actually improved the fit
  in test 1; the problem was its parent's DOF set, not its own.
- **OB stave/half-stave left at `all`/`all`** (matched DOF sets on both sides is the
  well-posed, documented pattern; no evidence so far implicates this level).
- **OB stave and half-stave ID fields written `[0-9]+`** (bug fix: the previous
  two-digit-only half-stave rule left staves 0-9 with no half-stave freedom at all).
  `[0-9]+` can't match `/`, so — unlike a bare `*`/`.*` — it's safe at the stave level
  too: it can never spill over into module/chip the way an unanchored wildcard would.
- **Chip rule layer range corrected to `[0-2]`** (was `[0-3]`, a no-op for layer 3 but
  misleading).

## A superset config to free every RB DOF (for `MilleData` production)

Separately from the restricted config above (meant for the final, constrained fit),
the earlier discussion on `MilleData` reuse concluded it's worth generating the Mille
binary once with every DOF free at every level, then exploring fixed/free/constraint
combinations downstream purely at the `pede`/`MilleSteer` stage (see the standalone
`MilleSteer` mode in `AlignmentSpec.h`). The simplest version of that superset is a
single catch-all rule — `rigidBody` is silently skipped wherever it doesn't apply
(the root, envelopes without their own geometry), so this is safe even though it's
unconditional:

```json
[ { "match": ".*", "rigidBody": "all" } ]
```

If the superset should be scoped to ITS only (leaving TRD/TOF/TPC/mean-vertex etc.
untouched), use one rule per level instead, covering both the IB chain (envelope →
half-barrel → stave → chip) and the OB chain (envelope → half-barrel → stave →
half-stave → module → chip):

```json
[
  { "match": "ITS", "rigidBody": "all" },
  { "match": "ITS/ITSULayer[0-6]/ITSUHalfBarrel[0-1]", "rigidBody": "all" },
  { "match": "ITS/ITSULayer[0-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+", "rigidBody": "all" },

  { "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+/ITSUHalfStave[0-9]+", "rigidBody": "all" },
  { "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+/ITSUHalfStave[0-9]+/ITSUModule[0-9]+", "rigidBody": "all" },
  { "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+/ITSUHalfStave[0-9]+/ITSUModule[0-9]+/ITSUChip[0-9]+", "rigidBody": "all" },

  { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+/ITSUChip[0-9]+", "rigidBody": "all" }
]
```

Note this is a config for `MilleData` generation, not for a fit to actually run with —
every chip rotation it frees is expected to be a weak/unobservable mode on its own
(small sensor, near-normal incidence); the point is only to make the corresponding
derivative columns exist in the binary so later tests can fix them (or constrain
them) cheaply, without a new `MilleData` pass, rather than to ever solve for them
directly.

## Running a restricted config over all-free Mille data: fix everything not mentioned

A restricted config applied to data generated with every DOF free has a trap:
`Volume::writeParameters` only writes parameter lines for volumes that *have* a
rigid-body DOF set. A volume matched by no rule has none, so its labels are absent from
`mp2param.txt` while their derivatives are still in the binary, and pede treats an
unlisted label as a variable parameter. In `ITS_RB_v0.json` this hits the OB modules
and OB chips (free in the all-free data, never mentioned).

`doc/ITS_RB_v0_explicit.json` is v0 with a catch-all `{"match": ".*", "rigidBody":
"fixed"}` as its **first** rule. Rules apply in order and a later match replaces the DOF
set, so the free set is exactly v0's; every volume the later rules don't mention now
carries an all-fixed DOF set and is written as `F` lines (presigma −1), which pede
honours. The rule is silently skipped for volumes that aren't rigid-body alignable (the
pattern contains `*`, so no warning). It stays correct whatever the superset freed,
unlike enumerating the missing levels (OB module, OB chip) by hand.

Check after regenerating the steering files: every label with derivatives in the data
should appear in `mp2param.txt`, e.g. compare its label count to the one in the
all-free run, and `pede.log` should report the expected `NVGB` (number of variable
parameters).

## 4. `solA1` — `ITS_RB_v0.json` actually run: what the fit really did

`/data/devAlg/tests/solA1` is the first *valid* run of the restricted config (`solA0` is void — it
is the pre-fix run, `NVGB = 146016`, i.e. pede varied every unlisted label). In `solA1`
`NVGB = 4134`, `NCGB = 1084`, all 157872 labels appear in `mp2param.txt` (4134 `V` + 153738 `F`),
and the constraint product matrix has full rank 1084. So the explicit-fixing recipe of the previous
section works, and the `std::regex` migration works in production (`.*`, `[0-2]`, `[0-9]+` all
matched as intended, including both one- and two-digit OB stave IDs: all 48 Lr6 staves are free).
`checkResid` nevertheless still degrades at OB. The cause is *not* the config bookkeeping.

### What is in the solution

All 1084 constraints are satisfied by `millepede.res` to ~1e-8 (verified term by term), so there is
no constraint bug. But decomposing the fitted stave-level corrections into a coherent
r-proportional term plus per-half-barrel azimuthal harmonics (ℓ = 0,1,2) leaves almost nothing:

| DOF (stave local frame) | coherent r-term | raw RMS Lr0→Lr6 [µm] | RMS after removing r-term + ℓ≤2 [µm] |
|---|---|---|---|
| TX (rφ, precise)      | −0.46 µm/cm | 8.4 → 79.4  | 0.4 → 36.0 |
| TY (sensor normal ≈ radial) | **−5.44 µm/cm** | 14.1 → 255.0 | 1.4 → 55.2 |
| TZ (beam)             | +0.20 µm/cm | 2.7 → 50.7  | 0.1 → 3.9 |

The local-frame axis assignment is fixed by the code, not guessed: `SensorITS::defineMatrixL2G`
shifts the measurement plane by `(0, 0.5*delta, 0)` to reach the epitaxial layer, and
`getFlatChipAlpha` states that the local X axis lies in the chip plane. So **local Y is the sensor
normal** (≈ radial), X is rφ and Z is beam — the two in-plane, precisely measured directions are
X and Z.

The TY column is the headline: the whole barrel is fitted with a **uniform radial scaling of
ε = −5.4·10⁻⁴** (−0.054%), perfectly linear in radius over all 7 layers — −13 µm at Lr0, −214 µm
at Lr6. Everything else is low-order azimuthal structure per half-barrel. The genuinely local,
stave-by-stave part is 1-2 µm at IB and 20-55 µm at OB (the latter is plausible as *real* OB stave
misalignment; OB mounting tolerances are O(100 µm)).

### Why these modes survive, and why they are the worst possible ones to leave free

1. **The hierarchical constraint is a rigid-body constraint.** Both
   `writeChildrenMeanConstraints` and `pinChildrenMean` impose that the weighted mean of a parent's
   nearest free descendants, *transported into the parent frame*, vanishes — 6 scalar equations per
   parent. That removes the net translation and rotation of each child set and nothing else. The
   stave level spans a whole function space over (r, φ, z); the modes that are *not* net rigid
   motions of any parent — radial scaling, ℓ=2 elliptical squashing, twist — carry no net rigid
   component and pass straight through. Verified numerically: a pure Δr = ε·r mode has only a
   small projection on the 22 half-barrel TY constraints, and the fit simply pays for it with a
   compensating φ-modulation that satisfies them exactly. Note also that a half-barrel is a
   *half*-shell, so an ℓ=0/1 φ-modulation of its staves' radii is a shape change, not a rigid
   motion of the half-barrel, and so is only partly reachable by its own pin. **No choice of
   `pinChildrenMean` can fix this — it is the wrong operator, not a wrong setting.**

2. **Radial scaling is an (almost exactly) null direction of the track χ².** A uniform dilation of
   the transverse plane maps straight lines to straight lines and circles to circles, so it is
   absorbed entirely by rescaling 1/pT and the impact parameter. The only things that break the
   degeneracy are a fixed external length scale — the beam spot, the B field map, and the material
   radii. Here the beam spot is *also* free and unpriored (see point 4), so even that handle is
   gone. Pede therefore takes a large amplitude in this direction for a negligible χ² gain, and
   its fitted value is set by systematics (material budget, cluster-shape bias, B map), not by
   geometry.

3. **It damages `checkResid` maximally at large radius.** A *pure* dilation would be invisible even
   to a leave-one-out estimator; what hurts is the φ-harmonic debris the pin forces in alongside
   it, which is not invisible, plus the coupling of the dilation to the (non-dilated) beam spot and
   momentum scale. This matches the plots: the IB (Lr0-2) panels of `checkResid_hman.pdf` are
   essentially unchanged, while Lr3-6 acquire a 20-30 µm stave-to-stave scatter in ΔY vs φ that was
   not there before, together with a systematic **B=0 vs B=+0.5T split** of ~15 µm. That split is
   the fingerprint of the momentum-scale degeneracy: the fit was given three binaries (B0, B−, B+)
   jointly, and a mode degenerate in one field setting but not another lands on a compromise value
   that suits neither.

4. **The mean vertex escaped the catch-all and is free with no prior.** `{"match": ".*",
   "rigidBody": "fixed"}` sets only a *rigid-body* DOF set; `PVT/meanVertex` carries a
   `MeanVertexDOFSet`, i.e. a *calibration* DOF set, which the clause never touches. Its 30 labels
   (6 time slots × X,Y,Z,SlopeX,SlopeY) are written as `V` with presigma `+0` and make up exactly
   the 30 by which `NVGB` exceeds the 4104 rigid-body DOFs. Fitted: X = +334 µm, Y = −1820 µm,
   Z = −1808 µm (RMS 865 µm over slots). Worth checking against the CCDB `MeanVertex` actually
   loaded — if it was, a 1.8 mm correction is a red flag; if the slots started at (0,0,0) it may
   simply be the real beam position. Either way it should be priored, since it is the last external
   anchor of the global scale.

5. **No diagnostics were available.** `mpSteer.txt` uses `method sparseMINRES 4 0.10`, which
   computes **no parameter errors** — the error column of `millepede.res` is `0.0000` throughout —
   and no eigenvalue spectrum. With `NFGB = 3050` this is an unnecessary choice. Note also that the
   headline "×5 χ²/ndf drop" is iteration 0 → 1 only (1.57e8 → 3.41e7), and the cut factor changes
   from 20 to 10 across exactly that step; iterations 2-4 buy a further 4% while the rejected-record
   count goes 118 → 55728. The ×5 is not a clean measure of geometry improvement.

6. **Weakly-measured DOFs are freed where mechanical knowledge is much better than the fit.** IB
   chip TY (the sensor normal) came out with 15-29 µm RMS and up to 84 µm, over chips sitting on one
   flat silicon ladder whose normal positions are known to a few µm. That DOF can only absorb noise,
   and it is the same channel that produced the `ITSBarFix` ramp (test 1).

### Suggested fixes, in order of value

**(a) Switch the solver and look at the spectrum — do this first.** In `mpSteer.txt` replace
`method sparseMINRES 4 0.10` by

```
  method inversion            5 0.10
```

`NFGB = 3050` inverts in seconds, and this yields real parameter errors and global correlation
coefficients. Run `method diagonalization 5 0.10` once as well: it prints the eigenvalue spectrum
and the eigenvectors, which will show the modes above directly as the smallest eigenvalues. Without
errors there is no way to tell a determined parameter from a weak-mode artefact, which is how a
−214 µm Lr6 shift went unnoticed.

**(b) Constrain the global deformation modes explicitly.** This is the real fix. It is now
**implemented** as the `globalModes` DOF-config clause — see
[Constraining Global Deformation Modes](DOFConfig_rules.md#constraining-global-deformation-modes) and
the subsection *"The `globalModes` clause as implemented"* below. The description of the stop-gap
Python harness is kept for reference.

`macro/add_global_mode_constraints.py` adds the equations
directly to an `mp2con.txt` written by `MilleSteer`, so the hypothesis can be tested *today* with no
rebuild:

```bash
macro/add_global_mode_constraints.py mp2param.txt mp2con.txt -m radial,curl,ellipse
pede mpSteer.txt
```

It appends, over all free stave-level DOFs, `Σ wᵢ rᵢ·TYᵢ = 0` (radial scaling),
`Σ wᵢ rᵢ·TXᵢ = 0` (curl) and `Σ wᵢ cos2φᵢ·TYᵢ = 0`, `Σ wᵢ sin2φᵢ·TYᵢ = 0` (elliptical). On the
`solA1` solution all four are violated at the 1e-3 cm level, i.e. they bite exactly where the
decomposition says they should. It uses nominal layer radii and infers φ from the stave index, which
is why it is a test harness and not the final answer.

The proper implementation is a sibling of `writeRigidBodyConstraints` — say
`writeGlobalModeConstraints` — driven by a new DOF-config clause, e.g.

```json
{ "match": "ITS/.*", "globalModes": ["radialScaling", "curl", "elliptical", "twist", "zScaling"] }
```

For each requested mode it walks the matched subtree and, for every free rigid-body DOF, writes the
projection of the analytic mode displacement field onto that DOF's global displacement direction,
taken from the volume's own `getL2G()` (and weighted by measurement count, as
`writeRigidBodyConstraints` already does). For radial scaling the coefficient of a volume's local-TY
DOF reduces to `r·(n̂·r̂) ≈ r` and of local-TX to `≈ 0`, which is why the hand-written version above
is a faithful stand-in. Modes worth having: radial scaling (Δr ∝ r), z scaling (Δz ∝ z), curl
(Δφ = const), twist/torsion (Δφ ∝ z), elliptical (Δr ∝ cos2φ, sin2φ), bowing (Δr ∝ z²).

#### The `globalModes` clause as implemented

`Volume::writeGlobalModeConstraints` (`src/Volume.cxx`) is the sibling of `writeRigidBodyConstraints`
proposed above; `writeRigidBodyConstraints` calls it once on the root, since a mode equation spans the
whole tree. The clause flags the participating volumes:

```json
{ "match": "ITS/ITSULayer[0-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+",
  "globalModes": ["radialScaling", "curl", "elliptical", "twist"] }
```

Coefficients come from each volume's own `getL2G()`: for every *free* rigid-body DOF, the projection
of the mode's analytic field **u**(x) (translations) or of the local rotation **ω**(x) = ½∇×**u** it
induces (rotations) onto the corresponding local axis. The mode list was later extended — see
[section 5](#5-the-inversion--global-constraint-runs-the-pin-itself-is-the-binding-problem) — and
`doc/ITS_RB_v2.json` now carries the full prescription rather than just `v1` plus this rule.

Verified against the real geometry through a standalone `MilleSteer` run on the `v0` free set:

- 5 equations written, 1084 → 1089 constraints, and the numerical rank of the full set is **1089**
  (the hierarchical 1084 alone rank 1084). The mode equations are therefore linearly independent of
  the mean constraints — direct confirmation that these modes lived in their null space.
- `radialScaling` gets coefficients exactly ∝ r on stave `TY`, constant within a layer
  (0.054 / 0.072 / 0.091 / 0.524 / 0.643 / 0.880 / 1.000 for Lr0→Lr6), *plus* `TX` terms on the IB
  staves only — the real IB stave tilt, which the Python harness could not know about.
- `twist` lands exclusively on the OB staves' `RY` (via **ω**; the IB staves have no free `RY`),
  i.e. the rotation part of the field is what carries this mode, not the translations.
- Evaluated on the `solA1` solution, Σcᵢpᵢ / Σ|cᵢ| is **+174 µm for `radialScaling`** (matching
  ε = −5.4·10⁻⁴ over a ~30 cm mean lever arm), −1.5 µm for `curl` (the Python harness gave
  −1.50 µm — independent cross-check), +50/−20 µm for the two elliptical phases, and −0.18 µm for
  `twist`. So the radial scaling is overwhelmingly the dominant violation, the ℓ=2 terms matter, and
  the twist was never a problem in `solA1`.
- `zScaling` and `bowing` are skipped with a warning at stave level, correctly: both project on a
  volume's local Z with a coefficient ∝ z, and an ITS stave is centred at z = 0. They need free DOFs
  distributed in z (module/chip level) to be observable at all.

**(c) `doc/ITS_RB_v1.json`** — `v0` with the two changes that are justified by mechanical knowledge
rather than by mode-suppression tuning, so the comparison stays attributable:

- `PVT/meanVertex` given an explicit `calib` clause freeing only X, Y, Z (dropping SlopeX/SlopeY)
  plus a 50 µm prior on X/Y and 1 mm on Z. Note the `measurement` must be its own rule, after the
  `calib` rule — a combined rule is fatal.
- **IB chip reduced to `TX, TZ`** (the two in-plane, measured directions), fixing the sensor-normal
  TY. This removes 432 parameters that only absorb noise and closes the test-1 ramp channel at its
  source. Expected `NVGB` ≈ 3690.

Everything else is left identical to `v0` on purpose. `v1` is meant to be run *together with* (b);
on its own it does not remove the radial scaling, because that mode can be realised by any free
radial DOF at stave level.

**(d) Do not expect priors to do (b)'s job.** A `measurement` prior at the honest mechanical
tolerance (O(100 µm) for an OB stave radius) is far too loose to cap a 214 µm mode, and tightening
it to the few-µm level needed would be fabricating knowledge. Priors belong where real external
knowledge exists — the IB half-barrel gap (30 µm, already in `v0`), the beam spot, survey results.

**(e) Optional, once (a)-(c) are in:** the OB stave and its two half-staves both at `"all"` is a
pure reparametrization — 6 stave DOFs + 12 half-stave DOFs reduced by the 6-equation mean
constraint. It is formally well-posed (and the fit does split them sensibly: stave RMS is large,
half-stave RMS is 8-19 µm), but it doubles the OB parameter count for no new freedom. If the
parameter count becomes a problem, keep one level and fix the other rather than both.

### Minor code issues (fixed)

- The warning suppression in `applyMeasurementConfig`, `applyPinChildrenMeanConfig` and the
  `rigidBody` handling now uses `isLiteralPattern()` (`src/Volume.cxx`): an inapplicable rule warns
  only if its pattern has no regex metacharacters, i.e. names one volume explicitly.
- `doc/ITS_RB_from_hstaves_IB.json` and `doc/ITS_RB_from_staves.json` were converted from glob `?`
  to regex `[0-9]`.
- `doc/DOFConfig_rules.md` now documents regex semantics (the README had no glob references).

## 5. The `inversion` + global-constraint runs: the pin itself is the binding problem

Four runs in `/data/devAlg/tests/resOp/` (`solA1Iv0`, `solA1Iv1`, and the same two with the
`macro/add_global_mode_constraints.py` constraints appended, `…GC`), all with
`method inversion 5 0.10`, so with real parameter errors for the first time. All four still degrade
`checkResid`.

### The constraints worked; they were just not the binding constraint

| | `Iv0` | `Iv0GC` |
|---|---|---|
| radial scaling ε | −5.4·10⁻⁴ | **−5.6·10⁻⁸** |
| u_r ℓ=0 at Lr6 | 212 µm | 8 µm |
| u_r ℓ=2 at Lr6 | 125 µm | 5 µm |
| curl θ_z | −71 µrad | −19 µrad |
| pede χ² | 3.2735·10⁷ | 3.2829·10⁷ (**+0.29 %**) |

Removing a 213 µm mode costs 0.29 % of the χ² — the quantitative proof that it was a null direction.
`v1` ≈ `v0` as predicted (Lr6 TY mean −212.8 vs −212.6 µm) at +1.0 % χ² for 444 fewer parameters.

### Root cause: the hierarchical mean constraint is paid off with a lever-arm cancellation

One pin equation, printed term by term — `DOF TY for ITS/ITSULayer6/ITSUHalfBarrel0/ITSUStave0`:

```
TY  ITSUHalfStave0  coef +0.5000   val   -119.39 um    c*p = -59.70 um
RZ  ITSUHalfStave0  coef +0.6450   val  +7523.80 urad  c*p = +48.53 um
TY  ITSUHalfStave1  coef +0.5000   val    -74.72 um    c*p = -37.36 um
RZ  ITSUHalfStave1  coef -0.6450   val  -7523.80 urad  c*p = +48.53 um
                                              sum = 0.0005 um
```

The two half-staves moved −119 and −75 µm (mean **−97 µm**) and the equation is satisfied to
5·10⁻⁴ µm because each took an RZ of **±7.5 mrad**, equal and opposite to 6 digits — the signature of
a null direction. The pin did not prevent the drift; it manufactured a 0.43° rotation to pay for it.

| pin equation | translation coef | rotation coef | cancellation |
|---|---|---|---|
| ITS envelope, TZ | TZ 0.01 | RX 0.28 | RX −8.5 ↔ TZ +8.5 µm |
| ITS envelope, TY | TY 0.01 | RZ 0.28 | RZ −13.5 ↔ TX +11.1, TY +2.5 µm |
| OB half-barrel, TY | TY 0.042 | RZ 1.724 (= r/N, r ≈ 41 cm) | RZ −43.2 ↔ TX +43.0 µm; mean TY only +0.24 µm |
| OB stave, TY | TY 0.50 | RZ 0.645 | RZ +97.1 ↔ TY −97.1 µm |
| IB stave → chip | TX,TY,TZ 0.11 | *none free* | **genuinely zero** |

RZ is the currency — rotation about the sensor normal, ~25 µm of unmeasured normal displacement over
the stave half-width, but weighted by the 41 cm lever arm in the constraint:

```
halfstave  RZ  median  971  p90 5136  max 8058 urad     stave  RZ  median 657  p90 2013  max 5571 urad
halfstave  RX  median  9.7  p90   37  max   89 urad     stave  RX  median  55  p90  128  max  194 urad
halfstave  RY  median  6.1  p90   20  max   46 urad     stave  RY  median  46  p90   87  max  139 urad
```

Wherever the children have free rotations the pin therefore imposes nothing on their mean
displacement. The only level where it bit is IB stave → chip, precisely because the chips have no
free rotations — which is why the chip mode appears as a *ramp* (1st moment) and not an offset.

### What survived, and the evidence that it is the fit and not mechanics

Per harmonic, ε = amplitude/r fitted per layer. A real deformation cannot be coherent in ε across the
IB and the OB — independently mounted structures spanning r = 2.3 → 39 cm. A weak mode must be.

| mode | ε [10⁻⁴] | spread over 7 layers | Lr6 amplitude | was constrained? |
|---|---|---|---|---|
| u_φ ℓ=2 (tangential elliptical) | 1.88–2.52 | **6 %** | **105 µm** | no |
| global rigid tilt θ_y | — | −155 ± 20 µrad on all 14 half-barrels | 69 µm in u_z | no |
| u_r ℓ=3 | 1.24–1.35 | **3 %** | 56 µm | no |
| u_r ℓ=4 | 0.74–1.38 | 14 % | 41 µm | no |
| u_φ ℓ=1 | 0.74–0.81 (OB) | 4 % | 34 µm | no |
| u_r ℓ=0, ℓ=2; u_φ ℓ=0 | 0.05–0.6 | 23–63 % (incoherent) | ≤ 8 µm | **yes ✓** |

Plus two localised modes. **Clamshell / half-shell gap**, the suspected defect but at the OB and
large: relative rigid motion of the two half-barrels ΔT_y = +36, +53, +62, +55 µm at Lr3-6 and
ΔR_x = −200 µrad, constant across all four OB layers (at the IB only +5 to +8 µm). A global tilt and a
clamshell both have zero mean over a half-shell (⟨cos φ⟩ = 0 on [0,π]), so no half-barrel pin can see
either. And the **IB chip dipole**, the `ITSBarFix` ramp back and bigger — chip TY averaged over all
staves of a layer, chip 0 → 8:

```
Lr0  -50.1 -24.4 -12.3  -3.7  +2.6  +8.6 +15.4 +24.0 +39.9   span  76 um
Lr1  -28.9 -17.4  -9.3  -3.1  +2.0  +6.0 +10.9 +16.0 +23.9   span  48 um
Lr2  -18.4  -9.7  -4.1  -0.2  +2.6  +4.1  +6.1  +7.7 +12.0   span  26 um
```

The 0th moment is honestly constrained; the **1st moment is not constrained at all**, and `v0`/`v1`
removed IB stave RX, so the mode has nowhere else to live. In `v1`, with chip TY fixed, it moved into
the **measured** directions as a z² bow: chip TZ +23.0 → −10.6 → +12.5 µm and chip TX +11.3 → −7.9 µm,
against ±2 µm in `v0`. That is why `v1` is not an improvement.

### What was implemented in response

1. **`pinChildrenMeanNoCross`** (`src/Volume.cxx`, `addChildrenMeanTerms(..., noCross)`): the same
   mean equation with the jacobian cross-terms dropped, so only the children's own same-named DOF
   enters and no other DOF of the child can pay for it. `pinChildrenMean` is unchanged and remains the
   right operator wherever the parent's DOF is a fit parameter.
   **Do not request both on the same volume** — over one contributor set with 5-6 free DOFs per child
   they span overlapping subspaces and the constraint set goes rank-deficient (measured: 14 deficient
   equations out of 849 when both were set at the ITS envelope and the OB half-barrels).
   Note what the no-cross form deliberately cannot see: a global rigid shift or tilt is a cos φ / sin φ
   pattern in the children's local DOFs, so its plain mean vanishes. That is what the cross-terms are
   for, and what `globalTilt` is for.
2. **Two generic harmonic families and `globalTilt`** added to `Volume::GlobalMode`:
   `radialHarmonic` (Δr = ε·r·cos/sin(ℓφ)) and `tangentialHarmonic` (rΔφ = ε·r·cos/sin(ℓφ)),
   ℓ = 1…4, and the two transverse partners of `curl`. The amplitude carries an explicit factor r
   because that is the measured shape (the ε table above). At ℓ = 0 the families reduce exactly to
   `radialScaling` and `curl`. The mode set is now indexed by a *flat equation index*
   (`Volume::NGlobalEq` = 25), which the `globalModes` clause flags individually, so a name may select
   a whole mode, one harmonic order (`radialHarmonic3`) or one equation (`globalTilt_y`).
3. **Propagated** to `src/WeakModes.cxx` (one diagnostic column per equation, resolved through
   `Volume::globalEqName`) and to `macro/drawWeakModes.C`, which now resolves the `uMode` columns by
   name from the file instead of by hard-coded index and has a new page for the harmonic families.
   `macro/add_global_mode_constraints.py` gained the same families and a note that it is superseded.
4. **`doc/ITS_RB_v2.json`**: IB chips fixed outright and IB stave RX/RY restored (the stave tilt is the
   physical parameter; the chips sit on one silicon ladder flat to a few µm); stave and half-stave RZ
   fixed (it buys ~25 µm of unmeasured displacement while carrying a 41 cm lever arm into every
   constraint); `pinChildrenMeanNoCross` at the ITS envelope and the OB half-barrels; the global-mode
   set above at the stave level. Verified: NVGB = 2427, 803 constraints of which 21 global-mode, and
   the constraint matrix is **full rank 803**, the global modes contributing 21 independent equations.

Validation of the new equations against the `solA1Iv0GC` solution, Σcᵢpᵢ/Σ|cᵢ|:

```
radialScaling                 +0.06 um     (already constrained there - cross-check of the new code)
radialHarmonic2_cos/sin2phi   -1.21/+0.59  (ditto, the harness's radial ellipse)
tangentialHarmonic2_cos/sin   -44.0/+49.4  <== the dominant survivor
globalTilt_y                  -41.1        <== the 155 urad tilt   (globalTilt_x only +5.4)
radialHarmonic3_sin3phi       +33.6        <== u_r l=3
radialHarmonic4_cos/sin4phi   -18.5/+18.6  <== u_r l=4
tangentialHarmonic1_cos/sin   +18.3/+11.2  <== u_phi l=1
twist                         -0.15        (never was a problem)
```

Every equation bites exactly where the decomposition said it should, and the two the harness had
already constrained come out at ~0 — an independent check of the geometry-based coefficients against
the index-based ones.

### Still open

- The **OB clamshell** (ΔT_y ≈ 55 µm, ΔR_x ≈ −200 µrad between the two half-shells) is not addressed
  by any of the above: it is a relative rigid motion of two half-shells, invisible to a mean over
  either of them and not a global deformation mode. It needs either a `measurement` prior on the OB
  half-barrel poses (freeing them first, as the IB ones are) or a dedicated clamshell constraint.
- `zScaling` and `bowing` remain unobservable at stave level (coefficient ∝ z of the volume, and a
  stave is centred at z = 0). They need the module/chip level free to be constrained at all, and with
  the IB chips now fixed in `v2` there is no level left that is distributed in z.
- The **IB chip 1st moment** is now moot in `v2` (chips fixed), but if chip DOFs are ever freed again
  the dipole needs its own constraint, not just the mean.

## 6. `solA1Iv2GC`: what checkResid still shows, and `ITS_RB_v3.json`

### The v2 residuals are better than V10n; the remaining structure is below the stave

- **Binning:** the apparent growth of the OB dy-vs-φ scatter in the first v2 PDF was a binning
  mismatch: the V10n hman files had 30 φ bins, the MP2 ones 100. With equal binning (RMS over φ bins
  of the truncated means, pp B0):
  - dy, Lr3-6: 11.8/11.1/28.1/7.8 → 8.5/8.5/19.7/6.5 µm;
  - dz, Lr3-6: 11.7/9.7/23.7/13.9 → 3.8/2.9/8.9/2.8 µm.
  OO B+ behaves the same. v2 ≈ v1GC, as expected: the modes it constrains are invisible to the
  leave-one-out residual.
- **Per sensor, from the trees:**
  - Per OB half-stave the v2 correction is anticorrelated with the V10n residual (−0.64…−0.88), and
    the half-stave mean RMS drops to 4-8 µm (dy) and 2-4 µm (dz).
  - Inside a half-stave the per-chip dy RMS is still 12-24 µm, dominated by an offset between the
    two chip rows of a module (rows 7-13 vs 0-6: mean −10…−22 µm, RMS 10-29 µm). It is reproduced
    in pp B0 and OO B+ (per-chip correlation 0.73-0.84), so it is real, and no v2 DOF can touch it.
- **Large-scale Z/tgλ structures (p.3/4/7/8):** these are identical for V10n and v2 and in both
  datasets.
  - In the IB, a per-sensor fit of dy, dz on {1, z_loc, tgλ, z_loc·tgλ} shows:
    - an in-plane chip rotation coherent over all staves and growing with z (Lr0: −3.5 … +2.5 µm/cm,
      chips 0 → 8), which is the 3 cm sawtooth in dy vs Z;
    - end-chip offsets of a few µm.
  - In the OB the structures sit at the module scale, on the outer modules.

  All need DOFs distributed in z.
- **Lr5 dz:** it is anomalously broad (RMS 350 µm against 220-245 µm in the other OB layers);
  still to be understood.

### The modes that DOFs distributed in z open: the polar families

For a track from the vertex, z ≈ r·t (t = tgλ), so fields r·f(t) (rφ, z, radial) are absorbed by
the track φ, tgλ and curvature, and fields f(t) by d0/z0 (broken only by the vertex). The existing
modes only cover the t-independent members, which sufficed while every free DOF sat at z ≈ 0.
`radialPolar`, `tangentialPolar` and `longitudinalPolar` (`Volume::GlobalMode`, 14 equations, see
[DOFConfig_rules.md](DOFConfig_rules.md#constraining-global-deformation-modes)) add the members
with p ∈ {0,1}, k ≤ 2. The IB chip-TY dipole of section 5 (span 76/48/26 µm for Lr0/1/2, i.e.
∝ 1/r) is `radialPolar_r0t1`, Δr ∝ z/r. `bowing` (Δr = ε·z²) is not a null direction; its null
counterpart is `radialPolar_r0t2`.

### `doc/ITS_RB_v3.json`

v2 plus DOFs distributed in z, with one representation per deformation:

- **IB chips:** `TX, TZ, RY` (in-plane shift and rotation, measured directly), plus `TY` with a
  5 µm prior. Normal positions are measured only at large |z|, through dz = Δn·tgλ.
- **OB modules:** the same, with a 30 µm prior on `TY`.
- **Tilt `RX` dropped** at the IB stave, the OB stave and the OB half-stave, and from the OB
  half-barrel pin: a tilt is now a `TY` ramp of the children, and keeping both is the ITSBarFix
  degeneracy.
- **`globalModes`** over all of `ITS/.+`: the v2 set plus `zScaling` and the three polar families.

MilleSteer on the all-free data:
- NVGB = 10458;
- 2010 constraints (1974 hierarchical + 36 global modes), numerically full rank, smallest
  normalised singular value 0.024;
- 2160 priors.

OB chip `TX` (the chip-row offsets) is the next step; it is left out here to keep the matrix at a
size `inversion` handles.

### pede crash with long constraints (Millepede-II 4fd638f67f)

pede aborted on the first v3 steering set with `free(): invalid next size` in LOOP1. A
bounds-checked build of the same commit locates it in `grpcon` (`pede.f90:1902`, `globalParCons(0)`).

**Why:** constraints are grouped in the order of their lowest parameter index. A constraint not yet
grouped joins the first group found through a shared parameter, and relabels only its *direct*
neighbours. When a constraint joins several groups formed before it, an earlier group can end up
empty, and the redundancy check then indexes (and in the optimised build writes) out of bounds.

**Seen in both runs:** a replica of the algorithm (pede's group count reproduced exactly) gives
- v2: 436 groups, one of them empty (it survived by luck);
- v3: 3 groups, two empty.

**Workaround:** `Volume::writeGlobalModeConstraints` gives every mode equation the lowest label of
any free rigid-body DOF of the tree, as an anchor term with coefficient 1e-12. The coefficient cannot
be 0: pede's text reader (`addItem`) silently skips zero factors. pede then processes the mode
equations first and nothing is merged afterwards (v3: 1 group, none empty). Against a largest
coefficient of 1, the anchor term is negligible for the equation and its rank.

### The v3 fit (`/data0/devAlg/tests/solA1Iv3`): the polar families are not enough

**Fit statistics:**
- NVGB = 10454, NFGB = 8444; constraint rank 2010/2010.
- Objective at iteration 3 (same cuts in both runs): 3.328·10⁷ against 3.424·10⁷ for v2 (−2.8 %).
- Final χ²/ndf: 1.025 against 1.12.

**The corrections are much larger than v2's.** Sensor displacement RMS per layer:

| | v2 | v3 |
|---|---|---|
| rφ, IB | 6-7 µm | 12-14 µm |
| rφ, OB | 20-45 µm | 33-85 µm |
| z, OB | 7-14 µm | 30-66 µm |

Lr5 module `TX` goes up to 864 µm. Per-chip residuals in checkResid are 2-4 µm (IB), so most of
this is in directions the tracks do not see.

**It is a ray mode.** Sensors are binned in (φ, t = z/r) cells, and the displacement divided by r is
compared between layers:
- **dz/r:** 87 % of the per-cell variance is common to all seven layers (Lr0-Lr6 correlation 0.85;
  70 % already in v2, at about a third of the amplitude). That is Δz = r·g(φ, t), a per-direction
  change of tgλ, invisible to any track from the vertex.
- **rφ/r:** common within the IB (correlation 0.99) and within the OB (0.97).

The polar families constrain only the axisymmetric low orders of g (t^k, k ≤ 2), and `globalTilt`
only cos φ, sin φ at t⁰. With DOFs distributed in z, the null space is effectively a function of
(φ, t) per displacement component, not a short list of global modes.

**Options, not yet tried:**
- **Ray constraints:** Σ_layers w·Δz/r = 0 and Σ_layers w·Δ(rφ)/r = 0 per (φ, t) cell, i.e. the
  general form of the polar/harmonic expansion.
- **Priors:** `measurement` priors on module/chip `TX`, `TZ` at the mechanical placement tolerance.
- **Data that breaks the ray degeneracy:** cosmics, tracks not from the vertex. The vertex z spread
  already breaks it partly at the IB, much less at the OB.

checkResid cannot arbitrate between these: by construction it is blind to the ray mode.

### `ITS_RB_v4.json`: ray constraints + priors (`/data0/devAlg/tests/solA1Iv4`)

What changed against v3:
- **`rayModes` on `ITS`:** 24 φ × 16 t cells over |t| ≤ 4, all three components. That gives 720
  equations; 120 cells hold a single layer only and are skipped. See
  [DOFConfig_rules.md](DOFConfig_rules.md#constraining-ray-modes).
- **Global modes reduced** to the ones outside the ray span: `twist` and the p = 0 polar members.
- **Priors on the in-plane DOFs:** IB chips TX, TZ 10 µm (TY 5 µm), OB modules TX, TY, TZ 30 µm.
  These are assumed tolerances, to be replaced by survey numbers if available.
- **IB stave `RY` dropped.** With it free, the fit realised a common in-plane rotation of a stave's
  chips about fixed centres as stave RY plus a compensating chip-TX ramp: the mean constraint forbids
  the same thing through chip RY. The result was a chip-TX RMS of 32 µm in violation of its 10 µm
  prior (chip TX vs RY·z correlation 0.88). Without stave RY the chip-TX RMS is 15 µm and chip RY
  carries the rotation; the χ² is slightly lower with 48 fewer parameters. The OB shows no such
  cancellation (module TX vs RY·z correlation −0.23).

Constraint set: 2655 constraints, rank 2655, one pede group; NFGB = 7751.

| | v2 | v3 | v4 |
|---|---|---|---|
| objective at iteration 3 (same cuts) | 3.424·10⁷ | 3.328·10⁷ | 3.348·10⁷ |
| final χ²/ndf | 1.12 | 1.025 | 1.030 |
| sensor displacement RMS, all ITS | 59 µm | 104 µm | 42 µm |
| rφ RMS, IB / OB | 6-7 / 20-45 µm | 12-14 / 33-85 µm | 7-8 / 13-20 µm |
| dz RMS, IB / OB | 1 / 7-14 µm | 7-9 / 30-66 µm | 7-10 / 17-24 µm |
| ray-common dz/r: fraction, RMS | 0.70, 0.34·10⁻⁴ | 0.87, 1.01·10⁻⁴ | < 0, 0.19·10⁻⁴ |
| OB module TX RMS (max) | fixed | 73 (864) µm | 22 (91) µm |

The ray constraints and priors cost 0.6 % of the χ² that v3 gained over v2, against a factor 2.5
smaller displacements: what they removed was nearly invisible to the tracks.

**What is left:**
- An anti-correlated Lr3 vs Lr6 pattern in rφ/r and dr (correlation −0.97, −0.80) at a small
  amplitude (common part 0.7·10⁻⁴ and 2 µm). It has zero mean per cell, so it passes the ray
  constraint: a sagitta-like (∝ r²) pattern per direction, absorbed by the curvature when B ≠ 0.
  Like the twist, the B = 0 data should break it; worth watching.
- The v4 geometry has not yet been run through checkResid.

### checkResid on v4: the IB gets worse, from IB chip DOFs the residuals reject

The OB improves, but several IB and vertex panels get worse than v2, the same in pp B0 and OO B+.
RMS over bins of the truncated means, pp, in µm:

| | v2 | v4 |
|---|---|---|
| Lr0 dy vs Z | 1.7 | 4.6 |
| Lr2 dy vs Z | 1.2 | 2.8 |
| Vtx dz vs φ | 1.8 | 5.4 |
| Vtx dz vs tgλ | 3.0 | 7.8 |
| IB dz vs φ, Lr0/1/2 | 1.4/0.8/1.2 | 2.1/1.4/2.0 |
| Lr6 dz vs Z | 8.7 | 15.8 |

What is established:
- **The Lr0 dy-vs-Z sawtooth is the new IB chip `RY`** (in-plane rotation). The fit put a common
  rotation per layer on the chips: +4.97, +1.15, −2.26·10⁻⁴ rad at Lr0/1/2, with errors of about
  0.4·10⁻⁴. The V10n in-chip dy slope that checkResid measures for those layers is +0.09, +0.02,
  −0.53 µm/cm, i.e. the reconstruction sees no such rotation. The chip-to-chip part of RY does
  anticorrelate with the V10n slopes (−0.46/−0.48 at Lr1/2).
- **The geometry is built as fitted.** `algdump.root` (pos0 = ideal) shows v4 − v2 in-plane chip
  rotations of +4.94/+1.13/−2.20·10⁻⁴ rad. Chip centres hardly move with z (≤ 0.18 µm/cm against
  5 µm/cm), so nothing compensates the rotation.
- **Not a convention error.** `AlignParam` rotations are right-handed and active, as in
  `RigidBodyDOFSet::fillDerivatives`; `computeJacobianL2T` moves the pivot correctly.
- **Not the IB part of the ray cells.** That part is not a null direction (σ(z_v) = 4.9/5.4 cm in
  pp/OO). But v5 (`ITS_RB_v5.json` = v4 with `rMin` 10 cm) gives the same χ² (3.3479 against
  3.3478·10⁷ at iteration 3) and the same IB chip parameters (Lr0 RY mean +4.85, changes ≤ 4 µm).

**Open:** the GBL fit on the Mille data and the KF of checkResid disagree at IB chip level. The fit
sees an in-chip slope in V10n that the reconstruction does not; the same chip-level DOFs also
failed in v1 (IB chip TX/TZ). Candidates:
- the workflow's IB measurement frame: `SensorITS::defineMatrixL2G` (chip node + epitaxial shift)
  and `getFlatChipAlpha` against the reconstruction's sensor matrix and tracking frame, together
  with how the cluster coordinates are carried into it;
- a difference between the Mille runs and the checkResid runs.

The discriminating test is a per-point residual dump of the alignment workflow on V10n (the debug
stream currently writes only `MisRes`), compared chip by chip with checkResid's in-chip dy-vs-z
slopes.

### `WeakModesDiag` fit B was void (also in v2)

`fitModes` used `Eigen::BDCSVD` on a design matrix that is exactly rank deficient by construction
(`curl` and `globalTilt` duplicate the rigid rotation columns). It returned NaN singular values:
- v2: rank 12, `conditioning` null, and the RMS unchanged after the rigid part;
- v3: rank 0.

With `JacobiSVD` the v3 fit gives rank 42/45, the three missing being exactly those duplicates.
Metric B values in earlier `weakModes.json` files should be disregarded; metric A was unaffected.

## Caveat / next step

The half-stave and chip DOF sets determine which derivative columns exist in the Mille
binary, so this requires **regenerating MilleData (stage 1)**, not just a new `pede`
run on the existing `mp2tree.txt` — in particular the one-digit-OB-stave half-stave
fix has no effect without it. After regenerating, re-run `checkResid` and check (a)
stave-level mean/RMS per layer stay flat and small (as in the pinned run), and (b) the
fine-grained IB dy/dz-vs-Z plots no longer show the chip-periodic ramp, before judging
the result against the initial (pre-alignment) geometry.
