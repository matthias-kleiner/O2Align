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

## Caveat / next step

The half-stave and chip DOF sets determine which derivative columns exist in the Mille
binary, so this requires **regenerating MilleData (stage 1)**, not just a new `pede`
run on the existing `mp2tree.txt` — in particular the one-digit-OB-stave half-stave
fix has no effect without it. After regenerating, re-run `checkResid` and check (a)
stave-level mean/RMS per layer stay flat and small (as in the pinned run), and (b) the
fine-grained IB dy/dz-vs-Z plots no longer show the chip-periodic ramp, before judging
the result against the initial (pre-alignment) geometry.
