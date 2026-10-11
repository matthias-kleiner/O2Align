# DOF Configuration Rules

`Volume::applyDOFConfig` traverses the hierarchy and applies matching rules to non-pseudo volumes. A `defaults` object, if present, is applied first to every volume; rules then run in order, and a later matching clause replaces the earlier DOF set of the same kind. Patterns are matched against symbolic volume names, as full-match ECMAScript regular expressions (not globs: `.` is any character, `?` and `*` are quantifiers), also retried with an implicit leading `.*`. Since `.*` can match `/`, use exact-depth patterns when targeting a hierarchy level.

## Rigid-Body DOFs

The `rigidBody` clause frees an arbitrary subset of a volume's six DOFs (`TX,TY,TZ,RX,RY,RZ`, in
the volume's local frame), all the others of that volume staying fixed; a volume not matched by any
rule keeps none free. Free several at once by listing them together, e.g. the two in-plane shifts
plus the in-plane rotation of every ITS stave:

```json
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]", "rigidBody": ["TX", "TZ", "RZ"] }
```

`"rigidBody": "all"` frees all six, `"rigidBody": "fixed"` (or `false`) fixes all six — the latter
is only needed to override a `defaults` clause or an earlier rule, since no DOF is free by default.
As with every other DOF-set clause, a later matching rule replaces the DOF set of an earlier one
rather than merging with it, so listing a subset is not cumulative across rules.

## Mean Vertex

The mean vertex volume is `PVT/meanVertex`. Its calibration DOFs are the position `X, Y, Z` and the
beam-line slopes `SlopeX, SlopeY` (dX/dZ, dY/dZ), all 5 defaulting to free. Each of the 5 is
independent of the others — any subset may be free while the rest stay fixed, with no all-or-nothing
grouping. A fixed DOF still enters the fit at its current prior value (from the `MeanVertexObject`
delivered by the CCDB for the active calibration slot), it is simply not adjusted.

Enable all 5 with:

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex" } }
```

Fix all 5 with:

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex", "fixed": true } }
```

Fix only the slopes, keeping X, Y, Z free:

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex", "fix": ["SlopeX", "SlopeY"] } }
```

As with any calibration DOF set, an external measurement of an individual DOF (e.g. a value derived
from a prior calibration) can be added via the `measurement` clause, see
[Measurements of Rigid-Body and Calibration DOFs](#measurements-of-rigid-body-and-calibration-dofs)
below; it goes in a separate rule following the one that frees/fixes the DOFs.

## Automatic vs. Optional Constraints at a Parent/Children Boundary

`writeRigidBodyConstraints` decides, for every rigid-body DOF of every volume, whether a constraint
over its children is needed at all, and whether it is written automatically or only on request. The
two cases that matter in practice:

1. **Parent DOF free, and some descendant has *any* rigid-body DOF free** — not necessarily the same
   one (a chain of free levels, with possibly several non-rigid-body/fixed levels looked through in
   between — see `writeChildrenMeanConstraints`). `addChildrenMeanTerms` builds the parent's DOF `iDOF`
   constraint from *every* free DOF `jDOF` of each contributing child, weighted by the full jacobian
   element `j(iDOF, jDOF)`, not only `j(iDOF, iDOF)`: a child's `RX`, say, contributes to the parent's
   `TY` constraint too, with a coefficient set by the lever arm between them. This is not an
   approximation — `getJP2L()`/`getJL2P()` is the same fixed 6x6 matrix per child used by
   `buildPointGlobals` to chain every leaf measurement up through the hierarchy, so a parent's free DOF
   and a differently-named free DOF of a child really can be exactly indistinguishable once propagated
   through it: a genuine degeneracy of the fit. The constraint "weighted mean of the children's
   movement (transported to the parent's frame) vanishes" is written **automatically**, with no DOF
   config needed beyond freeing both levels:
   ```json
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]", "rigidBody": ["TY"] },
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]", "rigidBody": ["TY"] }
   ```
   Here the half-barrel's `TY` and the mean `TY` of its staves are tied together automatically; a stave
   free only in, say, `RZ` instead would still contribute to that same constraint via `j(TY, RZ)`.

2. **Parent DOF fixed, descendants' DOF free.** There is no parameter at the parent level to be
   degenerate with, so **no constraint is written by default**: the children (and the hierarchy below
   them) are free to drift together in that DOF, which may be a genuinely unconstrained or only
   weakly constrained mode of the fit (e.g. nothing anchors the overall `TY` of a detector whose
   envelope is fixed but whose staves are all free in `TY`). Use `pinChildrenMean` to request the same
   mean-zero constraint explicitly even though the parent DOF is fixed — this **does not** turn the
   parent DOF into a fit parameter, it only pins the common mode of what's below it to the parent's
   (fixed) position:
   ```json
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]", "rigidBody": "fixed", "pinChildrenMean": ["TY"] },
   { "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]", "rigidBody": ["TY"] }
   ```
   Without the `pinChildrenMean` clause, this configuration leaves the mean `TY` of the staves
   unconstrained; with it, the mean is pinned to the (fixed) half-barrel position.

A third, unproblematic case is a free parent DOF with **no** free descendant DOF below it at all (see
[Aligning a Subtree as One Rigid Body](#aligning-a-subtree-as-one-rigid-body)): there is nothing to be
degenerate with, so no constraint is needed or written, automatically or otherwise, and the parent
DOF is still fully observable from the leaf measurements.

None of this applies to calibration DOFs (e.g. the mean vertex, TPC drift): they are never organized
in a parent/children hierarchy of their own, so there is no equivalent automatic constraint — only the
explicit `measurement` clause (previous section) adds a prior to a calibration DOF.

## Pinning the Common Mode of Children

`Volume::writeRigidBodyConstraints` forces, for every *free* rigid-body DOF of a parent, the
weighted mean of its active children's movement, transported to the parent's frame and projected
onto that DOF, to vanish — summed over *every* free DOF of each child (not only the same-named one),
weighted by the corresponding element of the child-to-parent jacobian; see
[Automatic vs. Optional Constraints](#automatic-vs-optional-constraints-at-a-parentchildren-boundary)
for why a child's differently-named free DOF is included on exactly the same footing. This breaks a
real degeneracy: moving the parent and moving all its children together in the same direction are
otherwise indistinguishable. A child without any free rigid-body DOF (no `rigidBody` rule, or a
`fixed` one) is looked through: its own children, with the jacobians chained, take its place, and so
on down to the nearest level having a free DOF. The walk stops at the first descendant with *any*
free DOF, so free DOFs of a grandchild along a direction its (partially free) parent keeps fixed
are not included in the mean — avoid such configurations. When the parent DOF is fixed instead, that degeneracy doesn't arise,
so no such constraint is written by default — the children are then free to collectively drift in
that DOF, which may be a genuinely unconstrained (or only weakly constrained) mode of the fit.

The `pinChildrenMean` clause asks for the same mean-zero constraint on a DOF of the parent even
though that DOF is fixed there, to remove such a mode explicitly. It takes a boolean (all DOFs), the
string `"all"`, or an array of DOF names (`TX,TY,TZ,RX,RY,RZ`), and — like `rigidBody`/`calib` — a
later matching rule replaces the flags of an earlier one rather than merging with them:

```json
{ "match": "TRD/sm[0-9][0-9]", "rigidBody": "fixed", "pinChildrenMean": ["TY", "RZ"] }
```

With the chambers free, this pins the mean `TY` shift and `RZ` rotation of the chambers of every
supermodule to the (fixed) supermodule position, without making the supermodule itself a fit
parameter. The flags don't depend on the volume's `rigidBody` clause: they also work on a volume with
no rigid-body DOF set at all. A pinned DOF that no free descendant DOF contributes to (nothing free below,
or none of the free DOFs feeds it through the chained jacobians) is skipped with a warning.
`pinChildrenMean` is ignored, with a warning, on a volume that is not rigid-body alignable (the root,
or an envelope without its own geometry such as `TRD_envelope`), so the common mode of the
supermodules of a detector can't be pinned this way.

### `pinChildrenMeanNoCross`: the same pin without the jacobian cross-terms

The transported form above sums *every* free DOF of each child weighted by a jacobian element, and
those elements carry the lever arm from the child to the parent. That is the correct operator when
the parent's DOF is a fit parameter — the degeneracy genuinely is "parent pose vs transported
children pose". It is the **wrong** operator for pinning a collective drift, because the fit can
satisfy it with whichever contributing DOF is cheapest in χ², and a weakly-observable rotation with a
large lever arm is by far the cheapest. Measured in the ITS (see
[AlignmentWeakModes_handover.md](AlignmentWeakModes_handover.md)): in `DOF TY for …/ITSUStave0` the two
half-staves moved −119 and −75 µm, and the equation was satisfied to 5·10⁻⁴ µm because each took an
RZ of ±7.5 mrad — rotation about the sensor normal, which displaces the measured coordinates by
almost nothing but enters the equation with a coefficient of 0.645 against 0.5 for TY. The pin did not
prevent the 97 µm drift; it manufactured a 0.43° rotation to pay for it. At the OB half-barrel the
coefficient ratio is 1.724 vs 0.042 — a factor of 41, the stave radius in cm.

`pinChildrenMeanNoCross` writes the same equation with the cross-terms dropped: the plain mean of the
children's *own* same-named DOF, which no other DOF of the child can pay for.

```json
{ "match": "ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-1]", "rigidBody": "fixed",
  "pinChildrenMeanNoCross": ["TX", "TY", "TZ", "RX", "RY"] }
```

It takes the same values as `pinChildrenMean` (a bool, `"all"`, or a list of DOF names), has the same
"replace, don't merge" semantics, and its flags are independent — so both clauses may appear on one
volume and each writes its own equation. **In practice, use one or the other per volume**: over the
same set of contributors the two span largely overlapping subspaces, and requesting both at a level
with 5-6 free DOFs per child makes the constraint set rank-deficient (pede then reports a rank below
the number of constraint equations). Keep the transported form where the parent's DOF is free, and
use the cross-term-free form to pin a fixed parent's children.

Two things it deliberately cannot do, because they are not collective *drifts* of any single DOF:

- a global rigid **shift** or **tilt** of the whole branch is a cos φ / sin φ pattern in each child's
  local DOFs, so its plain mean vanishes and this clause is blind to it (the transported form does see
  it — that is what the cross-terms are for);
- any deformation mode — use [`globalModes`](#constraining-global-deformation-modes).

A pinned DOF that no descendant has free is skipped with a warning, which is the expected message
when, say, the stave RZ below an OB half-barrel is fixed.

## Constraining Global Deformation Modes

`pinChildrenMean` and the automatic mean constraint are both *rigid-body* constraints: 6 scalar
equations per parent, removing the net translation and rotation of that parent's set of free
descendants and nothing else. A coherent deformation of the barrel as a whole — a uniform radial
scaling, an elliptical squashing of each layer, a twist — is not a net rigid motion of any single
parent, carries no net rigid component, and therefore passes straight through it, however
`pinChildrenMean` is set. Such modes are also (very nearly) null directions of the track χ², so
Millepede takes large amplitudes in them for a negligible χ² gain. This is a different operator, not
a different setting of the same one.

The `globalModes` clause flags the volumes whose free rigid-body DOFs take part in one global
equation per mode, written by `Volume::writeGlobalModeConstraints` (called once on the root from
`writeRigidBodyConstraints`), forcing the total amplitude of that mode over the whole flagged set to
vanish:

```json
{ "match": "ITS/ITSULayer[0-6]/ITSUHalfBarrel[0-1]/ITSUStave[0-9]+",
  "globalModes": ["radialScaling", "elliptical", "twist"] }
```

| mode | displacement field | equations |
|---|---|---|
| `radialScaling` | Δr = ε·r (uniform transverse dilation, absorbed by the 1/p<sub>T</sub> scale) | 1 |
| `zScaling`      | Δz = ε·z | 1 |
| `curl`          | Δφ = ε (rigid rotation about the beam) | 1 |
| `twist`         | Δφ = ε·z (torsion along the beam) | 1 |
| `elliptical`    | Δr = ε·cos2φ and Δr = ε·sin2φ (r-independent amplitude) | **2** |
| `bowing`        | Δr = ε·z² (radial sagitta along the beam) | 1 |
| `radialHarmonic` | Δr = ε·r·cos(ℓφ) and ε·r·sin(ℓφ), ℓ = 1…`GlobalHarmonicLMax` (4) | **8** |
| `tangentialHarmonic` | rΔφ = ε·r·cos(ℓφ) and ε·r·sin(ℓφ), same ℓ range | **8** |
| `globalTilt`    | rigid rotation of the branch about the global X and about the global Y axis | **2** |

The two harmonic families are the generic form of the azimuthal modes, with the amplitude explicitly
proportional to r because that is the measured shape: the fitted ε = amplitude/r of the surviving ITS
modes is coherent to 2-6 % across the IB and the OB, two independently mounted structures spanning
r = 2.3 → 39 cm. At ℓ = 0 they reduce exactly to `radialScaling` and `curl`, which is why their orders
start at ℓ = 1, and `tangentialHarmonic2` is the companion that `elliptical` — a purely radial field —
does not cover. `globalTilt` supplies the two transverse partners of `curl`; it is needed as an
explicit constraint because a tilt shows up as a child TZ pattern ∝ −r·cos φ whose mean vanishes both
over a half-shell and over the full circle, so no `pinChildrenMean*` of the TZ DOF can see it.

The clause takes a boolean (all equations), the string `"all"`, or an array of names, each of which may
be a mode name (all its equations), a single order of a harmonic family (`radialHarmonic3`), or a
single equation (`elliptical_cos2phi`, `tangentialHarmonic2_sin2phi`, `globalTilt_y`). Like
`rigidBody`/`calib`/`pinChildrenMean` a later matching rule replaces the flags of an earlier one.
It creates no DOF set and frees nothing: it only adds terms for DOFs that are already free, so it
must follow the rule that frees them. Flagging several nested levels at once is legitimate and is not
double counting: the equation then constrains the total amplitude summed over all of them, which is
the degenerate quantity, regardless of which level realises it.

The coefficients come from each volume's own `getL2G()`, not from nominal radii: for every free DOF,
the projection of the mode's analytic displacement field **u**(x) (for `TX,TY,TZ`) or of the local
rotation **ω**(x) = ½∇×**u** it induces (for `RX,RY,RZ`) onto the corresponding local axis. The
rotation part is not optional: it is what makes `twist` constrain the OB staves' `RY` rather than
their translations. Each equation is homogeneous, so it is normalised to a largest coefficient of 1,
and terms below 1e-9 of that are dropped as rounding residues of the L2G matrices.

A mode with no handle on the flagged level is skipped with a warning rather than silently. `zScaling`
and `bowing` project on a volume's local Z (beam) DOF with a coefficient ∝ z of that volume, which
vanishes for an ITS stave — a stave is centred at z = 0 — so these two need free DOFs that are
*distributed* in z (module or chip level) to be observable at all. `globalModes` is ignored, with a
warning, on a volume that is not rigid-body alignable.

`globalTilt` is nearly the same equation as a `pinChildrenMean` (transported form) of the parent's
`RX`/`RY`, so requesting both makes the constraint set rank-deficient. Prefer `globalTilt`: the pin's
contributor walk stops at the first descendant with any free DOF, which in the ITS excludes the 48 IB
staves, while `globalTilt` flagged at the stave level covers all 192.

Adding a mode to `Volume::GlobalMode` propagates automatically to the `WeakModesDiag` output and to
`macro/drawWeakModes.C`, which resolve the mode columns by name; a `static_assert` after the class
guards `Volume::NGlobalEq` against going out of step.

## Measurements of Rigid-Body and Calibration DOFs

The `measurement` clause supplies an external measurement (e.g. from a survey, or a value carried
over from a prior calibration) of individual DOFs of the matching volumes — rigid-body DOFs,
calibration DOFs, or both at once (each name is resolved against whichever of the volume's two DOF
sets has it; the two never share a name, e.g. `TX..RZ` vs. `X,Y,Z,SlopeX,SlopeY` for the mean
vertex). It is written to `mp2con.txt` as a Millepede `Measurement` record (`Measurement value
sigma`, followed by the label of the DOF with coefficient 1), i.e. a soft constraint `p = value ±
sigma` which, unlike a `Constraint`, does not fix the parameter. It is an object `DOF name ->
[value, sigma]` (or `{"value": v, "sigma": s}`), in the units of the DOF (cm for `TX,TY,TZ` and the
mean-vertex `X,Y,Z`; rad for `RX,RY,RZ`; dimensionless for `SlopeX,SlopeY`), relative to the geometry
used by the fit:

```json
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]", "rigidBody": ["TX", "TY"] },
{ "match": "ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]", "measurement": { "TX": [0.01, 0.005], "TY": {"value": -0.02, "sigma": 0.01} } }
```

```json
{ "match": "PVT/meanVertex", "calib": { "type": "meanvertex" } },
{ "match": "PVT/meanVertex", "measurement": { "SlopeY": [0.0002, 0.0001] } }
```

A `rigidBody`/`calib` clause *defines* a DOF set (listed DOFs free, all the others fixed), while
`measurement` only adds a prior to DOFs that are already free. A rule containing `measurement`
together with `rigidBody` or `calib` is therefore rejected (fatal): put the measurement in its own
rule, after the one defining the DOFs.

As for `pinChildrenMean`, no DOF set is created and a later matching rule replaces all the
measurements of both DOF sets of the volume (even if the new rule only touches one of them).
`sigma` must be positive. A measurement of a DOF which is not free at the time of writing, as well as
the clause on a volume having neither a rigid-body nor a calibration DOF set at all, is ignored with
a warning (under a literal pattern, without regex metacharacters) or silently skipped (under any other pattern, e.g.
a `measurement` rule meant for a different kind of volume).

## Aligning a Subtree as One Rigid Body

The opposite situation to `pinChildrenMean` is freeing a parent's rigid-body DOF while giving none
of its descendants a `rigidBody` rule at all. This is not auto-disabled and needs no mean-of-children
constraint: the parent-to-child jacobian is applied unconditionally down to every leaf measurement
regardless of whether a descendant carries a DOFSet, so the parent's free DOF is fully observable on
its own. A constraint is only needed to resolve the degeneracy of a DOF being free at *both* a parent
and a descendant level simultaneously (the case of the previous section); with nothing free below
the parent, at any depth, there is no such degeneracy. Use this to align, say, a half-barrel as a single rigid body
while keeping its staves/modules/chips fixed to it, without needing to give them any DOF config.

## TPC Drift

The TPC drift calibration belongs to `TPC_envelope` and is not enabled by default. Enable both `VDRIFT` and `DRIFTOFF` with:

```json
{ "match": "TPC_envelope", "calib": { "type": "tpcvdrift" } }
```

Fix both with:

```json
{ "match": "TPC_envelope", "calib": { "type": "tpcvdrift", "fixed": true } }
```

## ITS Alignable Match Strings

The patterns below are based on the output from the `macro/printAlignableVolumes.C` running over Run3 TGeometry root file. 
The ITS envelope is the `ITS` alignable. The layer name is a path component, not a separate alignable entry. In the match patterns, layer IDs use `[0-2]` for the inner barrel and `[3-6]` for the outer barrel. Half-barrel, half-stave, module, and chip IDs in these paths are single digits, matched with `[0-9]`. Inner-barrel stave IDs are single digits; outer-barrel stave IDs can be one or two digits, so both exact-width forms are listed there. Avoid a trailing `.*` when matching a specific level: it also matches `/` and descendants.

### Inner Barrel (Layers 0–2)

The inner-barrel hierarchy has half-barrels, staves, and chips:

| Alignable level | `match` string |
| --- | --- |
| ITS envelope | `ITS` |
| Half-barrel | `ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]` |
| Stave | `ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]` |
| Chip | `ITS/ITSULayer[0-2]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]/ITSUChip[0-9]` |

### Outer Barrel (Layers 3–6)

The outer-barrel hierarchy adds half-staves and modules between staves and chips:

| Alignable level | `match` string |
| --- | --- |
| ITS envelope | `ITS` |
| Half-barrel | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]` |
| Stave, one-digit ID | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]` |
| Stave, two-digit ID | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9][0-9]` |
| Half-stave, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]/ITSUHalfStave[0-9]` |
| Half-stave, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9][0-9]/ITSUHalfStave[0-9]` |
| Module, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]/ITSUHalfStave[0-9]/ITSUModule[0-9]` |
| Module, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9][0-9]/ITSUHalfStave[0-9]/ITSUModule[0-9]` |
| Chip, under one-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9]/ITSUHalfStave[0-9]/ITSUModule[0-9]/ITSUChip[0-9]` |
| Chip, under two-digit stave | `ITS/ITSULayer[3-6]/ITSUHalfBarrel[0-9]/ITSUStave[0-9][0-9]/ITSUHalfStave[0-9]/ITSUModule[0-9]/ITSUChip[0-9]` |

For example, replace a `match` value with one of the strings above in a rule such as `{ "match": "...", "rigidBody": ["TX", "TY", "TZ"] }`. ITS3 uses different names, such as `ITS3CarbonForm`, and is not selected by these patterns.

## TRD Alignable Match Strings

TRD has a dummy `TRD_envelope` (not a geometry alignable), followed by alignable supermodules and chambers. The geometry log uses two digits for supermodule IDs and one digit each for stack and plane IDs:

| Level | `match` string |
| --- | --- |
| O2Align envelope (not in geometry log) | `TRD_envelope` |
| Supermodule | `TRD/sm[0-9][0-9]` |
| Chamber | `TRD/sm[0-9][0-9]/st[0-9]/pl[0-9]` |

## TOF Alignable Match Strings

TOF has a dummy `TOF_envelope` (not a geometry alignable), followed by alignable supermodules and strips. Both IDs in the geometry names are formatted as two digits:

| Level | `match` string |
| --- | --- |
| O2Align envelope (not in geometry log) | `TOF_envelope` |
| Supermodule | `TOF/sm[0-9][0-9]` |
| Strip | `TOF/sm[0-9][0-9]/strip[0-9][0-9]` |

The envelopes are not rigid-body alignable; rigid-body clauses on them are ignored. Calibration clauses can be configured on these synthetic volumes.