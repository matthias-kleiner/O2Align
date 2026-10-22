# Handout: the IB end-chip dy-vs-Z pattern on Lr0, and a proposed fix

Written for: the person continuing the ITS alignment tuning (knows O2Align, Millepede-II/GBL, and checkResid).
State: 2026-10-22. Based on `solA2v7`, `solA2v8`, `solA2v8ns`, `solA2v9` in `/data0/devAlg/tests/`
(configs `doc/ITS_RB_v7.json`, `v8.json`, `v8ns.json`, `v9.json`), all against the same `mp2data.bin`
recorded once with a wide-open `dofConfigJson`. Nothing here needed a new `ProcessData` pass once that
data was recorded: `v8`/`v8ns`/`v9` only differ in how `mp2con.txt` is regenerated standalone.

## 1. The pattern

`checkResid`'s per-chip `dy_Lr0_*` histograms on IB stave0 (both half-barrels, every config tried) show
chip0, chip7 and chip8 — the two chips closest to the stave ends plus their immediate neighbour —
sitting off the flat baseline the middle chips (3-6) reach, by up to a few µm, well outside their own
statistical error. This survives every rigid-body/global-mode retuning done so far (`v5bow` through
`v8`), including the cluster-x fix that removed the original sawtooth (see
`handout_patterns_afterFixX.md`). It is not a closure-test artifact: `AlignmentSpec::applyMisalignment`
is a no-op for ITS (only ITS3), so nothing is injected here to recover.

Motivating numbers for chip0/stave0/Lr0/HB0: a raw (pre-fit) `dy` vs `z_loc` profile restricted to
`sens==0` shows a slope corresponding to roughly a **1 mrad**-scale RY-like rotation, on top of a small
(~6 µm) TX-like offset — but chip0 only has about **1 800 hits** against chip4's **~38 000** (a factor
~21), since only the high-`|tgλ|` tail of tracks reaches that far along the stave.

## 2. What it isn't: diluted by the hierarchical mean constraint

The working hypothesis going in was that the automatic parent/children mean constraint
(`Volume::writeChildrenMeanConstraints`, equal 1/N split by default) ties chip0's free RY too tightly
to its high-statistics neighbours, so the fit reports a deceptively small, deceptively confident value
for chip0 instead of the larger correction its own (scarce) hits call for. This looked plausible enough
to implement and test two different fixes for it, both against the *same* `mp2data.bin`:

- **`AlignParams.useStatWeightedMeanConstraints` + a per-rule `weightedMean` clause**
  (`doc/DOFConfig_rules.md#weighting-the-mean-by-statistics`): weight the same constraint by each
  contributor's accumulated track count instead of splitting equally. Verified in `mp2con.txt` to give
  chip0/chip8 ≈1.3% of the stave0/Lr0/HB0 equation instead of the equal 11.1% — correctly reflecting
  their hit-count ratio. Result (`v8` vs the unweighted `v8ns`, otherwise identical config): checkResid
  `dy` on chip0 got *worse* in both half-barrels (HB0 0.00110→0.00170, HB1 0.00231→0.00306), chip8 mixed
  (worse in HB0, better in HB1), chip7 worse in both.
- **Removing the constraint outright** (`v9.json`: drop `RY` from the IB stave's own `rigidBody`, so
  with the parent DOF no longer free there is nothing left to tie the chips together — no constraint is
  written at all, each chip's RY becomes fully independent). Result: chip0 got worse in *both*
  half-barrels (HB0 0.00110→0.00192, HB1 0.00231→0.00398, see table below), chip8 worse in HB0, better
  in HB1 — essentially the same mixed, non-improving picture, despite chip0's fitted RY moving by a
  large amount (`+0.256 mrad` in `v8ns` → `−0.277 mrad` in `v9`) and its formal pede error *not*
  blowing up (`44.6 µrad` → `47.4 µrad`, barely changed).

| | stave0 RY | chip0 RY | chip0 err | HB0 chip0 dy | HB1 chip0 dy |
|---|---|---|---|---|---|
| `v8ns` (equal-weight mean) | −0.089 mrad | +0.256 mrad | 44.6 µrad | 0.00110±0.00124 | 0.00231±0.00116 |
| `v8` (stat-weighted mean) | −0.035 mrad | −0.405 mrad | 50.1 µrad | 0.00170±0.00156 | 0.00306±0.00101 |
| `v9` (no constraint at all) | fixed (0) | −0.277 mrad | 47.4 µrad | 0.00192±0.00146 | 0.00398±0.00117 |

The `v9` result is the decisive one: if the constraint had been suppressing a real, larger correction,
removing it entirely should have let chip0 settle on a value that *reduces* the checkResid residual.
Instead it moved to a different, still-small value that fits chip0's own sparse hits slightly better in
a way that doesn't generalize — the signature of an **underdetermined** parameter (too few independent
hits to resolve RY reliably on its own), not a **suppressed** one. No reweighting, pinning or freeing of
the existing per-chip/per-stave RY DOFs can fix this, because none of them add new information — they
only relabel how the same, already-fitted amount is split among correlated parameters.

## 3. Proposed fix: a pooled, non-hierarchical "end-chip tilt" parameter

What's missing is genuine extra statistics, which can only come from combining many chips' hits into
one parameter — the same trick `PVT/meanVertex` already uses (many tracks across many frames feed one
shared position label, bypassing the normal chip→stave→…→root derivative chain). Sketch:

1. **A new shared calibration-like DOF**, `endChipTilt`, one amplitude `eps` per (IB layer, stave end) —
   start with 2 (one per end, pooling layers) or up to 6 (per layer × per end, if the two layers'
   bowing-amplitude difference noted in `handout_patterns_afterFixX.md` suggests the end effect differs
   by layer too). Declared on a small virtual pseudo-volume outside the real ITS geometry tree, the same
   way the mean vertex is, so it never interacts with the existing hierarchical mean/global-mode
   machinery.
2. **An extra derivative term, not a new constraint.** Wherever `buildPointGlobals` computes a chip's
   own RY derivative column today (`Volume::getLeafRigidBodyDerivatives`), add one more term whenever the
   point's sensor is chip0 or chip8 of an IB stave: the same derivative value, attributed to the shared
   `endChipTilt` label for that (layer, end) instead of the chip's own RY label.
3. **No constraint needed.** Many independent chips (one per stave, dozens of staves) feed this one
   parameter, so Millepede's normal equations combine their statistics directly; the formal error on
   `eps` reflects genuine pooled precision instead of an artificial tie to neighbours.
4. **Avoid reintroducing the same degeneracy**: fix each end chip's own individual RY (don't free it,
   unlike `v8ns`/`v9`), keep chips 1-7 free as before. The stave's own RY DOF keeps its automatic
   constraint against those middle chips only — a basis shape (uniform across chip index) linearly
   independent from the new parameter's basis (nonzero only at the two ends), so the two don't overlap.
5. **This is also a clean diagnostic.** If the pooled `eps` comes back small and well-determined, the
   `dy`-vs-`Z` pattern is *not* a common rigid end-chip tilt and the search should move elsewhere (a
   sensor-level response effect, not representable by any rigid-body DOF). If it comes back large and
   significant, compose it into the end chips' final `AlignParam` at `writeMillepedeResults`/
   `writeAlignParams`, the way the mean vertex's position already bypasses the normal rigid-body
   composition path.

This is a real code feature (new DOF-set type, a derivative-injection point in
`DetectorITS`/`AlignmentSpec::buildPointGlobals`, output composition for the new label) — comparable in
scope to the mean-vertex position DOF or the stat-weighted-constraint work already done, not a
`dofConfigJson` tweak. Not implemented yet; this handout is the plan to pick up from.
