# Handover: the OB dy patterns are a chip-row offset inside the OB modules (`solA2v10`)

Written for: the person continuing the ITS alignment tuning (knows O2Align, Millepede-II, and checkResid).
State: 2026-10-22. Everything below is read from `/data/devAlg/solA2v10` (config `doc/ITS_RB_v10.json`), compared
with `/data/devAlg/solA2v8ns` and with the older `res` trees in `/data/devAlg/tests/ITSBarFix` and
`/data/devAlg/tests/StavesOnly`. Nothing was re-run except checkResid post-processing and tree reading.
Sections 6–8 were added on 2026-10-23: the `solA2v11` result, the OB half-barrel weak mode, and the v12 config.

## 1. What v10 tried, and why it did nothing

v10 = v8ns + (a) `RZ` freed on the OB staves, with `RZ` added to the OB half-barrel `pinChildrenMeanNoCross`;
(b) per-layer `globalModes` blocks (`twist`, `bowing`, `zScaling`, the polar families) for `ITSULayer3..6`.
The aim was the large OB dy vs φ residuals and the dy vs Z structures at the OB layer edges.

Result: no change. Stat-weighted RMS of the per-bin medians (`macro/checkResidMedians.C`), OO / pp:

| histogram | v8ns | v10 |
|---|---|---|
| dy_Lr3 vs φ | 10.7 / 9.4 µm | 10.2 / 9.1 µm |
| dy_Lr4 vs φ | 10.9 / 9.0 | 10.5 / 8.9 |
| dy_Lr5 vs φ | 21.1 / 21.7 | 20.0 / 20.7 |
| dy_Lr6 vs φ | 4.6 / 6.0 | 4.4 / 6.0 |
| dy_Lr3..6 vs Z | 2.9 – 5.6 | identical to 0.05 µm |

- Pede did fit the new DOFs: the OB stave RZ come out at ~1.9 mrad RMS against errors of ~0.08 mrad, with no visible
  effect on dy. The fit used RZ to fake something else; it is not a real stave rotation.
- (b) was a configuration error: `zScaling` and the p = 1 polar modes lie in the span of `rayModes`, which bins the OB
  (`rMin` = 10 cm). `DOFConfig_rules.md` says not to combine them (near-degenerate constraint set).
- v10 and v8ns were **not** fitted on the same Mille data: v10 used `milleB{0,Minus,Plus}/mp2data_0.bin` (as v7), v8ns
  the 45-file `MilleV8/*` set. Their χ² (v10: 29.4 M / 31.9 M ndf = 0.92; v8ns: 26.0 M / 32.2 M = 0.81) are therefore not
  comparable, and the small differences in the table above mix config and data.

## 2. The OB dy vs φ pattern is not a stave pattern

With 100 φ bins (3.6°), the OB dy medians flip sign from bin to bin, e.g. at φ = 0.22, 0.28, 0.35:
Lr3 +33, −23, +11; Lr4 −12, +26, −15; Lr5 +52, −56, +33 µm (OO, v10). Lr3 and Lr5 are in phase, Lr4 opposite.
A harmonic analysis of these medians gives the same dominant harmonics in three layers with different stave counts:

| layer (staves) | strongest harmonics (amplitude µm) | harmonic n = nStaves |
|---|---|---|
| Lr3 (24) | n=40: 7.7, n=48: 7.0, n=12: 3.3 | 2.2 |
| Lr4 (30) | n=48: 10.8, n=1: 4.9, n=40: 3.4 | 0.8 |
| Lr5 (42) | n=48: 18.8, n=40: 16.3, n=12: 7.6 | 4.5 |
| Lr6 (48) | n=16: 3.5 | 0.6 |

So the structure is finer than a stave, about 1.5–2 cm in rφ (one chip width), aliased by the binning.

## 3. Root cause: the two chip rows of an OB module are offset in rφ

Per-sensor medians of dy from the `res` tree, grouped by `GeometryTGeo::getChipIdInModule` (0–6 = one row, 7–13 = the other):

| tree | layer | chips 0–6 | chips 7–13 |
|---|---|---|---|
| solA2v10, pp B0 (562862), all tracks | Lr3 | +6.4 … +8.7 | −7.0 … −9.0 |
| | Lr4 | +6.1 … +8.5 | −4.7 … −7.1 |
| | Lr5 | +2.5 … +5.3 | −3.3 … −5.6 |
| | Lr6 | +9.8 … +14.0 | −8.3 … −12.4 |
| ITSBarFix, pp B0, all tracks | Lr3 | +2.8 … +5.4 | −7.7 … −10.9 |
| | Lr6 | +4.2 … +8.9 | −7.5 … −12.0 |
| StavesOnly, pp B0, 500k tracks | Lr3 | −0.9 … +2.6 | −12.6 … −15.9 |

- The row gap is 10–23 µm, present in both half-staves, in every OB layer, and for three different alignments
  (ITSBarFix and StavesOnly are V10n-based, solA2v10 is MP2). It is therefore independent of the alignment.
- It is purely rφ: the same grouping in dz gives |row difference| < 2 µm in all OB layers.
- The two rows sit side by side in φ, so φ binning alternates between them, which produces section 2.
- The finest OB alignable in every config up to v10 is `ITSUModule` (14 chips, 2 rows). No module, half-stave or stave
  DOF can produce a row-against-row shift, which is why V10n and every MP2 variant show the same pattern.
- The gap is about half the ALPIDE rφ pitch (26.88 µm / 2 = 13.4 µm). Combined with the rows being mounted flipped by
  180°, this hints at a geometry or segmentation convention issue for the flipped row in O2 rather than a mechanical
  misalignment. Not investigated.

## 4. The dy vs Z edge structures are single chips in the end modules

Module index (`getModule`) × chip-in-module median dy, solA2v10 pp tree. Middle modules are flat within each row; the
end modules are not. Examples:

| layer, module | chip 3 (row A) | chip 10 (row B) | other end-module outliers |
|---|---|---|---|
| Lr6 m0 (z ≈ −63) | +27.6 | +5.3 | chip 4 +19.9, chip 13 −35.0 |
| Lr6 m6 (z ≈ +63) | +30.2 | +11.2 | chip 7 −17.0, chip 12 −19.3 |
| Lr3 m3 | +0.9 | −15.0 (chip 10) | chips 9–13 −3.8 … −15.0 |
| Lr5 m6 | — | — | chip 0 −9.2, chip 11 −15.4 |

The Lr6 spikes at z ≈ ±63 in the checkResid dy vs Z page are these chips. They are even in z (same sign at both
ends), consistent with the even-in-tgλ OB dy noted in `handover_patterns_afterFixX.md` (p.7). Whether they are
geometry or an incidence-angle effect (tgλ ≈ ±1.6 at these z) cannot be decided from the residuals; a fit with
chip-level DOFs is the test.

## 5. Next step: `doc/ITS_RB_v11.json`

v11 = v8ns (v8 without the `weightedMean` rule; `v8ns.json` itself is not committed) + OB chip TX:

```json
{ "match": "ITS/ITSULayer[3-6]/.../ITSUModule[0-9]+/ITSUChip[0-9]+", "rigidBody": ["TX"] },
{ "match": "ITS/ITSULayer[3-6]/.../ITSUModule[0-9]+/ITSUChip[0-9]+", "measurement": { "TX": [0.0, 0.002] } }
```

- Both v10 additions are dropped.
- TX is the chip's local x = rφ. The module TX is free, so the automatic children-mean constraint makes the chip TX
  deviations from the module mean; the 20 µm prior keeps the low-statistics end chips bounded.
- No new `MilleData` pass is needed: the data was recorded with `ITS_RB_all.json`, which frees all 6 DOFs of every OB
  chip (labels are present in `/data/devAlg/tests/mp2paramALLITS.txt`).
- Size: +23 688 parameters (6048 chips in Lr3/4, 17 640 in Lr5/6), nVarGlob ≈ 10.6 k → ≈ 34 k. With the current
  `matrix inversion` / full storage this is ~4.6 GB for the matrix and a much longer solve; switch to `minresqlp` if needed.
- Run it on the same Mille data as the run it is compared with.

What to check after the run:
1. The row gap of section 3 closes (rerun the per-chip grouping on the new `res` tree).
2. OB dy vs φ RMS of medians drops well below 10 µm and the n = 12/40/48 harmonics disappear.
3. The end-module chip spikes of section 4 go away (geometry) or stay (incidence-angle / cluster effect).
4. The fitted chip TX show the row structure (≈ ±7 µm between rows); if they do, report the half-pitch coincidence to
   the O2 ITS geometry/segmentation side.

## 6. v11 result (`/data0/devAlg/tests/solA2v11`, 2026-10-23)

pede ran on all 45 `milleB0/BPlus/BMinus` files (`mpSteerV11.txt`), while v10 used only `mp2data_0.bin` of each, so
v10 → v11 mixes config and data. The steering keyword `minresqlp` is not one of pede's method names
(`sparseMINRES`, `sparseMINRES-QLP`, …; pede accepts an 80 % fuzzy match). pede therefore ran its default, plain MINRES
with sparse storage (`METSOL = 4`, `MATSTO = 2` in `millepede.log`). MINRES reached its iteration limit of 2000
(`istop = 5`) in the last 4 pede iterations. The limit `itnlim=2000` is hard-coded in `pede.f90` and cannot be set from the steering file.

The checks of section 5:

| check | v10 | v11 |
|---|---|---|
| 1. checkResid row gap (pp) | 16–25 µm | 3–7 µm, same sign |
| 2. OB dy vs φ RMS of medians (OO / pp) | 9–21 µm | 3–6 µm |
| 2. harmonics n = 40 / 48 | 7–19 µm | ≤ 1.5 µm |
| 3. end-module dy spikes, e.g. Lr6 m6 chip 3 (pp / OO) | +46 / +42 µm | +4.5 / +6 µm, i.e. geometry |
| per-chip dy rms (n ≥ 150) | 15–34 µm | 7–16 µm |

4. Fitted chip TX:
   - In the geometry, row A (chips 0–6) has local x = −φ̂ and local z = −Z. Row B is flipped by 180° about the
     normal.
   - The children-mean constraint is written in local coordinates. It therefore gives mean(TX_A) = mean(TX_B) exactly,
     in every module.
   - Both rows fit the same local TX of about +4…+10 µm (median per layer). The row gap is thus a constant offset in
     the chip's own frame of about a quarter pitch, which appears as about half a pitch between the rows.
   - This points at a sensor/segmentation convention rather than at mounting. `SegmentationAlpide` has deliberately
     asymmetric edges (`PassiveEdgeReadOut` 1.2 mm, `PassiveEdgeTop` 37 µm), and `detectorToLocal` accounts for them.
     No bug was found there. This is the item to report to the O2 ITS geometry/segmentation side.

Problems left by v11:
- **TX gradient vs module RY.** A linear chip-TX gradient inside a module is degenerate with module RY (correlation
  −0.66…−0.88 in Lr4–6).
- **Stave-end chips.** These are Module0 chips 0/13, plus Module3 chips 6/7 in Lr3–4 and Module6 chips 6/7 in Lr5–6
  (|z| ≈ 40.7 / 72 cm). They take unphysical TX values up to 770 µm. The counts of chips with |TX| > 60 µm are
  Lr3 221, Lr5 112 and Lr6 463, all in the end modules. These chips are outside the checkResid acceptance.
- **dz vs Z, small-scale pattern (p.4).** These are random chip z offsets, because OB chip TZ is still fixed:
  - The within-module scatter is 5–17 µm MAD, with no common chip-in-module pattern (< 1.5 µm).
  - The offsets reproduce between pp and OO (r ≈ 0.75) and between v10 and v11 (r ≈ 0.93).
- **dz vs tgλ at large |tgλ| (p.8).** This is not a radial shift:
  - At fixed cluster Z (3 cm bins), Lr3 has d(dz)/d(tgλ) = +10…+29 µm per unit tgλ at |Z| > 22 cm, with the same sign at
    both ends, while dz/tgλ ≈ 0. A radial offset would make the two equal.
  - It is the same with V10n.
  - It is a nonlinear incidence (or z_v) effect, not alignment.
- **B+ vs B0 gap in dy vs φ (p.1).** In pp B0 only, there is an n = 2 term in every OB layer, with the same phase
  (cos ≈ +2…+3, sin ≈ −2.5…−4.7, i.e. 3–5.6 µm). It is absent in OO B+ and in the IB, and identical in v10 and v11. It
  depends on the dataset (field or period), not on the iteration. The test is separate pede fits on B0 and on B±. Lr4 OO
  also has an n = 1 term of about 4 µm.

Re-recording Mille data with the v11 alignment is not needed yet:
- Every open item above can be addressed on the existing data. The OB chip TZ derivatives are already in the
  `ITS_RB_all` recording.
- Linearization is not a limit at these sizes (tens of µm, ≲ 0.4 mrad).
- An incremental recording would turn the chip priors into priors on the increment, letting the badly constrained end
  chips drift further.
- Record once with the final alignment, as a convergence check; a converged fit should then return near-zero
  corrections.

## 7. The OB half-barrel "scissor" weak mode

The user's B=0 test extrapolates a 3-IB-point track to the OB clusters (lab X / Y vs φ). For v11 it shows a top/bottom
half-barrel pattern; for V10n it is flat.

The geometry diff vs V10n (mean shift of the OB sensor centres per layer and half-barrel; HB0 = top, y > 0) shows a
Y "scissor":

| v11 − V10n, dY (µm) | Lr3 | Lr4 | Lr5 | Lr6 |
|---|---|---|---|---|
| HB0 (top) | +5.3 | +4.0 | −0.9 | −4.5 |
| HB1 (bottom) | −6.2 | −2.8 | +1.0 | +4.0 |

- **Shape.** Top and bottom move in opposite directions, roughly linearly in R, with the sign changing at Lr5. The
  Lr3/4 halves open by about 10 µm while the Lr6 halves close by about 8.5 µm. X has a weaker version of the same
  pattern.
- **Growth.** The same pattern is about ±2.5 µm in v8ns, ±3–4 µm in v10 and ±4.5–6 µm in v11. It grows as more OB DOFs
  are freed.
- **Why `rayModes` misses it.** With `rMin` = 10, each (φ, t) cell equation only fixes Σ_s (e·d_s)/r_s over all OB
  layers together. For this pattern the sum is about 1% of Σ|·| (e.g. −13 against 1291 for dY on the top HB).
- **Why NoCross misses it.** `pinChildrenMeanNoCross` is blind to it by construction: a half-barrel translation is a
  cos φ / sin φ pattern of the stave DOFs (see `DOFConfig_rules.md`).
- **Why the data barely constrain it.** An OB-only rφ field a + b·R is absorbed by each track's d0 and φ. Only the IB
  plus vertex lever constrains it, and that lever carries little weight in the global fit. It is exactly what the B=0
  test measures.
- **Possible link to section 6.** A top/bottom-opposite Y shift projects mostly onto n = 2 in dy vs φ. It may be related
  to the B0-only n = 2 term of section 6; not verified.
- **Caveat.** Projecting the geometry diffs through an exact 3-point IB line gives changes of tens of µm, larger than
  the user's plots show. The IB's own 1–3 µm rφ changes between versions are amplified about 10× at Lr6. The user's
  test probably includes the vertex or weights the points; get its definition before comparing quantitatively.

## 8. Next step: `doc/ITS_RB_v12.json`

v12 = v11 plus the following changes. No new Mille data is needed; use the same 45 files (`mpSteerV11.txt`).

| change | rule |
|---|---|
| OB half-barrels: transported pin for the translations (section 7) | `"pinChildrenMean": ["TX","TY"]` + `"pinChildrenMeanNoCross": ["TZ","RY"]` |
| OB module RY fixed (TX-gradient degeneracy) | module `rigidBody: ["TX","TY","TZ"]` |
| OB chip TZ freed (dz pattern, p.4) | chip `rigidBody: ["TX","TZ"]`, priors TX and TZ 20 µm |
| stave-end chips: tighter priors | TX and TZ 5 µm on Module0 chips 0/13, Module3 chips 6/7 (Lr3–4), Module6 chips 6/7 (Lr5–6) |

- Using the transported pin is safe here. OB stave RZ is fixed and the staves are centred at z = 0. In the generated
  equations the stave TZ and RY coefficients are about 1e-5 and 5e-4, so no cheap DOF can pay for a half-barrel shift.
- `DOFConfig_rules.md` warns that putting both pin forms on one volume can make the constraints rank-deficient. Here
  they cover different DOFs, but check `pede.log` for a reported rank.
- A later `measurement` rule replaces all earlier priors of a volume, so the end-chip rules repeat TX and TZ.
- Generated steering, checked:
  - free parameters 34 146 → 56 142 (+23 688 chip TZ, −1 692 module RY);
  - RY constraints 488 → 200 (no module RY below the half-staves);
  - TZ constraints 495 → 2 187 (each module's TZ is now the mean of its chips).
- Expect pede to need about 25 GB (v11 peaked at 8.25 GB for about 28.6k parameters).

What to check after the run:
1. The section 7 scissor is gone from the geometry diff vs V10n, and the B=0 IB-extrapolation test is flat again.
2. Per-chip dz scatter in the OB drops; the small-scale dz vs Z pattern (p.4) disappears.
3. The stave-end chip TX/TZ stay within a few × 5 µm.
4. The row gap and the dy vs φ improvements of section 6 survive.

## 9. Tools

Per-chip grouping of the `res` tree (compiled with ACLiC in the O2 env; `CheckResidTypes.h` lacks the
`GlobalTrackID.h` include, add it before including the header):

```cpp
o2::base::GeometryManager::loadGeometry("o2sim_geometry-aligned.root"); // any ITS geometry: only numbering is used
auto g = o2::its::GeometryTGeo::Instance();
g->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
for (auto& p : tr->points) {
  if (p.lr < 3 || p.sens < 0) continue;
  h[p.lr]->Fill(g->getModule(p.sens) * 14 + g->getChipIdInModule(p.sens), p.dy); // then per-bin median
}
```

Use medians (`TH1::GetQuantiles`), not means; ignore bins with < 150 entries. The full pp tree (2 M tracks) takes a few
minutes. Harmonic analysis of the medians: a plain DFT of the 100 φ-bin medians of `checkResidMedians.C` output.

The scripts for sections 6–7 were kept only in a session scratchpad and are probably gone. They are easy to redo:
- per-sensor median dy/dz: `TTree::Draw("points.dy:points.sens")` into a TH2F, then per-bin quantiles;
- chip frames and centres: `getMatrixL2G(s)` applied to (0,0,0) and to the unit axes, for each geometry;
- geometry diffs: differences of those dumps between two `o2sim_geometry-aligned.root` files;
- labeled chip TX: parsed from `millepede_labeled.res`.
