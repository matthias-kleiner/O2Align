# Handover: the OB dy patterns are a chip-row offset inside the OB modules (`solA2v10`)

Written for: the person continuing the ITS alignment tuning (knows O2Align, Millepede-II, and checkResid).
State: 2026-10-22. Everything below is read from `/data/devAlg/solA2v10` (config `doc/ITS_RB_v10.json`), compared
with `/data/devAlg/solA2v8ns` and with the older `res` trees in `/data/devAlg/tests/ITSBarFix` and
`/data/devAlg/tests/StavesOnly`. Nothing was re-run except checkResid post-processing and tree reading.

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

## 6. Tools

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
