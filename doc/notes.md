
# ITS_RB_v2.json : more or less good results but too much constrained.
Results are in /data0/devAlg/tests/solA1Iv2GC

Observations:
These modulations are the same for V10n and v2, and the same in pp B=0 and OO B+. So the fit neither made them nor tried to remove them: v2 has nothing free that is distributed along z. Below the stave and half-stave level everything is fixed, and a rigid stave can only produce structure that is constant or linear in z. To remove the shapes you need DOFs distributed along z. Each such DOF is measured locally, but freeing them opens a new family of weak modes, and those need their own global constraints.

## What the structures are

I fitted, for every sensor, the mean dy and dz on the basis {1, z_loc, tgλ, z_loc·tgλ}, using the pp B0 trees of V10n and v2. Then I averaged each coefficient over all staves at a given z.

**IB (well measured):**
- **dy vs Z, p.3:** the 3 cm sawtooth is an in-plane rotation of each chip (about the sensor normal), the same on all staves. At Lr0, dy/dz_loc goes −3.5, −1.5, −1.0, −0.8, −0.6, +0.2, +1.6, +2.5, +2.0 µm/cm for chips 0→8, so it grows with z. The statistical error is about 0.35 µm/cm per chip.
- **End chips:** they carry the deviations at the stave ends.
  - dy offsets reach −2.8 µm (Lr1) and −4.9 µm (Lr2) on chip 8.
  - dz offsets are +5.6 µm (chip 0) and +2.4 µm (chip 8) at Lr0.
  - On p.4 these appear as the ripple with the chip period plus the rise at |Z| > 10.
- **dz vs Z and vs tgλ:** a displacement normal to the sensor, Δn, gives dz ≈ Δn·tgλ, and tgλ reaches ±5 at the Lr0 end chips. So 2 µm of normal displacement there becomes 10 µm in dz. Chip normal positions are therefore well measured at the stave ends and unmeasured at z ≈ 0.

**OB (only from the PDF):**
- I can't quote numbers here: there are only about 200 hits per chip, and the outer modules never see tgλ ≈ 0, so the per-chip fit is ill-conditioned.
- The features sit at the module scale (21 cm) and mostly on the outer modules: |Z| > 21 cm at Lr3/4, |Z| > 52 cm at Lr5/6.
- dy rises linearly across the outer modules (Lr3: about 25 µm over 14 cm), which is a module rotation about the normal.
- dz vs tgλ is odd, up to +45 µm at tgλ = 2 on Lr3: a normal or z displacement growing towards the ends, i.e. a bow, or the module pitch/length.

## DOFs that address them

Use one representation per physical deformation, and only directions the residuals measure. "RY" below is the rotation about local Y, the sensor normal.

| level | free | why |
|---|---|---|
| IB chip | `TX, TZ, RY` | in-plane shifts and in-plane rotation, all measured directly (5 µm precision over 3 cm) |
| IB chip | `TY` with a soft `measurement` of about 5 µm | measured only at large \|z\| via dz = Δn·tgλ; the prior covers the central chips. If you free chip TY, **fix IB stave `RX`**: a stave tilt is then a TY ramp of its chips, and having both is the ITSBarFix degeneracy |
| OB module | `TX, TZ, RY`, plus `TY` with a prior of about 20–30 µm | same reasoning at the 21 cm scale; this is also where the per-layer z scale and module pitch live |
| OB half-stave | drop `RX` | for the same reason as IB stave `RX` |

The chip-row offsets from the previous analysis need OB chip `TX`; that can come in the same pass or the next one.

Each child DOF has a parent DOF of the same name, so the automatic mean constraint handles the 0th moment. The 1st moment along z is safe for the in-plane DOFs, because the in-chip and in-module slope is measured, so a TX ramp is distinguished from a parent RY. For TY it isn't, hence the prior and dropping `RX`.

## The weak modes this opens, and the constraints they need

Take a straight track from the vertex: at radius r it is at z = r·tgλ. Any displacement of the form below is absorbed by the track's own parameters.

| field | absorbed by |
|---|---|
| Δx_rφ = r·f(φ, tgλ) | the track φ |
| Δz = r·g(φ, tgλ) | the track tgλ |
| Δr = r·h(φ, tgλ) | curvature when B ≠ 0; invisible when B = 0 |
| Δr = a(tgλ), Δz = d(tgλ) | track z0; only the vertex breaks them |

The current mode set (`radialScaling`, `curl`, the harmonics, `globalTilt`) covers only the members that don't depend on tgλ. That was enough while every free DOF sat at stave level, centred at z = 0. Once DOFs are distributed in z, the tgλ-dependent members become accessible. The IB chip-TY dipole seen in v0 was one of them: its span went 76, 48, 26 µm for Lr0, Lr1, Lr2, roughly falling with radius, i.e. Δr ∝ z/r.

The modes I would add to `Volume::GlobalMode`, flagged at chip and module level:
- `zScaling`, Δz = εz. It already exists; it was skipped only because staves are centred at z = 0.
- `telescope`, Δz = εr.
- Radial modes Δr = ε·z/r and Δr = ε·z²/r². The second should **replace `Bowing` (Δr = εz²)**: εz² with one ε for all layers is not a null direction, while z²/r² is.
- Tangential modes r·Δφ = εz and r·Δφ = εz²/r.
- Optionally the cos φ / sin φ variants of each.

I derived this family by hand, for straight tracks plus curvature. Before trusting the set, run `method diagonalization` once on the extended configuration and check that every near-zero eigenvector is one of these modes. These modes are also invisible to checkResid, so the structures on p.3/4/7/8 lie outside them: constraining the modes cannot block the fix.

## Caveats

- Some patterns are the same on all staves, such as the Lr0 in-plane rotation growing with z and the OB end-module effects. Before freeing about 2000 module parameters, rule out the ideal geometry (chip and module pitch in the O2 model) and the cluster dictionary at large incidence. In Lr0 the in-chip slope and the tgλ dependence separate well, thanks to the vertex z spread, which favours geometry. In the OB that separation is weak.
- Lr5 dz is anomalously broad: RMS 350 µm against 220–245 µm on the other OB layers, and only a third of its points survive |dz| < 200 µm. This deserves a separate look, because it also inflates the Lr5 modulations.
- The all-free Mille data already contain the chip and module columns, so all of this needs only MilleSteer + pede.
==========================

v3 prel

The new modes and `ITS_RB_v3.json` are in place and the v3 fit runs. The result isn't usable yet: the 14 new constraints remove only part of a larger family of weak modes, and the fit fills the rest.

**What I added**
- **Modes:** three families in `Volume::GlobalMode`, `radialPolar`, `tangentialPolar` and `longitudinalPolar`. Each is a displacement r^p·(z/r)^k with p = 0, 1 and k ≤ 2, leaving out the members that already exist, which gives 14 equations. A finite-difference check confirms the rotation terms match the fields.
- **Diagnostics:** `WeakModesDiag` and `drawWeakModes.C` handle them; the macro has a new page plotting them against z/r. The docs are updated, `DOFConfig_rules.md` and section 6 of the handover.
- **`ITS_RB_v3.json`:**
  - IB chips and OB modules free in TX, TZ and RY (in-plane rotation), plus TY with priors of 5 µm (chips) and 30 µm (modules).
  - The RX tilt is dropped at stave and half-stave level.
  - Global modes apply over all of `ITS/.+`.
  - With the all-free Mille data this gives 10458 free parameters and 2010 constraints, numerically full rank.

**Two existing bugs found along the way**
- **pede crash:** with these long constraints (up to about 9000 terms) pede aborted with `free(): invalid next size`. A bounds-checked build of the same pede commit traced it to its constraint grouping (`grpcon`, `pede.f90:1902`). Its merging of constraint groups can leave a group empty and then write out of bounds. v2 already had one empty group and survived by luck. Workaround: every global-mode equation now carries the lowest free label in the tree with coefficient 1e-12, so pede processes those equations first. The coefficient can't be 0, because pede's reader drops zero terms. pede then forms a single group and runs.
- **`WeakModesDiag` fit B was void, in v2 as well:** the SVD it used (Eigen's `BDCSVD`) returned NaN, because `curl` and `globalTilt` exactly duplicate the rigid rotations. I switched it to `JacobiSVD`, which gives rank 42/45. The fit-B values (metric B) in earlier `weakModes.json` files should be ignored; metric A was correct.

**The v3 fit** (in a new directory, `/data0/devAlg/tests/solA1Iv3`)
- χ²/ndf is 1.025 against 1.12, and at equal cuts the objective drops 2.8%.
- The corrections are much larger than in v2:

| sensor displacement RMS | v2 | v3 |
|---|---|---|
| rφ, IB | 6–7 µm | 12–14 µm |
| rφ, OB | 20–45 µm | 33–85 µm |
| z, OB | 7–14 µm | 30–66 µm |

  Lr5 module TX reaches 864 µm.
- **Most of it is a "ray" mode.** I binned sensors in (φ, tgλ) cells and divided the displacement by r. 87% of dz/r is common to all seven layers, which is Δz = r·g(φ, tgλ): a change of tgλ per direction that no track from the vertex can see. rφ/r is likewise common within the IB and within the OB. v2 already had the z part, at about a third of the amplitude.
- The polar equations remove only the smooth, axisymmetric pieces of g. Once DOFs are spread along z, the invisible displacements form a whole function of (φ, tgλ), not a short list of global modes. checkResid can't settle this, because it is blind to these modes by construction. So I haven't run it on v3.

**Decision for you:** how to suppress that family. Options:
1. **Ray constraints:** per (φ, tgλ) cell, require the layer-averaged Δz/r and Δ(rφ)/r to vanish. This is the general form of what the polar equations do, and I can implement it as a new clause.
2. **Priors:** soft `measurement` priors on module and chip TX/TZ at the mechanical placement tolerance. This needs only a config change.
3. **Data that breaks the degeneracy:** cosmics or tracks not from the vertex.

Options 1 and 2 can be combined. Which way do you want to go?

Nothing is committed. All the scratch work (macros, the bounds-checked pede, the grouping replica) is in the session scratchpad.

===========================

v4 prel

I implemented both: per-cell ray constraints and priors on the fine DOFs. With them (config `doc/ITS_RB_v4.json`) the large invisible corrections of v3 are gone, at a cost of 0.6% in χ². I haven't run the v4 geometry through checkResid yet.

**1. Ray constraints: a new `rayModes` clause**
- It goes on the detector envelope, here `ITS`. All sensors below it are binned in cells of φ and t = z/r.
- Each cell gets one equation per component (rφ, z, r) saying the cell's mean angular displacement, seen from the vertex, is zero. Every free DOF on a sensor's path to the root enters through that sensor's linearised displacement.
- Each DOF's effect on a sensor is computed by finite differences of the same `AlignParam` matrix the corrections are applied with, so the rotation convention matches the fit by construction.
- Cells whose sensors span less than 1 cm in radius (a single layer) are skipped, since moving one layer is measurable.
- Every equation carries the pede grouping workaround from last time.
- With 24 φ × 16 t cells this gives 720 equations.
- **Don't combine it with global modes that already lie in the ray span** (`curl`, `radialScaling`, the harmonic families, `globalTilt`, the p = 1 polar terms): the constraint set becomes nearly degenerate. v4 keeps only `twist` and the p = 0 polar terms.
- It's documented in `DOFConfig_rules.md` under "Constraining Ray Modes".

**2. Priors (config only)**
- IB chips: TX, TZ 10 µm; TY 5 µm.
- OB modules: TX, TY, TZ 30 µm.

These widths are my assumptions, not survey numbers; please replace them if you have better ones.

**One extra config change: IB stave RY dropped.** In the first v4 pass the IB chip TX had an RMS of 32 µm against its 10 µm prior. That was the stave's in-plane rotation compensated by a TX ramp of its chips, leaving the chips rotated about fixed centres (correlation 0.88). This is the same double-representation problem as RX before. Without stave RY, chip TX drops to 15 µm RMS, chip RY carries the rotation, and χ² is slightly lower with 48 fewer parameters. The OB shows no such effect.

**Results** (`/data0/devAlg/tests/solA1Iv4`; the first pass is kept in `runA/`)

| | v2 | v3 | v4 |
|---|---|---|---|
| objective at iteration 3 | 3.424e7 | 3.328e7 | 3.348e7 |
| final χ²/ndf | 1.12 | 1.025 | 1.030 |
| sensor displacement RMS, all ITS | 59 µm | 104 µm | 42 µm |
| rφ RMS, OB | 20–45 µm | 33–85 µm | 13–20 µm |
| ray-common dz/r | 0.34e-4 | 1.01e-4 | 0.19e-4 |
| OB module TX RMS (max) | fixed | 73 (864) µm | 22 (91) µm |

pede accepts all 2655 constraints at full rank, as a single group.

**What's left:** a small pattern where Lr3 and Lr6 move opposite to each other in rφ/r and r (correlations −0.97 and −0.80), with a common part of about 0.7e-4 in rφ/r and 2 µm radially. It has zero mean per cell, so the ray constraint can't see it. It looks like a curvature-type mode that B ≠ 0 data absorb and B = 0 data should break; worth watching.

**Next step for you:** run `makeCDB` + checkResid on `solA1Iv4`. Because checkResid is blind to the ray modes, it measures only whether the newly freed module and chip DOFs fix the Z and tgλ modulations; it can't confirm the ray-mode fix itself.

Nothing is committed. The changes are in `Volume.h/.cxx`, `WeakModes.cxx`, `drawWeakModes.C`, the two docs, and the new `ITS_RB_v4.json`.


==========================

## ITS_RB_v4.json


The stronger IB patterns come from the IB chip DOFs that v4 freed: the fit asks for chip corrections that checkResid's residuals don't support. I haven't found the root cause yet, but I've ruled out several candidates.

**What got worse.** The same panels degrade in pp B0 and OO B+ (RMS over bins of the truncated means, pp):

| µm | V10n | v2 | v4 |
|---|---|---|---|
| Vtx dz vs φ (your p.2) | 4.3 | 1.8 | 5.4 |
| Vtx dz vs tgλ | 2.5 | 3.0 | 7.8 |
| IB dz vs φ, Lr0/1/2 | 2.0/1.3/2.5 | 1.4/0.8/1.2 | 2.1/1.4/2.0 |
| Lr0 dy vs Z | 1.8 | 1.7 | 4.6 |
| Lr6 dz vs Z | 9.0 | 8.7 | 15.8 |

The OB gets better (dy vs Z at Lr3–5, dz vs φ at Lr5), and so does Vtx dy vs Z (3.9 → 1.3 µm).

**The clearest symptom** is the 3 cm sawtooth in Lr0 dy vs Z on every chip (p.3). It comes from the newly freed chip in-plane rotation RY. The fit gave each IB layer's chips a common rotation:

| | fitted mean chip RY | in-chip slope in V10n |
|---|---|---|
| Lr0 | +5.0e-4 rad (5 µm/cm) | +0.09 µm/cm |
| Lr1 | +1.2e-4 rad | +0.02 µm/cm |
| Lr2 | −2.3e-4 rad | −0.53 µm/cm |

The checkResid data show no such rotation, so applying it creates the sawtooth. The chip-to-chip variation of RY does follow the measured slopes; only the common part is unsupported.

**Ruled out:**
- **Wrong rotation sign:** `AlignParam` rotations are right-handed and active, the same as the derivatives.
- **Wrong rotation pivot:** the tracking-to-chip jacobian moves the pivot correctly.
- **Wrong application:** the geometry dump `algdump.root` shows the chips rotated exactly as fitted (+4.94/+1.13/−2.20 against the fitted +4.97/+1.15/−2.26), with nothing compensating.
- **The IB part of the ray constraints:** that part really isn't a null direction, because the vertex z spread (4.9 cm pp, 5.4 cm OO) smears t across about four cells at the IB. I added `rMin`/`rMax` to `rayModes` and ran `ITS_RB_v5.json` (ray cells only beyond 10 cm, in `/data0/devAlg/tests/solA1Iv5`). χ² and the IB chip parameters came out essentially unchanged, so it isn't the cause.

**What remains:** the fit (GBL on the Mille data) and the reconstruction (the KF in checkResid) disagree about the IB at chip level. The fit sees in-chip slopes in V10n that the reconstruction doesn't. v1, which freed IB chip TX/TZ, failed in a similar way. The likely suspects:
- how the alignment workflow builds the IB measurement frame (`SensorITS::defineMatrixL2G` and `getFlatChipAlpha`) compared with the reconstruction's sensor matrix and tracking frame, and how the cluster coordinates enter it;
- the Mille data coming from different runs than the checkResid samples.

The workflow can't currently show its own per-point residuals; its debug stream writes only `MisRes`.

**Decision for you:**
1. **Diagnose:** add a per-point residual dump to the alignment workflow, run it on V10n for a few TFs, and compare the in-chip dy-vs-z slopes chip by chip with checkResid. That settles which side is wrong.
2. **Workaround:** keep v5's OB part and fix the IB chips as in v2, so you can check with checkResid that the IB recovers while the OB gains stay. Lr6 dz vs Z also got worse, so the OB module DOFs may need the same scrutiny.

Which do you want? The two can run in parallel.

The `rMin` change and these findings are in `DOFConfig_rules.md` and the handover. Nothing is committed.
