
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

