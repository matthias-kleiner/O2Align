# TRD tracklet points

A TRD point is the tracklet of one chamber (layer) of a track, defined as in the Kalman update of O2
(`TRDGlobalTrackingSpec.cxx`, `refitTRDTrack`), see `DetectorTRD::makeTrackletPoint`.

| quantity | definition |
|---|---|
| x, alpha | chamber T2L of the (possibly pre-aligned) volume applied to the calibrated tracklet, sector angle |
| y | `y_trklt − tilt·(z_row − z_trk)`, the tilt correction applied only if `|z_row − z_trk| < l` and `σz_trk² < l²/12`; then `− pull·getCorrYDy(snp)` if `TRDUseAngularPull >= 3` |
| z | `z_row + c·tgl`, `z_row` = pad-row centre in the tracking frame (moves with the chamber), `c` = non-row-crossing coefficient |
| covariance | `GPUTRDRecoParam::recalcTrkltCov(tilt, snp, l, pull)`: σz² = l²/12 (up to the tilt mixing), y-z correlation from the tilt |

`z_trk` is the z of the track at the tracklet X, `l` the pad-row length, `tilt = tan(pad tilting angle)`.

Parameters (`AlignParams`):
- `TRDNonRCCorrDzDtgl`: `c`; the default (< −990) is the O2 value `GPUTRDRecoParam::getZCorrCoeffNRC()` (1.4). Any other value overrides it.
- `TRDUseAngularPull`: as O2 `rec.trd.useAngularPull` (default 3): ≠ 0 increases σy with |pull|, ≥ 3 also shifts y. 4 is the same as 3 here.
- `TRDAngularPullNSigma`: tracklets with |pull| above it are not used; `TRDCorrDVT`, `extraClsErrYTRD/ZTRD`: ad-hoc corrections / errors.
- The tracklet error parameters of the reconstruction (`trkltResRPhiIdeal`, `trkltResVsTanPhiMisalign`, pile-up ranges) are the
  `GPU_rec_trd.*` configurable params: pass the same `--configKeyValues` as the reconstruction used.

The angular pull is `(dy_trklt + tilt·tgl·cdrHght − 3·snp/√(1−snp²)) / √getDyRes(snp)`. The y shift depends on the track
snp in the chamber frame, i.e. on the chamber orientation: it can partly absorb a chamber RZ, treat it as part of the
point definition, `TRDUseAngularPull=0` switches it (and the pull term in σy) off.

Not reproduced from O2: the occupancy of the chamber in σy and in the pull, the pile-up y correction and error.

## Derivatives

The z row of the derivatives is the usual one (z_row moves with the chamber). The tilt-corrected y also depends on `z_row`:
a displacement changing the z residual by dz changes the y residual by `−tilt·dz`, hence `der.row(0) −= tiltDz·der.row(1)`
in `buildPointGlobals`.

## Approximation: the track z in the tilt correction is a constant

Physically the tracklet y is the y of the hit at the pad-row centre of a tilted pad,
`y_trklt = y_trk + tilt·(z_row − z_trk)`, i.e. the prediction of `y_trklt` depends on the track z as well as on the track y.
The exact GBL measurement would be 2-dimensional with a projection matrix mixing the two local track offsets
(`y_pred = u − tilt·v`, `z_pred = v`), and the measurement `y_trklt − tilt·z_row`.

What is done instead (as in the Kalman update of O2): the corrected y, `y_trklt − tilt·(z_row − z_trk,KF)`, is a measurement of
the track y alone, with `z_trk,KF` the z of the track from the refit preceding the GBL fit, a fixed number. Consequences:
- the dependence of the prediction on the fitted z of the GBL track, `−tilt·δz` with `δz = z_GBL − z_trk,KF`, is not modelled;
  the error on y is `tilt·δz` (tilt ≈ 0.05 → 50 µm for δz = 1 mm), to be compared with σy of the tracklet;
- the dependence on the chamber position through `z_row` is modelled (`tiltDz`), so the alignment derivatives are right
  for the chamber displacement, but not for the way the displacement changes the fitted z.
`δz` is not measured; it is expected to be well below σy since the z of the track at the TRD is constrained by the inner
detectors and by the TRD z itself. The exact treatment needs the projection matrix above in `fillGBLPoints`.
