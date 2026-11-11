// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_ALIGN_DBG_H
#define O2_ALIGN_DBG_H

#include <vector>
#include <cmath>
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "O2Align/AlignmentTypes.h"

namespace o2::alignrs
{

/// Compact debug record of a GBL point, written for the tracks selected by Params::debugTracksFraction.
/// Residuals are measurement - track in the tracking frame (cm), the kinks are in the tracking basis (rad).
struct DbgPoint {
  enum : uint16_t {
    HasMeas = 0x1 << 0,   // the point carries a measurement
    HasKink = 0x1 << 1,   // the point carries a scatterer (the kinks are defined)
    IsVertex = 0x1 << 2,  // the point is the primary vertex
    ITSOverlap = 0x1 << 3 // the point is on an ITS overlap
  };
  int8_t det = -1;   // detector index (Detector::DetIdx, 0=PVT), -1 for the scatterer-only points
  int16_t sens = -1; // sensor ID within the detector
  int8_t lr = -1;    // detector-specific layer-like index
  uint16_t flags = 0;
  float x = 0.f;     // X of the point in the tracking frame
  float alpha = 0.f; // rotation angle of the tracking frame
  float y = 0.f;     // measured Y in the tracking frame
  float z = 0.f;     // measured Z
  float sigY = 0.f;  // cluster errors
  float sigZ = 0.f;
  float sigYZ = 0.f; // cluster Y-Z covariance (cm^2)
  float snp = 0.f;   // track snp at the point
  float dyPre = 0.f; // residuals wrt the linearization track, as fed to GBL (include the injected misalignment, if any)
  float dzPre = 0.f;
  float dyFit = 0.f; // residuals wrt the GBL fit
  float dzFit = 0.f;
  float dyUnb = 0.f; // unbiased (leave-one-out) residuals: wrt the GBL fit without this measurement
  float dzUnb = 0.f;
  float sigYUnb = 0.f; // errors of the unbiased residuals (cluster and track errors), 0 if not computed
  float sigZUnb = 0.f;
  float kinkY = 0.f; // GBL kinks at the scatterer (snp-like, tgl-like)
  float kinkZ = 0.f;
  float pullKinkY = 0.f; // kink / its fitted error
  float pullKinkZ = 0.f;
  // lab-frame coordinates of the measured point (x, y) of the tracking frame rotated by alpha, Z is the same in both frames
  float getLabX() const { return x * std::cos(alpha) - y * std::sin(alpha); }
  float getLabY() const { return x * std::sin(alpha) + y * std::cos(alpha); }
  float getLabR() const { return std::sqrt(x * x + y * y); }
  bool hasMeas() const { return flags & HasMeas; }
  bool hasKink() const { return flags & HasKink; }
  ClassDefNV(DbgPoint, 2)
};

/// Debug record of the cosmic the track was made of (cosmic-ext input): the quantities of the input the selection can cut on.
/// The reference is the track the legs were prepared with: `polished` if available, otherwise the matcher's track.
struct DbgCosmic {
  bool valid = false;          // the track is a cosmic and this record is filled
  bool hasTOFTime = false;     // a top-bottom TOF flight pair gave the time
  bool polished = false;       // the reference is the refit at the TOF time
  uint8_t flags = 0;           // TrackCosmics flags: 1 = NearBeamBottom, 2 = NearBeamTop
  float scoreTOFPair = -1.f;   // TrackCosmicsExtended::scoreTOFPair
  float scoreTOFReversed = -1.f;
  float chi2MatchPolished = -1.f;
  float chi2Match = 0.f;       // matcher: top/bottom match
  float chi2Refit = 0.f;       // matcher: global refit
  float timeTOFMUS = -1.f;     // TOF time (valid if hasTOFTime)
  float timeMUS = 0.f;         // matcher's time and its error
  float timeErrMUS = 0.f;
  float p = 0.f;               // reference track: momentum, pt, tgl, phi, distance to the beam line and z there
  float pt = 0.f;
  float tgl = 0.f;
  float phi = 0.f;
  float d0 = 0.f;
  float z0 = 0.f;
  float tpcOccBef = -1.f;      // TPC occupancy before / after the cosmic (-1: not in the input)
  float tpcOccAft = -1.f;
  // input per leg (0 = bottom, 1 = top): attached TPC clusters (all, on the A side, rows >= 63, 97, 127), ITS / TRD / TOF hits
  short nTPC[2] = {0, 0};
  short nTPCA[2] = {0, 0};
  short nTPCRow63[2] = {0, 0};
  short nTPCRow97[2] = {0, 0};
  short nTPCRow127[2] = {0, 0};
  short nITS[2] = {0, 0};
  short nTRD[2] = {0, 0};
  short nTOF[2] = {0, 0};
  short nPntTop = 0;           // points of the top leg in the fit (the first ones of DbgTrack::points with a detector)
  ClassDefNV(DbgCosmic, 1)
};

/// Debug record of a track: its parameters at the reference point (the innermost one) and all its GBL points
struct DbgTrack : public o2::track::TrackPar {
  o2::dataformats::GlobalTrackID gid{};
  FitInfo kfFit{};       // KF refit
  FitInfo gblFit{};      // GBL fit (of the whole composed trajectory for the multi-track PV constraint)
  bool accepted = false; // the GBL fit passed the chi2 cut and the track went to the output
  bool composed = false; // the track was fitted in a composed trajectory with the common vertex
  float bz = 0.f;        // nominal field, kG
  uint32_t tfCounter = 0;
  uint32_t firstTForbit = 0;
  int runNumber = 0;
  std::vector<DbgPoint> points;
  DbgCosmic cosm;        // cosmic-ext input only (cosm.valid)
  ClassDefNV(DbgTrack, 2)
};

} // namespace o2::alignrs
#endif
