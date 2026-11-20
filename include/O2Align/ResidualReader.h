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

/// \file ResidualReader.h
/// \brief Reader of the TPC unbinned residuals of the async / calibration pass (o2tpc_residuals_*.root, TrackInterpolation): tree
/// "unbinnedResid" (residuals, per-track TrackDataCompact, CTPLumi) and tree "trackData" (TrackData with the ITS track), one entry per TF
/// in both. Per track the measured points = reference + residual: the TPC clusters (IDEAL transformation with the file's reference
/// vDrift / t0, i.e. no correction: the frame the correction maps are made in), then nExtDetResid points of the other detectors: TRD
/// (rows 160-165, skipped), TOF (row 170, sec = TOF sector, channel = pad in the sector), ITS (rows 180-186, channel = sensor), the
/// primary vertex (row 190: alpha = tgSlp * pi, x = channel * MaxVtxX / 0x7fff). The TF time is the CTPLumi orbit - 1 (the
/// residuals' timeMS is built from the TF counter). Decoding as alice-tpc-notes staticMapCreatorITSHelpers.C and gblfit's
/// io/CollisionReader (o2_macros), whose port this is.

#ifndef O2_ALIGN_RESIDUALREADER_H
#define O2_ALIGN_RESIDUALREADER_H

#include <array>
#include <memory>
#include <string>
#include <vector>
#include "ReconstructionDataFormats/Track.h"

class TFile;
class TTree;

namespace o2::alignrs
{

/// a TPC cluster of a residual track, in its sector frame, IDEAL (uncorrected) position
struct ResidualTPCCluster {
  int sector = 0;     ///< 0-35
  int row = 0;        ///< pad row 0-151
  float x = 0.f;      ///< the pad row's x [cm]
  float y = 0.f;      ///< measured y [cm] (reference + residual) = the pad-plane y
  float z = 0.f;      ///< measured z [cm]
  bool rejected = false; ///< flagged by TrackInterpolation::validateTrack
  int flags = -1;        ///< ClusterNative flags (UnbinnedResid::channel of a TPC residual: split pad 0x1, single pad 0x8, shared 0x10)
};

/// an ITS cluster, a TOF cluster or the primary vertex: a point in the frame (alpha, x); for ITS the frame is the sensor's (set by
/// DetectorITS from GeometryTGeo, not here), for TOF only alpha (the pad defines the position)
struct ResidualPoint {
  double alpha = 0.;
  double x = 0.;
  double y = 0.; ///< measured (reference + residual)
  double z = 0.;
  int sector = -1; ///< TOF: sector
  int channel = -1; ///< ITS: sensor; TOF: pad in the sector
  int layer = -1;   ///< ITS: layer 0-6
};

struct ResidualTrack {
  unsigned int orbit = 0; ///< orbit of the TF (CTPLumi orbit - 1)
  int sourceId = -1;      ///< GlobalTrackID source of the track the residuals were made of
  int filterFlag = -1;    ///< TrackInterpolation::validateTrack: 0 validated, > 0 the reason it failed, -1 not validated
  o2::track::TrackParCovD itsOuter; ///< TrackData::par = the ITS track at its OUTER end (r ~ 40 cm), with covariance
  int nClustersITS = 0;
  int nClustersTPC = 0;
  std::vector<ResidualTPCCluster> tpc;
  std::vector<ResidualPoint> its;
  std::vector<ResidualPoint> tof;
  bool hasVertex = false;
  ResidualPoint vertex;
  std::array<float, 4> multStack{}; ///< TPC occupancy per stack IROC, OROC1-3 (TrackDataCompact::getMultStack; 0 if not stored)
};

class ResidualReader
{
 public:
  ResidualReader();
  ~ResidualReader();
  /// the residual files: one .root file or a text file with one per line
  static std::vector<std::string> inputFiles(const std::string& input);
  /// open a file; false if it cannot be read (both trees with the same number of entries needed)
  bool open(const std::string& file);
  long getEntries() const { return mEntries; }
  /// all tracks of TF `entry` of the open file; false if it cannot be read
  bool readTimeFrame(long entry, std::vector<ResidualTrack>& tracks);
  float getRefVDrift() const { return mRefVDrift; }
  float getRefT0() const { return mRefT0; }

 private:
  struct Branches;
  std::unique_ptr<TFile> mFile;
  std::unique_ptr<Branches> mBranches;
  long mEntries = 0;
  float mRefVDrift = 0.f; ///< the open file's drift reference [cm/us]
  float mRefT0 = 0.f;     ///< [us]
};

} // namespace o2::alignrs

#endif
