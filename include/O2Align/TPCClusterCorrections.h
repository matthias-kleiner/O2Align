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

/// \file TPCClusterCorrections.h
/// \brief Per-cluster TPC corrections of the map chain on top of the base map (the same ones the chain's DCA validation applies,
/// o2_macros gblfit collisionDcaSeries): the sector-edge fluctuation term, the spacer-cross (centre) table by the track position and the
/// IDC of a TF for the derivative map.

#ifndef O2_ALIGN_TPCCLUSTERCORRECTIONS_H
#define O2_ALIGN_TPCCLUSTERCORRECTIONS_H

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace o2::tpc
{
class TPCScaler;
}

namespace o2::alignrs
{

/// The sector-edge fluctuation term: per sector two edge-shape maps (its high-y boundary s|s+1 = "R", its low-y boundary s-1|s = "L";
/// voxResTree DS in cm per <scale> um of the boundary's dipole) times the boundary's dipole at the TF time minus its mean over the map
/// window (the window maps hold the window's mean edge state):
///   delta = shapeR (v_R(t) - <v_R>_w) / scale + shapeL (v_L(t) - <v_L>_w) / scale,   x, y, z -= delta (measured -> true),
/// applied to the cluster after the map. The shapes are interpolated bilinearly in u = y / (half width of the pad row's pad plane) and the
/// side-oriented z/x between the bin centres (u: nY uniform bins over [-1, 1], z/x: nZ uniform bins over [0, maxZ2X]; from the file's user
/// info), constant beyond the outermost centres. Dipole table "time_ms idcA idcC v_0 .. v_35" (nan = masked -> 0); a table row within 2 s
/// of the TF time is needed, else no term. Windows "<index> <start ms> ..." (a window ends where the next starts, the last after 5 min).
class TPCEdgeTerm
{
 public:
  /// false if a map, the table or the windows cannot be read (error filled)
  bool load(const std::string& mapR, const std::string& mapL, const std::string& table, const std::string& windows, double scaleUM, std::string& error);
  bool isLoaded() const { return !mNodes[0].empty(); }
  /// the 36 boundaries' (v - <v>_window) [um] at the TF time [ms] (all 0 without a table row or a window)
  std::array<float, 36> scales(double timeMS) const;
  /// x, y, z (the cluster after the map, in its sector frame) corrected in place; rowX = the pad row's x
  void apply(int sector, int row, float rowX, float& x, float& y, float& z, const std::array<float, 36>& scales) const;
  const std::string& describe() const { return mDescription; }

 private:
  bool loadMap(const std::string& file, int term, std::string& error);
  int window(double timeMS) const;
  std::array<double, 3> shape(int term, int sector, int row, float rowX, float y, float z) const;
  int mNX = 152, mNY = 20, mNZ = 40;
  double mMaxZ2X = 2.;
  std::array<std::vector<float>, 2> mNodes; ///< [term][((sector * nX + row) * nY + iy) * nZ + iz][3]
  std::vector<double> mHalfWidth;           ///< per pad row: half width of the pad plane [cm]
  std::vector<double> mTimeMS;              ///< table rows
  std::vector<std::array<float, 36>> mV;    ///< table: v per boundary [um], nan = masked
  std::vector<double> mWindowStartMS;
  std::vector<std::array<float, 36>> mWindowMean;
  static constexpr double WindowLengthMS = 300000.;
  double mScaleUM = 100.;
  std::string mDescription;
};

/// The spacer-cross (sector centre) table, by the TRACK position: a y shift away from the sector centre line per 0.25-pad bin of the
/// track's distance from y = 0 (24 bins, out to 6 pads), lines "sector side row bin shift_um" (side 0 = -y, 1 = +y; used first where
/// present) and "group class bin shift_um" (group 0 rows 0-9, 1 10-24, 2 25-62, 3 63-151; class -1 = all, 0 normal, 1 split pad / shared,
/// 2 single pad on the pad centre; a missing class falls back to -1). The cluster gets y -= sign(track y) x shift. The optional drift /
/// occupancy / class x inclination / angle terms of the chain's tables (lines starting with D, O, C, A) are not supported (counted).
class TPCCentreTable
{
 public:
  bool load(const std::string& file);
  bool isLoaded() const { return mNGroupEntries + mNRowEntries > 0; }
  /// shift [cm] away from the centre line at |track y| = padsFromCentre [pads] on side (0 = -y, 1 = +y) of sector, cluster class cls
  double shift(int sector, int side, int row, int cls, double padsFromCentre) const;
  int nGroupEntries() const { return mNGroupEntries; }
  int nRowEntries() const { return mNRowEntries; }
  int nIgnoredTerms() const { return mNIgnored; }

 private:
  static constexpr int NGroups = 4;
  static constexpr int NClasses = 4; // class -1 (all), 0, 1, 2
  static constexpr int NBins = 24;
  static constexpr int NSectors = 36;
  static constexpr int NRows = 152;
  static int group(int row) { return row < 10 ? 0 : (row < 25 ? 1 : (row < 63 ? 2 : 3)); }
  static size_t rowIndex(int sector, int side, int row, int bin) { return ((size_t(sector) * 2 + side) * NRows + row) * NBins + bin; }
  std::array<double, NGroups * NClasses * NBins> mShift{};
  std::array<char, NGroups * NClasses * NBins> mValid{};
  std::vector<float> mRowShift; ///< [sector][side][row][bin], cm
  std::vector<char> mRowValid;
  int mNGroupEntries = 0;
  int mNRowEntries = 0;
  int mNIgnored = 0;
};

/// The IDC of a TF for the derivative map, as the production (TPCScalerSpec, --lumi-type 2): the mean of the A and C scalers of
/// TPC/Calib/Scaler at the TF time with TPC/Calib/ScalerWeights (ion drift time capped at 500 ms). The scaler objects (one per ~5 min,
/// overlapping) are fetched from the CCDB when the TF time leaves the ones in memory.
class TPCScalerIDC
{
 public:
  TPCScalerIDC();
  ~TPCScalerIDC();
  /// the IDC at timeMS [ms]; false if no scaler object covers it (the production falls back to the CTP lumi: skip the TF)
  bool idcAt(double timeMS, float& idc);

 private:
  void loadScaler(double timeMS);
  std::vector<std::unique_ptr<o2::tpc::TPCScaler>> mScalers;
  size_t mCurrent = 0;
  std::vector<long> mTried; ///< 10-s keys of the time stamps already asked for
  bool mWeightsInit = false;
  struct Weights;
  std::unique_ptr<Weights> mWeights;
};

} // namespace o2::alignrs

#endif
