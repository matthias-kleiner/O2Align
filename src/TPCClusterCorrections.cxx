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

#include "O2Align/TPCClusterCorrections.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <fmt/format.h>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/Logger.h"
#include "SpacePoints/TrackResiduals.h"
#include "TPCBase/Mapper.h"
#include "TPCCalibration/TPCScaler.h"

namespace o2::alignrs
{

// ---------------------------------------------------------------------------------------------------------------- TPCEdgeTerm

bool TPCEdgeTerm::loadMap(const std::string& file, int term, std::string& error)
{
  std::unique_ptr<TFile> input(TFile::Open(file.c_str()));
  auto* tree = input ? input->Get<TTree>("voxResTree") : nullptr;
  auto info = [&](const char* key) -> std::string {
    auto* object = tree ? tree->GetUserInfo()->FindObject(key) : nullptr;
    return object ? object->GetTitle() : "";
  };
  if (!tree || info("y2xBinning").empty() || info("z2xBinning").empty() || info("maxZ2X").empty()) {
    error = fmt::format("{} has no voxResTree with the binning in its user info", file);
    return false;
  }
  const std::string y2x = info("y2xBinning"), z2x = info("z2xBinning");
  if (y2x.find_first_not_of("0123456789") != std::string::npos || z2x.find_first_not_of("0123456789") != std::string::npos) {
    error = fmt::format("{}: only uniform binnings (a number of bins), got y/x \"{}\", z/x \"{}\"", file, y2x, z2x);
    return false;
  }
  if (info("y2xEdge") != "padplane") {
    error = fmt::format("{}: y/x grid edge \"{}\", expected \"padplane\"", file, info("y2xEdge"));
    return false;
  }
  mNY = std::stoi(y2x);
  mNZ = std::stoi(z2x);
  mMaxZ2X = std::stod(info("maxZ2X"));
  using TR = o2::tpc::TrackResiduals;
  TR::VoxRes* result = nullptr;
  tree->SetBranchAddress("voxRes", &result);
  auto& nodes = mNodes[term];
  nodes.assign(size_t(36) * mNX * mNY * mNZ * 3, 0.f);
  long nSet = 0;
  for (long entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    const int row = result->bvox[TR::VoxX], iy = result->bvox[TR::VoxF], iz = result->bvox[TR::VoxZ];
    if (result->bsec >= 36 || row >= mNX || iy >= mNY || iz >= mNZ || !(result->flags & TR::SmoothDone)) {
      continue;
    }
    const size_t node = ((size_t(result->bsec) * mNX + row) * mNY + iy) * mNZ + iz;
    for (int c = 0; c < 3; ++c) {
      nodes[3 * node + c] = result->DS[c];
    }
    ++nSet;
  }
  tree->ResetBranchAddresses();
  delete result;
  LOGP(info, "TPC sector-edge term: {}: {} smoothed nodes (grid {} x {} x {}, maxZ2X {})", file, nSet, mNX, mNY, mNZ, mMaxZ2X);
  if (nSet == 0) {
    error = fmt::format("{}: no smoothed node", file);
  }
  return nSet > 0;
}

bool TPCEdgeTerm::load(const std::string& mapR, const std::string& mapL, const std::string& table, const std::string& windows, double scaleUM, std::string& error)
{
  mScaleUM = scaleUM;
  if (!loadMap(mapR, 0, error) || !loadMap(mapL, 1, error)) {
    return false;
  }
  const auto& mapper = o2::tpc::Mapper::instance();
  mHalfWidth.assign(mNX, 0.);
  for (int row = 0; row < mNX; ++row) {
    mHalfWidth[row] = 0.5 * mapper.getNumberOfPadsInRowSector(row) * mapper.getPadRegionInfo(o2::tpc::Mapper::REGION[row]).getPadWidth();
  }
  std::ifstream in(table);
  if (!in) {
    error = fmt::format("cannot read the dipole table {}", table);
    return false;
  }
  for (std::string line; std::getline(in, line);) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream fields(line);
    double t = 0.;
    std::string token;
    if (!(fields >> t)) {
      continue;
    }
    fields >> token >> token; // idcA, idcC
    std::array<float, 36> v{};
    bool ok = true;
    for (int b = 0; b < 36 && ok; ++b) {
      ok = static_cast<bool>(fields >> token);
      v[b] = ok ? std::strtof(token.c_str(), nullptr) : NAN;
    }
    if (ok) {
      mTimeMS.push_back(t);
      mV.push_back(v);
    }
  }
  std::ifstream win(windows);
  if (!win) {
    error = fmt::format("cannot read the windows {}", windows);
    return false;
  }
  for (std::string line; std::getline(win, line);) {
    std::istringstream fields(line);
    int index = 0;
    double start = 0.;
    if (fields >> index >> start) {
      mWindowStartMS.push_back(start);
    }
  }
  std::sort(mWindowStartMS.begin(), mWindowStartMS.end());
  mWindowMean.assign(mWindowStartMS.size(), {});
  std::vector<std::array<int, 36>> counts(mWindowStartMS.size(), std::array<int, 36>{});
  std::vector<std::array<double, 36>> sums(mWindowStartMS.size(), std::array<double, 36>{});
  for (size_t i = 0; i < mTimeMS.size(); ++i) {
    const int w = window(mTimeMS[i]);
    if (w < 0) {
      continue;
    }
    for (int b = 0; b < 36; ++b) {
      if (std::isfinite(mV[i][b])) {
        sums[w][b] += mV[i][b];
        ++counts[w][b];
      }
    }
  }
  for (size_t w = 0; w < mWindowStartMS.size(); ++w) {
    for (int b = 0; b < 36; ++b) {
      mWindowMean[w][b] = counts[w][b] ? static_cast<float>(sums[w][b] / counts[w][b]) : 0.f;
    }
  }
  mDescription = fmt::format("maps {} (R) / {} (L), dipole table {} ({} rows), windows {} ({}), scale {} um", mapR, mapL, table, mTimeMS.size(),
                             windows, mWindowStartMS.size(), mScaleUM);
  if (mTimeMS.empty() || mWindowStartMS.empty()) {
    error = "empty dipole table or window list: " + mDescription;
    return false;
  }
  return true;
}

int TPCEdgeTerm::window(double timeMS) const
{
  const auto it = std::upper_bound(mWindowStartMS.begin(), mWindowStartMS.end(), timeMS);
  if (it == mWindowStartMS.begin()) {
    return -1;
  }
  const int w = static_cast<int>(it - mWindowStartMS.begin()) - 1;
  const double end = w + 1 < static_cast<int>(mWindowStartMS.size()) ? mWindowStartMS[w + 1] : mWindowStartMS[w] + WindowLengthMS;
  return timeMS < end ? w : -1;
}

std::array<float, 36> TPCEdgeTerm::scales(double timeMS) const
{
  std::array<float, 36> s{};
  if (mTimeMS.empty()) {
    return s;
  }
  auto it = std::lower_bound(mTimeMS.begin(), mTimeMS.end(), timeMS); // nearest table row
  size_t i = it == mTimeMS.end() ? mTimeMS.size() - 1 : static_cast<size_t>(it - mTimeMS.begin());
  if (i > 0 && std::abs(mTimeMS[i - 1] - timeMS) < std::abs(mTimeMS[i] - timeMS)) {
    --i;
  }
  const int w = window(timeMS);
  if (std::abs(mTimeMS[i] - timeMS) > 2000. || w < 0) {
    return s;
  }
  for (int b = 0; b < 36; ++b) {
    s[b] = std::isfinite(mV[i][b]) ? mV[i][b] - mWindowMean[w][b] : 0.f;
  }
  return s;
}

std::array<double, 3> TPCEdgeTerm::shape(int term, int sector, int row, float rowX, float y, float z) const
{
  std::array<double, 3> value{};
  if (sector < 0 || sector >= 36 || row < 0 || row >= mNX || !(rowX > 0.f) || !(mHalfWidth[row] > 0.)) {
    return value;
  }
  // bin centres: u = -1 + (i + 0.5) 2 / nY, z/x = (j + 0.5) maxZ2X / nZ; constant beyond the outermost centres
  const double u = y / mHalfWidth[row];
  const double zx = (sector >= 18 ? -z : z) / rowX;
  const double fu = std::clamp((u + 1.) * mNY / 2. - 0.5, 0., double(mNY - 1));
  const double fz = std::clamp(zx / mMaxZ2X * mNZ - 0.5, 0., double(mNZ - 1));
  const int iu = std::min(static_cast<int>(fu), mNY - 2), iz = std::min(static_cast<int>(fz), mNZ - 2);
  const double wu = fu - iu, wz = fz - iz;
  const auto& nodes = mNodes[term];
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 2; ++b) {
      const double w = (a ? wu : 1. - wu) * (b ? wz : 1. - wz);
      const size_t node = ((size_t(sector) * mNX + row) * mNY + iu + a) * mNZ + iz + b;
      for (int c = 0; c < 3; ++c) {
        value[c] += w * nodes[3 * node + c];
      }
    }
  }
  return value;
}

void TPCEdgeTerm::apply(int sector, int row, float rowX, float& x, float& y, float& z, const std::array<float, 36>& scales) const
{
  if (sector < 0 || sector >= 36) {
    return;
  }
  const int bR = sector;                                       // boundary s|s+1 (the sector's high-y edge)
  const int bL = (sector / 18) * 18 + (sector % 18 + 17) % 18; // boundary s-1|s (its low-y edge)
  const auto r = shape(0, sector, row, rowX, y, z);
  const auto l = shape(1, sector, row, rowX, y, z);
  const double fr = scales[bR] / mScaleUM, fl = scales[bL] / mScaleUM;
  x -= static_cast<float>(r[0] * fr + l[0] * fl);
  y -= static_cast<float>(r[1] * fr + l[1] * fl);
  z -= static_cast<float>(r[2] * fr + l[2] * fl);
}

// ---------------------------------------------------------------------------------------------------------------- TPCCentreTable

bool TPCCentreTable::load(const std::string& file)
{
  std::ifstream in(file);
  for (std::string line; std::getline(in, line);) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (line[0] == 'D' || line[0] == 'O' || line[0] == 'C' || line[0] == 'A') { // optional terms of the chain's tables: not supported
      ++mNIgnored;
      continue;
    }
    std::istringstream f(line);
    std::vector<double> t;
    for (double x; f >> x;) {
      t.push_back(x);
    }
    if (t.size() == 4) { // group class bin shift_um
      const int g = static_cast<int>(t[0]), c = static_cast<int>(t[1]), b = static_cast<int>(t[2]);
      if (g < 0 || g >= NGroups || c < -1 || c > 2 || b < 0 || b >= NBins) {
        continue;
      }
      const size_t k = (size_t(g) * NClasses + (c + 1)) * NBins + b;
      mShift[k] = t[3] * 1e-4;
      mValid[k] = 1;
      ++mNGroupEntries;
    } else if (t.size() == 5) { // sector side row bin shift_um
      const int s = static_cast<int>(t[0]), sd = static_cast<int>(t[1]), r = static_cast<int>(t[2]), b = static_cast<int>(t[3]);
      if (s < 0 || s >= NSectors || sd < 0 || sd > 1 || r < 0 || r >= NRows || b < 0 || b >= NBins) {
        continue;
      }
      if (mRowShift.empty()) {
        mRowShift.assign(size_t(NSectors) * 2 * NRows * NBins, 0.f);
        mRowValid.assign(mRowShift.size(), 0);
      }
      mRowShift[rowIndex(s, sd, r, b)] = static_cast<float>(t[4] * 1e-4);
      mRowValid[rowIndex(s, sd, r, b)] = 1;
      ++mNRowEntries;
    }
  }
  return isLoaded();
}

double TPCCentreTable::shift(int sector, int side, int row, int cls, double padsFromCentre) const
{
  const int b = static_cast<int>(padsFromCentre / 0.25);
  if (b < 0 || b >= NBins) {
    return 0.;
  }
  if (!mRowValid.empty() && sector >= 0 && sector < NSectors && row >= 0 && row < NRows && mRowValid[rowIndex(sector, side, row, b)]) {
    return mRowShift[rowIndex(sector, side, row, b)];
  }
  const int g = group(row);
  const size_t kc = (size_t(g) * NClasses + (cls + 1)) * NBins + b;
  if (cls >= 0 && mValid[kc]) {
    return mShift[kc];
  }
  const size_t ka = (size_t(g) * NClasses) * NBins + b;
  return mValid[ka] ? mShift[ka] : 0.;
}

// ---------------------------------------------------------------------------------------------------------------- TPCScalerIDC

struct TPCScalerIDC::Weights {
  o2::tpc::TPCScalerWeights weights;
};

TPCScalerIDC::TPCScalerIDC() : mWeights(std::make_unique<Weights>()) {}

TPCScalerIDC::~TPCScalerIDC() = default;

void TPCScalerIDC::loadScaler(double timeMS)
{
  for (size_t i = 0; i < mScalers.size(); ++i) {
    if (timeMS >= mScalers[i]->getStartTimeStampMS() && timeMS <= mScalers[i]->getEndTimeStampMS(o2::tpc::Side::A)) {
      mCurrent = i;
      return;
    }
  }
  const long key = static_cast<long>(timeMS) / 10000;
  if (std::find(mTried.begin(), mTried.end(), key) != mTried.end()) { // no object covers it: keep the current one (idcAt rejects the TF)
    return;
  }
  mTried.push_back(key);
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  if (!mWeightsInit) {
    mWeightsInit = true;
    if (auto* w = ccdb.getForTimeStamp<o2::tpc::TPCScalerWeights>("TPC/Calib/ScalerWeights", static_cast<long>(timeMS))) {
      mWeights->weights = *w;
    }
  }
  auto* tree = ccdb.getForTimeStamp<TTree>("TPC/Calib/Scaler", static_cast<long>(timeMS));
  if (!tree) {
    return;
  }
  auto scaler = std::make_unique<o2::tpc::TPCScaler>();
  scaler->setFromTree(*tree);
  const auto& weights = mWeights->weights;
  if (weights.isValid()) { // as TPCScalerSpec with --ion-drift-time -1, max time for the weights 500 ms
    scaler->setScalerWeights(weights);
    scaler->useWeights(true);
    const float tInt = std::min(std::abs(weights.mFirstTimeStampMS), 500.f);
    if (tInt > 0) {
      scaler->setIonDriftTimeMS(tInt);
    }
  }
  for (size_t i = 0; i < mScalers.size(); ++i) {
    if (mScalers[i]->getStartTimeStampMS() == scaler->getStartTimeStampMS()) {
      mCurrent = i;
      return;
    }
  }
  LOGP(info, "TPC scaler {:.0f} - {:.0f} ms (TF {:.0f} ms), integration {:.3f} ms, ion drift time {:.0f} ms, weights {}", scaler->getStartTimeStampMS(),
       scaler->getEndTimeStampMS(o2::tpc::Side::A), timeMS, scaler->getIntegrationTimeMS(), scaler->getIonDriftTimeMS(), weights.isValid());
  mScalers.push_back(std::move(scaler));
  mCurrent = mScalers.size() - 1;
}

bool TPCScalerIDC::idcAt(double timeMS, float& idc)
{
  loadScaler(timeMS);
  if (mScalers.empty()) {
    return false;
  }
  const auto& scaler = *mScalers[mCurrent];
  // no real object covers the time (end of run: a default object starting at 0): the production takes the CTP lumi, not available here
  if (scaler.getStartTimeStampMS() <= 0 || timeMS < scaler.getStartTimeStampMS() - 1e3 || timeMS > scaler.getEndTimeStampMS(o2::tpc::Side::A) + 1e3) {
    return false;
  }
  idc = 0.5f * (scaler.getMeanScaler(timeMS, o2::tpc::Side::A) + scaler.getMeanScaler(timeMS, o2::tpc::Side::C));
  return true;
}

} // namespace o2::alignrs
