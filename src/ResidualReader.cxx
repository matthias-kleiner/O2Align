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

#include <cmath>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"
#include "DataFormatsCTP/LumiInfo.h"
#include "SpacePoints/TrackInterpolation.h"
#include "SpacePoints/SpacePointsCalibParam.h"
#include "O2Align/ResidualReader.h"

namespace o2::alignrs
{

namespace
{
constexpr int RowVertex = 190;
float decodeResidual(short value) { return value * o2::tpc::param::MaxResid / 0x7fff; }
float decodeY(short value) { return value * o2::tpc::param::MaxY / 0x7fff; }
float decodeZ(short value) { return value * o2::tpc::param::MaxZ / 0x7fff; }
} // namespace

struct ResidualReader::Branches {
  TTree* residuals = nullptr;
  TTree* trackData = nullptr;
  std::vector<o2::tpc::UnbinnedResid>* res = nullptr;
  std::vector<o2::tpc::TrackDataCompact>* info = nullptr;
  o2::ctp::LumiInfo* lumi = nullptr;
  std::vector<o2::tpc::TrackData>* trk = nullptr;
};

ResidualReader::ResidualReader() = default;
ResidualReader::~ResidualReader() = default;

std::vector<std::string> ResidualReader::inputFiles(const std::string& input)
{
  std::vector<std::string> files;
  if (input.size() > 5 && input.substr(input.size() - 5) == ".root") {
    files.push_back(input);
    return files;
  }
  std::ifstream list(input);
  if (!list) {
    LOGP(fatal, "ResidualReader: cannot open the file list {}", input);
  }
  for (std::string line; std::getline(list, line);) {
    if (!line.empty() && line[0] != '#') {
      files.push_back(line);
    }
  }
  return files;
}

bool ResidualReader::open(const std::string& file)
{
  mBranches.reset();
  mFile.reset(TFile::Open(file.c_str()));
  mEntries = 0;
  if (!mFile || mFile->IsZombie()) {
    LOGP(error, "ResidualReader: cannot open {}", file);
    mFile.reset();
    return false;
  }
  auto br = std::make_unique<Branches>();
  br->residuals = mFile->Get<TTree>("unbinnedResid");
  br->trackData = mFile->Get<TTree>("trackData");
  if (!br->residuals || !br->trackData || br->residuals->GetEntries() != br->trackData->GetEntries()) {
    LOGP(error, "ResidualReader: no unbinnedResid / trackData with the same entries in {}", file);
    mFile.reset();
    return false;
  }
  if (br->residuals->SetBranchAddress("res", &br->res) < 0 || br->residuals->SetBranchAddress("trackInfo", &br->info) < 0 ||
      br->residuals->SetBranchAddress("CTPLumi", &br->lumi) < 0 || br->trackData->SetBranchAddress("trk", &br->trk) < 0) {
    LOGP(error, "ResidualReader: missing branch in {}", file);
    mFile.reset();
    return false;
  }
  mEntries = br->residuals->GetEntries();
  mBranches = std::move(br);
  mRefVDrift = 0.f, mRefT0 = 0.f;
  if (auto* records = mFile->Get<TTree>("records"); records && records->GetEntries() > 0) {
    records->SetBranchAddress("TPCVDriftRef", &mRefVDrift);
    records->SetBranchAddress("TPCDriftTimeOffsetRef", &mRefT0);
    records->GetEntry(0);
    records->ResetBranchAddresses();
  }
  LOGP(info, "ResidualReader: {} with {} TFs, reference vDrift {} cm/us, t0 {} us", file, mEntries, mRefVDrift, mRefT0);
  return true;
}

bool ResidualReader::readTimeFrame(long entry, std::vector<ResidualTrack>& tracks)
{
  tracks.clear();
  if (!mBranches || entry < 0 || entry >= mEntries) {
    return false;
  }
  auto& br = *mBranches;
  if (br.residuals->GetEntry(entry) <= 0 || br.trackData->GetEntry(entry) <= 0) {
    return false;
  }
  const auto& residuals = *br.res;
  const auto& infos = *br.info;
  const auto& trackData = *br.trk;
  if (infos.size() != trackData.size()) {
    return false;
  }
  const unsigned int orbit = br.lumi->orbit - 1; // as the map creation: the residuals' timeMS is built from the TF counter
  tracks.resize(infos.size());
  for (size_t i = 0; i < infos.size(); ++i) {
    const auto& info = infos[i];
    const auto& data = trackData[i];
    auto& track = tracks[i];
    track.orbit = orbit;
    track.sourceId = info.sourceId;
    track.filterFlag = info.filterFlag;
    track.nClustersITS = data.nClsITS;
    track.nClustersTPC = data.nClsTPC;
    {
      std::array<double, 5> par{};
      std::array<double, 15> cov{};
      for (int k = 0; k < 5; ++k) {
        par[k] = data.par.getParam(k);
      }
      for (int k = 0; k < 15; ++k) {
        cov[k] = data.par.getCov()[k];
      }
      track.itsOuter = o2::track::TrackParCovD(data.par.getX(), data.par.getAlpha(), par, cov);
    }
    const unsigned int first = info.idxFirstResidual;
    const unsigned int lastTPC = first + info.nResiduals;
    const unsigned int last = lastTPC + info.nExtDetResid;
    if (last > residuals.size()) {
      return false;
    }
    for (unsigned int k = first; k < lastTPC; ++k) {
      const auto& r = residuals[k];
      ResidualTPCCluster& cl = track.tpc.emplace_back();
      cl.sector = r.sec;
      cl.row = r.row;
      cl.x = r.getX();
      cl.y = decodeY(r.y) + decodeResidual(r.dy);
      cl.z = decodeZ(r.z) + decodeResidual(r.dz);
      cl.rejected = r.rejected;
    }
    for (unsigned int k = lastTPC; k < last; ++k) {
      const auto& r = residuals[k];
      ResidualPoint point;
      point.y = decodeY(r.y) + decodeResidual(r.dy);
      point.z = decodeZ(r.z) + decodeResidual(r.dz);
      if (r.row == RowVertex) {
        point.alpha = r.tgSlp * M_PI / 0x7fff;
        point.x = r.channel * o2::tpc::param::MaxVtxX / 0x7fff;
        track.vertex = point;
        track.hasVertex = true;
      } else if (r.isITS()) { // the frame (alpha, x) of the sensor is set by the ITS detector (GeometryTGeo)
        point.channel = r.channel;
        point.layer = r.getITSLayer();
        track.its.push_back(point);
      } else if (r.isTOF()) { // the position is that of the pad (the TOF detector)
        point.alpha = o2::math_utils::sector2Angle(r.sec % 18);
        point.sector = r.sec % 18;
        point.channel = r.channel;
        track.tof.push_back(point);
      }
    }
  }
  return true;
}

} // namespace o2::alignrs
