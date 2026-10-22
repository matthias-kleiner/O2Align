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

#include <algorithm>

#include <TGeoManager.h>

#include "Framework/Logger.h"
#include "DetectorsBase/GeometryManager.h"
#include "O2Align/AlignGeomUtils.h"

namespace o2::alignrs::geomutils
{

VolumePattern::VolumePattern(const std::string& pattern)
  : mDirect(pattern), mAnyPrefix(".*(?:" + pattern + ")")
{
}

bool VolumePattern::matches(const std::string& symName) const
{
  return std::regex_match(symName, mDirect) || std::regex_match(symName, mAnyPrefix);
}

std::string getAlignablePath(const std::string& symName)
{
  const auto* pne = gGeoManager->GetAlignableEntry(symName.c_str());
  if (pne == nullptr) {
    LOGP(fatal, "Symbolic volume '{}' has no corresponding alignable entry!", symName);
  }
  return pne->GetTitle();
}

int getPathDepth(const std::string& path)
{
  return static_cast<int>(std::count(path.begin(), path.end(), '/'));
}

const TGeoHMatrix& getGeometryMatrix(const std::string& symName)
{
  const auto* m = o2::base::GeometryManager::getMatrix(symName.c_str());
  if (m == nullptr) {
    LOGP(fatal, "Cannot get the matrix of {}", symName);
  }
  return *m;
}

TGeoHMatrix getOnPath(const PathMatrices& matrices, std::string path, bool includeSelf)
{
  auto stripLast = [&path]() {
    const auto pos = path.find_last_of('/');
    path.resize(pos == std::string::npos ? 0 : pos);
  };
  if (!includeSelf) {
    stripLast();
  }
  while (!path.empty()) {
    if (auto it = matrices.find(path); it != matrices.end()) {
      return it->second;
    }
    stripLast();
  }
  return {};
}

TGeoHMatrix getInitialGlobalDelta(const o2::detectors::AlignParam& par)
{
  auto delta = par.createMatrix();
  if (par.isGlobal()) {
    return delta;
  }
  const auto& gp = getGeometryMatrix(par.getSymName());
  return gp * delta * gp.Inverse();
}

PathMatrices collectInitialCumulativeDeltas(const std::vector<o2::detectors::AlignParam>* initial)
{
  std::vector<std::pair<std::string, const o2::detectors::AlignParam*>> byDepth;
  if (!initial) {
    return {};
  }
  for (const auto& par : *initial) {
    byDepth.emplace_back(getAlignablePath(par.getSymName()), &par);
  }
  std::stable_sort(byDepth.begin(), byDepth.end(), [](const auto& a, const auto& b) { return getPathDepth(a.first) < getPathDepth(b.first); });
  PathMatrices byPath;
  for (const auto& [path, par] : byDepth) {
    if (byPath.contains(path)) {
      LOGP(fatal, "Initial alignment contains {} more than once", par->getSymName());
    }
    byPath[path] = getInitialGlobalDelta(*par) * getOnPath(byPath, path, false);
  }
  return byPath;
}

} // namespace o2::alignrs::geomutils
