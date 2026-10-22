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

// Helpers shared between Detector::MP2AlignParams and the weak-modes diagnostics to turn CCDB
// o2::detectors::AlignParam objects into the cumulative global delta they represent relative to
// the ideal geometry, see the derivation in Detector.cxx.
#ifndef O2ALIGN_ALIGNGEOMUTILS_H
#define O2ALIGN_ALIGNGEOMUTILS_H

#include <string>
#include <regex>
#include <unordered_map>
#include <vector>

#include <TGeoMatrix.h>

#include "DetectorsCommonDataFormats/AlignParam.h"

namespace o2::alignrs::geomutils
{

/// A DOF-config-style symbolic-volume pattern: std::regex_match against either the full name or
/// the name with an implicit arbitrary prefix (so a suffix pattern can match at any hierarchy depth).
class VolumePattern
{
 public:
  explicit VolumePattern(const std::string& pattern);
  bool matches(const std::string& symName) const;

 private:
  std::regex mDirect;
  std::regex mAnyPrefix;
};

/// matrices keyed by the TGeo path of the alignable volume they belong to
using PathMatrices = std::unordered_map<std::string, TGeoHMatrix>;

/// TGeo path of an alignable volume, used to navigate the TGeo hierarchy by path prefixes
std::string getAlignablePath(const std::string& symName);

/// TGeo level of a path, used to order volumes parents first
int getPathDepth(const std::string& path);

/// L2G matrix of an alignable volume in the loaded (pre-aligned) geometry
const TGeoHMatrix& getGeometryMatrix(const std::string& symName);

/// Matrix of the TGeo volume `path` (if includeSelf) or else of its nearest TGeo ancestor present
/// in the map, the identity if there is none.
TGeoHMatrix getOnPath(const PathMatrices& matrices, std::string path, bool includeSelf);

/// global delta D^p equivalent to an initial AlignParam of either convention
TGeoHMatrix getInitialGlobalDelta(const o2::detectors::AlignParam& par);

/// cumulative global delta P^p of the initial alignment for every volume it contains, the objects
/// being accumulated parents first, like GeometryManager::applyAlignment applies them. An empty
/// map (identity) for every volume if initial is null.
PathMatrices collectInitialCumulativeDeltas(const std::vector<o2::detectors::AlignParam>* initial);

} // namespace o2::alignrs::geomutils

#endif
