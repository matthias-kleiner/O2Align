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
#include <cmath>
#include <format>
#include <sstream>

#include "MathUtils/Utils.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/DOFSet.h"
#include "O2Align/Detector.h"
#include "O2Align/Label.h"
#include "O2Align/TPCMap.h"

O2ParamImpl(o2::alignrs::TPCMapParams);

namespace o2::alignrs
{

bool TPCMapGrid::parse(const std::string& specification)
{
  if (specification.rfind("nodes/", 0) != 0) {
    return false;
  }
  mXNodes.clear();
  mZNodes.clear();
  std::stringstream tokens(specification.substr(6));
  for (std::string token; std::getline(tokens, token, '/');) {
    const auto equal = token.find('=');
    if (equal == std::string::npos) {
      continue;
    }
    const std::string key = token.substr(0, equal);
    const std::string value = token.substr(equal + 1);
    if (key == "ps") {
      mPerSector = std::stoi(value) != 0;
    } else if (key == "ny") {
      mNY = std::stoi(value);
    } else if (key == "zabs") {
      mZAbsolute = std::stoi(value) != 0;
    } else if (key == "x" || key == "z") {
      std::stringstream values(value);
      for (std::string node; std::getline(values, node, ',');) {
        (key == "x" ? mXNodes : mZNodes).push_back(std::stod(node));
      }
    }
  }
  return mXNodes.size() >= 2 && mZNodes.size() >= 2 && mNY >= 1;
}

void TPCMapGrid::basis(double x, double yx, double absZ, int nodes[8], double weights[8]) const
{
  auto cellUniform = [](double value, double low, double high, int n, int& index, double& fraction) {
    const double u = std::clamp((value - low) / (high - low) * (n - 1), 0., n - 1 - 1e-9);
    index = n > 1 ? static_cast<int>(u) : 0;
    fraction = n > 1 ? u - index : 0.;
  };
  auto cellNodes = [](double value, const std::vector<double>& grid, int& index, double& fraction) {
    if (value <= grid.front()) {
      index = 0;
      fraction = 0.;
    } else if (value >= grid.back()) {
      index = static_cast<int>(grid.size()) - 2;
      fraction = 1.;
    } else {
      index = static_cast<int>(std::upper_bound(grid.begin(), grid.end(), value) - grid.begin()) - 1;
      fraction = (value - grid[index]) / (grid[index + 1] - grid[index]);
    }
  };
  int ix = 0;
  int iy = 0;
  int iz = 0;
  double fx = 0.;
  double fy = 0.;
  double fz = 0.;
  cellNodes(x, mXNodes, ix, fx);
  cellUniform(yx, -mYXMax, mYXMax, mNY, iy, fy);
  cellNodes(mZAbsolute ? absZ : absZ / x, mZNodes, iz, fz);
  int n = 0;
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 2; ++b) {
      for (int c = 0; c < 2; ++c) {
        nodes[n] = node(std::min(ix + a, nX() - 1), std::min(iy + b, mNY - 1), std::min(iz + c, nZ() - 1));
        weights[n++] = (a ? fx : 1 - fx) * (b ? fy : 1 - fy) * (c ? fz : 1 - fz);
      }
    }
  }
}

bool TPCMap::init(const std::string& specification)
{
  return mGrid.parse(specification) && IDOffset + mGrid.nSlots() * mGrid.nNodes() < Label::ID_MAX;
}

void TPCMap::fillGlobals(const FrameInfoExt& frame, const DerivativeContext& ctx, int* labels, Eigen::Ref<Eigen::MatrixXd> der) const
{
  constexpr int NSectorsPerSide = 18;
  const int sector = static_cast<int>(frame.label.id());
  const double sectorAlpha = o2::math_utils::detail::sector2Angle<double>(sector % NSectorsPerSide);
  const double sign = std::cos(frame.alpha - sectorAlpha) > 0. ? 1. : -1.; // -1: the sector frame rotated by pi
  const double csp = std::sqrt((1. - ctx.snp) * (1. + ctx.snp));
  const double dydx = ctx.snp / csp;        // unchanged by the rotation
  const double dzdx = sign * ctx.tgl / csp; // in the sector frame
  const double x = sign * ctx.measX;        // the track at the measurement plane, in the sector frame
  const double y = sign * ctx.trkY;
  int nodes[8];
  double weights[8];
  mGrid.basis(x, y / x, std::abs(ctx.trkZ), nodes, weights);
  const int slotOffset = IDOffset + mGrid.slot(sector) * mGrid.nNodes();
  for (int component = 0; component < 3; ++component) {
    for (int k = 0; k < 8; ++k) {
      const int column = component * 8 + k;
      labels[column] = Label(Detector::DetTPC, slotOffset + nodes[k], false, true).rawGBL(component);
      der(0, column) = 0.;
      der(1, column) = 0.;
      if (component == 0) { // dX: the measurement plane moves, seen through the track slopes
        der(0, column) = -sign * dydx * weights[k];
        der(1, column) = -dzdx * weights[k];
      } else if (component == 1) {
        der(0, column) = sign * weights[k];
      } else {
        der(1, column) = weights[k];
      }
    }
  }
}

void TPCMap::writeParameters(std::ostream& os, double presigma) const
{
  static constexpr const char* ComponentNames[3] = {"dX", "dY", "dZ"};
  for (int slot = 0; slot < mGrid.nSlots(); ++slot) {
    for (int node = 0; node < mGrid.nNodes(); ++node) {
      const Label lbl(Detector::DetTPC, IDOffset + slot * mGrid.nNodes() + node, false, true);
      for (int component = 0; component < 3; ++component) {
        os << std::format("{:<10} {:>+15g} {:>+15g} ! V {} TPC map slot {} node {}\n", lbl.raw(component), 0.0, presigma,
                          ComponentNames[component], slot, node);
      }
    }
  }
}

bool TPCMap::decodeLabel(uint32_t raw, int& slot, int& component, int& node) const
{
  const uint32_t det = (raw & Label::DET_MASK) >> Label::DET_SHIFT;
  const bool calib = (raw & Label::CALIB_MASK) != 0;
  const bool sens = (raw & Label::SENS_MASK) != 0;
  const int64_t id = static_cast<int64_t>((raw & Label::ID_MASK) >> Label::ID_SHIFT) - 1 - IDOffset;
  if (det != Detector::DetTPC || !calib || sens || id < 0 || id >= mGrid.nSlots() * mGrid.nNodes()) {
    return false;
  }
  slot = static_cast<int>(id / mGrid.nNodes());
  node = static_cast<int>(id % mGrid.nNodes());
  component = static_cast<int>(raw & Label::DOF_MASK);
  return component < 3;
}

} // namespace o2::alignrs
