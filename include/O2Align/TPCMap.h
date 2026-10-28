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

#ifndef O2_ALIGN_TPCMAP_H
#define O2_ALIGN_TPCMAP_H

#include <array>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2::alignrs
{
struct FrameInfoExt;
struct DerivativeContext;

/// Configuration of the TPC correction map residual (TPCMap)
struct TPCMapParams : public o2::conf::ConfigurableParamHelper<TPCMapParams> {
  std::string grid = "";  // node grid "nodes/ps=<0|1>/ny=<n>/zabs=<0|1>/x=<x nodes>/z=<z nodes>" (see TPCMapGrid); empty: no map
  float presigma = 0.1f; // presigma of every map parameter [cm], i.e. a ridge towards the correction applied to the clusters
  bool gauge = false;    // constrain the global modes straight tracks of the TPC alone cannot see (TPCMap::writeConstraints)
  // shift of the TPC clusters after the transformation (TPCMap::clusterShift): a fitted map and / or global modes for tests
  std::string applyRes = ""; // millepede.res with a map on the same grid, applied to the clusters scaled by applyScale
  float applyScale = -1.f;   // -1: correct the clusters by that map, +1: inject it
  float injectRotZ = 0.f;    // [rad] a global rotation about z added to the clusters (dX = -t y, dY = t x in the sector frame)
  float injectScaleR = 0.f;  // a radial scale added to the clusters (dX = e x, dY = e y)
  float injectShiftZ = 0.f;  // [cm] a z shift added to the clusters
  float injectScaleZ = 0.f;  // a z scale added to the clusters (dZ = s z, global z)

  O2ParamDef(TPCMapParams, "TPCMapParams");
};

/// Node grid of the map residual: one set of nodes per TPC side (the sectors merged) or per sector (ps=1), non-uniform x nodes,
/// ny uniform y/x nodes over +-tan(10 deg), non-uniform |z| (zabs=1) or |z|/x nodes; trilinear basis.
class TPCMapGrid
{
 public:
  /// parse a "nodes/..." specification; false if invalid
  bool parse(const std::string& specification);

  int nX() const { return mXNodes.size(); }
  int nY() const { return mNY; }
  int nZ() const { return mZNodes.size(); }
  int nNodes() const { return nX() * nY() * nZ(); }
  int nSlots() const { return mPerSector ? 36 : 2; }
  int nParameters() const { return nSlots() * 3 * nNodes(); }
  int slot(int sector) const { return mPerSector ? sector : (sector >= 18 ? 1 : 0); }
  /// index of a node within its slot
  int node(int ix, int iy, int iz) const { return (ix * mNY + iy) * nZ() + iz; }
  double xNode(int ix) const { return mXNodes[ix]; }
  double yxNode(int iy) const { return mNY > 1 ? -mYXMax + 2 * mYXMax * iy / (mNY - 1) : 0.; }
  double zNode(int iz) const { return mZNodes[iz]; }
  bool zAbsolute() const { return mZAbsolute; }
  /// the 8 nodes (indices within the slot) and weights at (x, y/x, |z|); nodes beyond the grid are clamped to its edges
  void basis(double x, double yx, double absZ, int nodes[8], double weights[8]) const;

 private:
  bool mPerSector = false;
  bool mZAbsolute = true;
  int mNY = 3;
  double mYXMax = 0.17632698070846498; // tan(10 deg)
  std::vector<double> mXNodes;
  std::vector<double> mZNodes;
};

/// The residual of the TPC correction map on top of the correction applied to the clusters: the distortion d = (dX, dY, dZ) in the
/// sector frame, measured = true + d, trilinear on the nodes of TPCMapGrid. Its parameters are not DOFs of a volume: a cluster depends
/// on the 8 nodes around it only, hence the derivatives are supplied sparsely per point (fillGlobals) with labels of their own,
/// Label(DetTPC, IDOffset + slot * nNodes + node, sens = false, calib = true), the DOF index being the component (0 dX, 1 dY, 2 dZ).
class TPCMap
{
 public:
  static constexpr uint32_t IDOffset = 1u << 18; // above the IDs of the drift time slots (the envelope's calibration labels)
  static constexpr int NGlobals = 24;           // map columns of a point: 8 nodes x 3 components

  /// false if the specification is invalid
  bool init(const std::string& specification);
  const TPCMapGrid& getGrid() const { return mGrid; }

  /// Labels and derivatives of the prediction of a TPC point (rows: y, z of its measurement frame) wrt the map parameters, written
  /// to the columns [col, col + NGlobals) of der. The basis is evaluated at the track (ctx), not at the measurement, whose noise
  /// would otherwise bias the edge nodes. The measurement frame is the sector frame, or the sector frame rotated by pi for the
  /// upper leg of a cosmic (DetectorTPC::prepareCosmicExt): y and x flip, dy/dx is unchanged and dz/dx flips.
  void fillGlobals(const FrameInfoExt& frame, const DerivativeContext& ctx, int* labels, Eigen::Ref<Eigen::MatrixXd> der) const;

  /// Millepede "Parameter" records of all the map parameters (appended to the parameter file, without its header)
  void writeParameters(std::ostream& os, double presigma) const;

  /// Millepede "Constraint" records fixing the global affine modes of the map that straight TPC tracks do not see, common to both
  /// sides: rotation about z (dX = -y, dY = x), radial scale (dX = x, dY = y), z shift (dZ = 1) and z scale (dZ = z, global z),
  /// with the node coordinates in the sector frame. Needed when nothing outside the TPC (ITS, TOF) anchors the fit. The relative
  /// modes of the two sides are left free: cosmics crossing the central electrode see them.
  void writeConstraints(std::ostream& os) const;

  /// the map parameter of a raw Millepede label: false if the label is not one of the map
  bool decodeLabel(uint32_t raw, int& slot, int& component, int& node) const;

  /// load the node values of a fitted map on the same grid from a millepede.res; false if it has no parameter of this map
  bool loadValues(const std::string& resFile);
  /// value [cm] of component (0 dX, 1 dY, 2 dZ) of the loaded map at (x, y, z) of a sector frame
  double value(int sector, int component, double x, double y, double z) const;
  /// true if TPCMapParams asks for a cluster shift (applyRes loaded, or an injected global mode)
  bool hasClusterShift() const;
  /// the shift of a TPC cluster at (x, y, z) of its sector frame: applyScale x the loaded map + the injected global modes
  std::array<double, 3> clusterShift(int sector, double x, double y, double z) const;

 private:
  TPCMapGrid mGrid;
  std::vector<double> mValues; // node values [cm] of a loaded map, index (slot * 3 + component) * nNodes + node; empty if none
};

} // namespace o2::alignrs

#endif
