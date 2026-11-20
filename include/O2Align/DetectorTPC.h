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

#ifndef O2_ALIGN_DETECTORTPC_H
#define O2_ALIGN_DETECTORTPC_H

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <memory>
#include <string>
#include <vector>

#include "DataFormatsTPC/VDriftCorrFact.h"
#include "O2Align/Detector.h"
#include "O2Align/Volume.h"

namespace o2::tpc
{
struct ClusterNative;
}

namespace o2::gpu
{
class GPUParam;
class TPCFastTransformPOD;
} // namespace o2::gpu

namespace o2::tpc::cluster_bias
{
struct ClusterBiasCorrection;
} // namespace o2::tpc::cluster_bias

namespace o2::alignrs
{
class TPCMap;
class TPCEdgeTerm;
class TPCCentreTable;
struct ResidualTPCCluster;

class DetectorTPC final : public Detector
{
 public:
  static constexpr int NStacks = 4; // IROC, OROC1, OROC2, OROC3

  DetectorTPC() : Detector(DetTPC) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData) final;
  bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) final;
  int prepareCosmicLeg(const o2::dataformats::TrackCosmicsExtended& cosm, int leg, const o2::track::TrackParCovD& ref, std::vector<FrameInfoExt>& frames) final;
  bool prepareResidualTrack(const ResidualTrack& trk, Track& resTrack) final;
  int getMinCosmicMeasurements() const final;

  /// both must be provided by the caller before the 1st prepareTrack call
  void setCorrMaps(const o2::gpu::TPCFastTransformPOD* maps) { mCorrMaps = maps; }
  void setTPCParam(const o2::gpu::GPUParam* par) { mTPCParam = par; }
  const o2::gpu::TPCFastTransformPOD* getCorrMaps() const { return mCorrMaps; }
  /// a shift added to every cluster after the transformation (a fitted map to correct with, or injected modes for tests)
  void setClusterShift(const TPCMap* shift) { mClusterShift = shift; }
  const o2::gpu::GPUParam* getTPCParam() const { return mTPCParam; }
  /// counts of the per-cluster bias correction (TPCMapParams.clusterBias*): corrected, rejected (d_pred below the cut), no prediction
  std::string getClusterBiasStats() const;
  /// the per-cluster bias correction is configured (loads the table at the first call)
  bool hasClusterBias()
  {
    if (!mClusterBiasInit) {
      initClusterBias();
    }
    return mClusterBias != nullptr;
  }
  /// a per-cluster table by the track is configured (bias / edge table or centre table): the tracks are prepared in two passes
  bool hasTrackTables() { return hasClusterBias() || mCentreTable != nullptr; }
  /// cosmics with the correction are prepared in two passes (AlignmentSpec::prepareCosmicTrack): 1 = no correction and no cluster
  /// measured within ClusterBiasPassOneEdgePads of the sector edge (the KF fit for the prediction), 2 = correction + rejection by the
  /// prediction (the references = the fitted states of pass 1); 0 = outside a cosmic (as 2)
  void setClusterBiasPass(int pass) { mClusterBiasPass = pass; }
  static constexpr float ClusterBiasPassOneEdgePads = 2.f;
  /// counts of the centre table: clusters shifted, no prediction (kept uncorrected)
  std::string getCentreTableStats() const;

  /// the IDC derivative map and its scale for the TF ((IDC - IDC of the base map) / IDC of the derivative), added to the base map's
  /// correction at the nominal position (TPCMapParams.derivativeMap); nullptr: off. Not owned
  void setDerivative(const o2::gpu::TPCFastTransformPOD* der, float scale)
  {
    mDerivative = der;
    mDerivativeScale = scale;
  }
  /// the sector-edge fluctuation term with the 36 boundary scales of the TF (TPCMapParams.edgeShape*); nullptr: off. Not owned
  void setEdgeTerm(const TPCEdgeTerm* term, const std::array<float, 36>& scales)
  {
    mEdgeTerm = term;
    mEdgeScales = scales;
  }
  /// the spacer-cross table by the track (TPCMapParams.centreTable); nullptr: off. Not owned
  void setCentreTable(const TPCCentreTable* table) { mCentreTable = table; }
  /// residual-ext with a table by the track: the reference track of the lookup (its state anywhere along the track, propagated to the
  /// clusters without material); nullptr = the reference fit itself (the clusters measured near the sector edge, residualTPCEdgePads, and
  /// near the centre line, tpcCentrePads, left out, no table). Not owned
  void setResidualReference(const o2::track::TrackParCovD* ref) { mResidualRef = ref; }

  /// The drift calibration already applied to the correction maps: the fitted VDRIFT correction is
  /// relative to it, hence it is needed to express the result wrt the reference drift velocity.
  void setVDrift(const o2::tpc::VDriftCorrFact& v) { mVDrift = v; }

  /// the drift frame of the base maps (TPCMapParams.driftNorm) for the TF being processed: common reference vD0 [cm/us], t00 [us], the
  /// window's vD slope per side (A, C); on = false: the clusters as the input transforms them
  struct DriftNorm {
    bool on = false;
    double vD0 = 0.;
    double t00 = 0.;
    std::array<double, 2> b{};
  };
  void setDriftNorm(const DriftNorm& n) { mDriftNorm = n; }
  const DriftNorm& getDriftNorm() const { return mDriftNorm; }
  /// residual-ext: the drift reference of the open residual file [cm/us], [us]
  void setResidualDriftReference(float vD, float t0)
  {
    mResRefVDrift = vD;
    mResRefT0 = t0;
  }
  /// z on the frame of n: from a cluster z made with the reference refVD [cm/us], refT0 [us]
  static float driftNormZ(const DriftNorm& n, float z, bool sideC, double refVD, double refT0)
  {
    const double s = sideC ? -1. : 1.;
    const double drift = 250. - std::abs(z);
    return static_cast<float>(z + s * ((1. - n.vD0 / refVD) * drift + n.vD0 * (n.t00 - refT0) - n.b[sideC ? 1 : 0] * drift));
  }
  const o2::tpc::VDriftCorrFact& getVDrift() const { return mVDrift; }

  const std::string& getTimeSlotsJson() const final;

  /// Millepede label of the drift calibration DOFs of a given time slot. The drift is a property of
  /// the whole detector, hence the DOFs belong to the envelope and the slot ID plays the role of
  /// its volume ID: every slot is calibrated independently.
  static Label getVDriftLabel(int slotID = 0) { return Label(DetTPC, slotID, false, true); }
  Volume* getEnvelopeVolume() const { return mEnvelope; }

  void initCalib();

  /// index of the stack the padrow belongs to, -1 if the padrow is out of range
  static int getStack(int padrow)
  {
    for (int i = 0; i < NStacks; ++i) {
      if (padrow <= sStackMinMaxRow[i].second) {
        return i;
      }
    }
    return -1;
  }

  /// distance (in padrows) to the closest min/max padrow of the stack containing the padrow
  static int getDistanceToStackEdge(int padrow)
  {
    const int st = getStack(padrow);
    if (st < 0) {
      return -999;
    }
    return std::min(padrow - sStackMinMaxRow[st].first, sStackMinMaxRow[st].second - padrow);
  }

 protected:
  Volume::Ptr buildHierarchy(Volume::SensorMapping& sensorMap) final;
  void onSlotChange(int slotID) final;

 private:
  bool isClusterAccepted(const o2::tpc::ClusterNative& cl, int row) const;
  /// x, y, z (nominal, on the maps' drift frame) -> the base map + the derivative x its scale (both at the nominal position) + the
  /// sector-edge term; rowX = the pad row's x
  void correctCluster(int sector, int row, float rowX, float& x, float& y, float& z) const;
  /// the centre table's shift [cm, signed: subtract from y] for the track at padY [cm, pad plane] of the cluster of class cls
  double centreShift(int sector, int row, int cls, double trackPadY) const;
  /// residual-ext lookup pass: dy [cm] = the centre table + the bias / edge table at the reference track trkAt (in the cluster's sector
  /// frame, at its x) for cluster c of track trk, corrected y yCorr, nominal z zNominal; false = dropped (the track is within
  /// clusterBiasRejectPads of the pad-plane edge)
  bool residualTables(const ResidualTrack& trk, const ResidualTPCCluster& c, float yCorr, float zNominal, const o2::track::TrackParCovD& trkAt, float& dy);

  static constexpr std::array<std::pair<int, int>, NStacks> sStackMinMaxRow = {
    std::pair<int, int>{0, 62}, std::pair<int, int>{63, 96}, std::pair<int, int>{97, 126}, std::pair<int, int>{127, 151}};

  const o2::gpu::TPCFastTransformPOD* mCorrMaps{nullptr};
  const o2::gpu::GPUParam* mTPCParam{nullptr};
  o2::tpc::VDriftCorrFact mVDrift{}; // drift calibration accounted by mCorrMaps
  Volume* mEnvelope{nullptr};        // fictitious top volume, owner of the drift calibration DOFs
  std::vector<Volume*> mSensors;     // sector volumes, indexed by sector
  const TPCMap* mClusterShift{nullptr}; // see setClusterShift, not owned
  DriftNorm mDriftNorm{};                // see setDriftNorm
  float mResRefVDrift{0.f};             // see setResidualDriftReference
  float mResRefT0{0.f};
  // per-cluster bias correction of the cosmics (TPCMapParams.clusterBias*), loaded at the first cosmic
  void initClusterBias();
  bool mClusterBiasInit{false};
  int mClusterBiasPass{0};
  std::shared_ptr<const o2::tpc::cluster_bias::ClusterBiasCorrection> mClusterBias; // null: off
  std::array<float, NStacks> mClusterBiasMult{-1.f, -1.f, -1.f, -1.f};               // multStack per stack for the occupancy term (< 0: off)
  std::array<uint64_t, 3> mClusterBiasCount{};                                         // corrected, rejected (d_pred below the cut), no prediction
  // map-chain corrections on top of the base map, see the setters (not owned)
  const o2::gpu::TPCFastTransformPOD* mDerivative{nullptr};
  float mDerivativeScale{0.f};
  const TPCEdgeTerm* mEdgeTerm{nullptr};
  std::array<float, 36> mEdgeScales{};
  const TPCCentreTable* mCentreTable{nullptr};
  std::array<uint64_t, 2> mCentreTableCount{}; // shifted, no prediction
  const o2::track::TrackParCovD* mResidualRef{nullptr};
};

} // namespace o2::alignrs

#endif
