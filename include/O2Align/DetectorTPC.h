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
  /// cosmics with the correction are prepared in two passes (AlignmentSpec::prepareCosmicTrack): 1 = no correction and no cluster
  /// measured within ClusterBiasPassOneEdgePads of the sector edge (the KF fit for the prediction), 2 = correction + rejection by the
  /// prediction (the references = the fitted states of pass 1); 0 = outside a cosmic (as 2)
  void setClusterBiasPass(int pass) { mClusterBiasPass = pass; }
  static constexpr float ClusterBiasPassOneEdgePads = 2.f;

  /// The drift calibration already applied to the correction maps: the fitted VDRIFT correction is
  /// relative to it, hence it is needed to express the result wrt the reference drift velocity.
  void setVDrift(const o2::tpc::VDriftCorrFact& v) { mVDrift = v; }
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

  static constexpr std::array<std::pair<int, int>, NStacks> sStackMinMaxRow = {
    std::pair<int, int>{0, 62}, std::pair<int, int>{63, 96}, std::pair<int, int>{97, 126}, std::pair<int, int>{127, 151}};

  const o2::gpu::TPCFastTransformPOD* mCorrMaps{nullptr};
  const o2::gpu::GPUParam* mTPCParam{nullptr};
  o2::tpc::VDriftCorrFact mVDrift{}; // drift calibration accounted by mCorrMaps
  Volume* mEnvelope{nullptr};        // fictitious top volume, owner of the drift calibration DOFs
  std::vector<Volume*> mSensors;     // sector volumes, indexed by sector
  const TPCMap* mClusterShift{nullptr}; // see setClusterShift, not owned
  // per-cluster bias correction of the cosmics (TPCMapParams.clusterBias*), loaded at the first cosmic
  void initClusterBias();
  bool mClusterBiasInit{false};
  int mClusterBiasPass{0};
  std::shared_ptr<const o2::tpc::cluster_bias::ClusterBiasCorrection> mClusterBias; // null: off
  std::array<float, NStacks> mClusterBiasMult{-1.f, -1.f, -1.f, -1.f};               // multStack per stack for the occupancy term (< 0: off)
  std::array<uint64_t, 3> mClusterBiasCount{};                                         // corrected, rejected (d_pred below the cut), no prediction
};

} // namespace o2::alignrs

#endif
