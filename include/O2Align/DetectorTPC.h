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
#include <vector>

#include "DataFormatsTPC/VDriftCorrFact.h"
#include "O2Align/Detector.h"
#include "O2Align/Volume.h"

namespace o2::gpu
{
class GPUParam;
class TPCFastTransformPOD;
} // namespace o2::gpu

namespace o2::dataformats
{
struct TrackCosmicsExtended;
} // namespace o2::dataformats

namespace o2::alignrs
{

class DetectorTPC final : public Detector
{
 public:
  static constexpr int NStacks = 4; // IROC, OROC1, OROC2, OROC3

  DetectorTPC() : Detector(DetTPC) {}

  void prepareData(o2::globaltracking::RecoContainer* recoData) final;
  bool prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack) final;

  /// Seed resTrack and append the TPC points of both legs of a cosmic, ordered along the motion of the muon, fitting them
  /// with the KF on the way. Each cluster is measured in its sector frame or, where the muon moves inwards there, in the sector
  /// frame rotated by pi, so the motion always points towards +x as the track parameterisation requires. All the clusters are
  /// transformed with the time of the cosmic (timeMUS, since the start of the TF), so the two legs share one z reference. In
  /// gaps (between the legs) material-only points are added. ptB0 > 0 is imposed as the pT of the seed (scattering errors w/o
  /// field). On success resTrack.track is the state at the last point; on failure resTrack is left in an undefined state and the
  /// reason is returned.
  enum class CosmicStatus : uint8_t { OK,
                                      TooFewClusters,
                                      PropagationFailed,
                                      UpdateFailed };
  CosmicStatus prepareCosmicExt(const o2::dataformats::TrackCosmicsExtended& cosm, float timeMUS, float ptB0, Track& resTrack);
  /// clusters of the cosmics skipped by prepareCosmicExt (propagation failed or |snp| above Params::cosmicsMaxClusterSnp)
  uint64_t getCosmicClustersSkipped() const { return mCosmicClustersSkipped; }

  /// both must be provided by the caller before the 1st prepareTrack call
  void setCorrMaps(const o2::gpu::TPCFastTransformPOD* maps) { mCorrMaps = maps; }
  void setTPCParam(const o2::gpu::GPUParam* par) { mTPCParam = par; }
  const o2::gpu::TPCFastTransformPOD* getCorrMaps() const { return mCorrMaps; }
  const o2::gpu::GPUParam* getTPCParam() const { return mTPCParam; }

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
  /// padrow / pad cuts of Params (pad-row range, stack edges, sector edges)
  static bool acceptClusterPosition(int row, float pad);

  /// In a gap between clusters of a cosmic: propagate trk towards the next cluster (global position xGloTo, yGloTo), appending
  /// material-only points every Params::cosmicsMaterialStep along the direction of the track, in the frame rotated to that
  /// direction (snp = 0 there). The material near the beam line and of the inner field cages is then lumped where it is, not
  /// all at the next cluster. False if the propagation failed.
  bool addCosmicMaterialPoints(double xGloTo, double yGloTo, o2::track::TrackParCovD& trk, o2::track::TrackParD* refLin, Track& resTrack) const;

  static constexpr std::array<std::pair<int, int>, NStacks> sStackMinMaxRow = {
    std::pair<int, int>{0, 62}, std::pair<int, int>{63, 96}, std::pair<int, int>{97, 126}, std::pair<int, int>{127, 151}};

  const o2::gpu::TPCFastTransformPOD* mCorrMaps{nullptr};
  const o2::gpu::GPUParam* mTPCParam{nullptr};
  o2::tpc::VDriftCorrFact mVDrift{}; // drift calibration accounted by mCorrMaps
  Volume* mEnvelope{nullptr};    // fictitious top volume, owner of the drift calibration DOFs
  std::vector<Volume*> mSensors; // sector volumes, indexed by sector
  uint64_t mCosmicClustersSkipped{0}; // see getCosmicClustersSkipped
};

} // namespace o2::alignrs

#endif
