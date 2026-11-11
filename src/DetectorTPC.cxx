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
#include <cstdlib>
#include <format>
#include <memory>

#include <optional>
#include <sstream>
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include <algorithm>
#include <optional>
#include "DataFormatsGlobalTracking/TrackCosmicsExtended.h"
#include "DataFormatsTPC/ClusterNative.h"
#include "DataFormatsTPC/Constants.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "TPCBase/ParameterDetector.h"
#include "DataFormatsTPC/WorkflowHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/ClusterBiasCorrection.h"
#include "O2Align/DetectorTPC.h"
#include "O2Align/Params.h"
#include "O2Align/SensorTPC.h"
#include "O2Align/TPCMap.h"
#include "O2Align/TrackFit.h"
#include "TPCFastTransformPOD.h"

#include "GPUO2ExternalUser.h"
#include "GPUParam.h"
#include "GPUParam.inc"

namespace o2::alignrs
{

void DetectorTPC::prepareData(o2::globaltracking::RecoContainer* /*recoData*/)
{
  // Nothing can be precomputed per TF: the cluster transformation needs the track time and both
  // the supercluster building and the cluster errors need the track angles, hence everything is
  // done in prepareTrack.
}

const std::string& DetectorTPC::getTimeSlotsJson() const
{
  return Params::Instance().VDTimeSlotsJson;
}

void DetectorTPC::onSlotChange(int slotID)
{
  // the drift parameters of every slot are independent: the labels of the derivatives follow the slot
  mEnvelope->setActiveCalibSlot(slotID);
  LOGP(info, "TPC drift calibration slot {}", slotID);
}

void DetectorTPC::initCalib()
{
  auto* calib = mEnvelope ? dynamic_cast<TPCVDriftDOFSet*>(mEnvelope->getCalib()) : nullptr;
  if (!calib) {
    return; // the drift calibration was not requested by the DOF configuration
  }
  calib->setZLength(o2::tpc::ParameterDetector::Instance().TPClength);
}

Volume::Ptr DetectorTPC::buildHierarchy(Volume::SensorMapping& sensorMap)
{
  // The TPC has no alignable entries in the geometry: both the envelope and the sector volumes
  // are fictitious and define their frames themselves.
  uint32_t gLbl{0};
  const uint32_t det = mDetIdx;
  auto root = std::make_unique<EnvelopeTPC>("TPC_envelope", gLbl++, det, false, true);
  // w/o a geometry counterpart the envelope cannot be shifted, but it may own calibration DOFs
  root->setRigidBodyAllowed(false);
  // the drift is calibrated per time slot, all of them sharing the single DOF set of the envelope
  if (mTimeSlots) {
    root->setCalibSlots(*mTimeSlots);
  }
  mEnvelope = root.get();
  mSensors.assign(o2::tpc::constants::MAXSECTOR, nullptr);
  for (int isec = 0; isec < o2::tpc::constants::MAXSECTOR; ++isec) {
    const Label lbl(det, isec, true);
    const auto symName = std::format("TPC/sec{:02d}", isec);
    auto* sector = root->addChild<SensorTPC>(symName.c_str(), lbl, true);
    sector->setSensorId(isec);
    sensorMap[lbl] = sector;
    mSensors[isec] = sector;
  }
  return root;
}

namespace
{
/// TPC cluster accepted for the refit, transformed to the sector frame
struct SelectedCluster {
  uint8_t sector{0}, row{0};
  float qTot{0.f};
  int16_t state{0};
  float x{0.f}, y{0.f}, z{0.f};
  /// Y of the cluster in the lab frame, used to order the clusters of a cosmic along its flight
  float yLab() const
  {
    constexpr int NSectorsPerSide = o2::tpc::constants::MAXSECTOR / 2;
    const double alpha = o2::math_utils::detail::sector2Angle<double>(sector % NSectorsPerSide);
    return static_cast<float>(x * std::sin(alpha) + y * std::cos(alpha));
  }
};

/// the TPCMap cluster shift (DetectorTPC::setClusterShift) added to a transformed cluster, if any
void applyClusterShift(const TPCMap* shift, SelectedCluster& sel)
{
  if (shift) {
    const auto d = shift->clusterShift(sel.sector, sel.x, sel.y, sel.z);
    sel.x += d[0];
    sel.y += d[1];
    sel.z += d[2];
  }
}

/// charge-weighted merge of the clusters of a track on neighbouring pad-rows of the same sector
struct SuperCluster {
  explicit SuperCluster(const SelectedCluster& c) : sector(c.sector), firstRow(c.row) { add(c); }
  void add(const SelectedCluster& c)
  {
    sumX += c.qTot * c.x;
    sumY += c.qTot * c.y;
    sumZ += c.qTot * c.z;
    sumRow += c.qTot * c.row;
    charge += c.qTot;
    state |= c.state;
    ++nClusters;
  }
  std::array<float, 3> position() const { return {static_cast<float>(sumX / charge), static_cast<float>(sumY / charge), static_cast<float>(sumZ / charge)}; }
  uint8_t meanRow() const { return static_cast<uint8_t>(std::lround(sumRow / charge)); }
  double alpha() const { return o2::math_utils::detail::sector2Angle<double>(sector % (o2::tpc::constants::MAXSECTOR / 2)); }

  uint8_t sector{0}, firstRow{0};
  int16_t state{0};
  int nClusters{0};
  double sumX{0.}, sumY{0.}, sumZ{0.}, sumRow{0.}, charge{0.};
};

/// Merge the clusters provided by next() (std::optional<SelectedCluster>, ordered along the track) into the superclusters
/// of at most maxRowsCombined rows and hand them over to visit(const SuperCluster&), which returns false to stop.
/// npntCut, the number of points required, is relaxed by each merging. Returns the number of the accepted superclusters.
template <typename Next, typename Visit>
int forEachSuperCluster(Next&& next, int maxRowsCombined, int& npntCut, Visit&& visit)
{
  int npoints = 0;
  auto pending = next();
  while (pending) {
    SuperCluster sc(*pending);
    while ((pending = next()) && pending->sector == sc.sector && std::abs(pending->row - sc.firstRow) < maxRowsCombined) {
      sc.add(*pending);
    }
    npntCut -= sc.nClusters - 1;
    if (!visit(sc)) {
      break;
    }
    ++npoints;
  }
  return npoints;
}

/// errors (yy, yz, zz) of the supercluster for the track of the given angles
std::array<float, 3> superClusterCov(const o2::gpu::GPUParam& tpcParam, const SuperCluster& sc, float z, double snp, double tgl)
{
  const auto& params = Params::Instance();
  std::array<float, 3> cov{0.f, 0.f, 0.f};
  // TODO: this disables the occupancy / charge components of the error estimation
  tpcParam.GetClusterErrors2(sc.sector, sc.meanRow(), z, snp, tgl, -1.f, 0.f, 0.f, cov[0], cov[2]);
  tpcParam.UpdateClusterError2ByState(sc.state, cov[0], cov[2]);
  if (sc.nClusters > 1) { // conservative: the merged clusters are not independent measurements
    const float fact = 1.f / std::sqrt(static_cast<float>(sc.nClusters));
    cov[0] *= fact;
    cov[2] *= fact;
  }
  cov[0] += params.extraClsErrYTPC * params.extraClsErrYTPC;
  cov[2] += params.extraClsErrZTPC * params.extraClsErrZTPC;
  return cov;
}

FrameInfoExt superClusterFrame(int detIdx, const SuperCluster& sc, const std::array<float, 3>& pnt3D, const std::array<float, 3>& cov)
{
  FrameInfoExt pnt;
  pnt.detIdx = detIdx;
  pnt.lr = static_cast<int8_t>(DetectorTPC::getStack(sc.meanRow()));
  pnt.label = Label(detIdx, sc.sector, true);
  pnt.x = pnt3D[0];
  pnt.alpha = static_cast<float>(sc.alpha());
  pnt.cluster = o2::BaseCluster<float>(static_cast<int16_t>(sc.sector), pnt3D[0], pnt3D[1], pnt3D[2], cov[0], cov[2], cov[1]);
  return pnt;
}
} // namespace

/// cuts on the position of a cluster, besides the maximum pad-row which stops the outward refit
bool DetectorTPC::isClusterAccepted(const o2::tpc::ClusterNative& cl, int row) const
{
  const auto& params = Params::Instance();
  if (row < params.minTPCPadRow) {
    return false;
  }
  if (params.discardEdgePadrows > 0 && getDistanceToStackEdge(row) < params.discardEdgePadrows) {
    return false;
  }
  auto padFromEdge = cl.getPad();
  int npads = o2::gpu::GPUTPCGeometry::NPads(row);
  if (padFromEdge > npads / 2) {
    padFromEdge = npads - 1 - padFromEdge;
  }
  return padFromEdge >= params.discardEdgePadDepth;
}

bool DetectorTPC::prepareTrack(o2::globaltracking::RecoContainer* recoData, const GlobalIDSet& ids, Track& resTrack)
{
  const auto& params = Params::Instance();
  const auto gid = ids[GTrackID::TPC];
  if (!gid.isIndexSet()) {
    return false;
  }
  if (!mCorrMaps || !mTPCParam) {
    LOGP(fatal, "TPC correction maps ({}) and GPUParam ({}) must be set before processing TPC tracks",
         (void*)mCorrMaps, (void*)mTPCParam);
  }
  const auto& trk = recoData->getTrack<o2::tpc::TrackTPC>(gid);
  const int nClus = trk.getNClusters();
  if (nClus < params.minTPCClusters) {
    return false;
  }
  const auto clusterIdxStruct = recoData->getTPCTracksClusterRefs();
  const auto clusterNativeAccess = recoData->inputsTPCclusters->clusterIndex;
  const auto shMap = recoData->clusterShMapTPC;

  float trackTime{0.f}, trackTimeErr{0.f};
  recoData->getTrackTime(resTrack.gid, trackTime, trackTimeErr);
  const float tOffset = trackTime / (o2::constants::lhc::LHCBunchSpacingMUS * 8);

  auto prop = o2::base::PropagatorD::Instance();
  // we continue the fit of the seed prepared by the inner detectors outward, w/o resetting its cov.matrix
  auto trkParam = resTrack.track;
  o2::track::TrackParD trkRef, *refLin = nullptr;
  if (params.useStableRef) {
    refLin = &(trkRef = trkParam);
  }
  float chi2 = 0.f;
  const size_t nPointsIni = resTrack.info.size();

  // Accepted clusters in the outward direction: the clusters of the track are ordered from the
  // outermost to the innermost, hence they are traversed from the last one on.
  int iNext = nClus - 1;
  auto nextCluster = [&]() -> std::optional<SelectedCluster> {
    while (iNext >= 0) {
      uint8_t sector = 0, row = 0;
      const auto& cl = trk.getCluster(clusterIdxStruct, iNext--, clusterNativeAccess, sector, row);
      if (row > params.maxTPCPadRow) { // outward refit: all following clusters will have a larger padrow
        iNext = -1;
        break;
      }
      if (!isClusterAccepted(cl, row)) { // the following clusters still have a chance to be accepted
        continue;
      }
      SelectedCluster sel{.sector = sector, .row = row, .qTot = static_cast<float>(cl.getQtot()), .state = shMap[&cl - clusterNativeAccess.clustersLinear]};
      mCorrMaps->Transform(sector, row, cl.getPad(), cl.getTime(), sel.x, sel.y, sel.z, tOffset);
      applyClusterShift(mClusterShift, sel);
      return sel;
    }
    return std::nullopt;
  };

  // the min number of points is relaxed by each merging of clusters into a supercluster
  int npntCut = params.minTPCClusters;
  const int npoints = forEachSuperCluster(nextCluster, params.maxTPCRowsCombined, npntCut, [&](const SuperCluster& sc) {
    const auto pnt3D = sc.position();
    if (!prop->propagateToAlphaX(trkParam, refLin, sc.alpha(), pnt3D[0], false, params.maxSnp, params.maxStep, 1, params.corrType)) {
      return false;
    }
    const auto cov = superClusterCov(*mTPCParam, sc, pnt3D[2], trkParam.getSnp(), trkParam.getTgl());
    const std::array<double, 2> pos{pnt3D[1], pnt3D[2]};
    const std::array<double, 3> covD{cov[0], cov[1], cov[2]};
    chi2 += static_cast<float>(trkParam.getPredictedChi2Quiet(pos, covD));
    if (!trkParam.update(pos, covD)) {
      return false;
    }
    if (refLin) { // displace the reference to the last updated cluster
      refLin->setY(pos[0]);
      refLin->setZ(pos[1]);
    }
    resTrack.info.push_back(superClusterFrame(mDetIdx, sc, pnt3D, cov));
    return true;
  });
  if (params.verbose > 0) {
    LOGP(info, "Prepare {}: Nclus={} -> Npoints={}, chi2 = {}", gid.asString(), nClus, npoints, chi2);
  }
  if (npoints < npntCut) {
    resTrack.info.resize(nPointsIni);
    return false;
  }
  resTrack.track = trkParam; // the seed is replaced only by a successful fit
  resTrack.kfFit.chi2 += chi2;
  return true;
}

// The attached clusters of the leg, transformed with the time of the cosmic (the TOF one if there is any) common for the 2 legs,
// are ordered along the Y of the lab frame and merged to the superclusters. Their errors are evaluated for the angles of the reference track.
// The state of the shared clusters (not stored) is not accounted. Returns the number of the clusters used.
int DetectorTPC::prepareCosmicLeg(const o2::dataformats::TrackCosmicsExtended& cosm, int leg, const o2::track::TrackParCovD& ref, std::vector<FrameInfoExt>& frames)
{
  const auto& params = Params::Instance();
  if (!mCorrMaps || !mTPCParam) {
    LOGP(fatal, "TPC correction maps ({}) and GPUParam ({}) must be set before processing TPC tracks", (void*)mCorrMaps, (void*)mTPCParam);
  }
  const float timeMUS = cosm.hasTOFTime() ? cosm.timeTOFMUS : cosm.cosmic.getTimeMUS().getTimeStamp();
  const float tOffset = timeMUS / (o2::constants::lhc::LHCBunchSpacingMUS * 8);
  std::vector<SelectedCluster> selected;
  if (!mClusterBiasInit) {
    initClusterBias();
  }
  // the reco's per-cluster bias correction (TPCMapParams.clusterBias*), decided by the prediction = the reference track at the cluster's
  // row in its sector frame: the table's inward bias at the predicted pad phase and distance d_pred from the nearest pad-plane edge,
  // moved outward on the prediction's side; false = rejected (d_pred below the cut)
  const auto& mapParams = TPCMapParams::Instance();
  auto correctClusterBias = [&](const o2::tpc::ClusterNative& cl, SelectedCluster& sel) {
    namespace cb = o2::tpc::cluster_bias;
    constexpr int NSectorsPerSide = o2::tpc::constants::MAXSECTOR / 2;
    o2::track::TrackParCovD trk;
    if (!propagateCosmicRef(ref, o2::math_utils::detail::sector2Angle<double>(sel.sector % NSectorsPerSide), sel.x, trk)) {
      ++mClusterBiasCount[2]; // kept uncorrected: the supercluster's own propagation decides on it
      return true;
    }
    // the table is defined on the PAD PLANE (where the charge landed): the prediction there = the cluster's ideal (uncorrected) y +
    // (predicted - measured), the map being smooth over the residual; in PbPb the map moves clusters by mm, near the edges by ~1 pad
    float xIdeal = 0.f, yIdeal = 0.f, zIdeal = 0.f;
    mCorrMaps->TransformIdeal(sel.sector, sel.row, cl.getPad(), cl.getTime(), xIdeal, yIdeal, zIdeal, tOffset);
    const double padWidth = o2::gpu::GPUTPCGeometry::PadWidth(sel.row);
    const double yPred = yIdeal + (trk.getY() - sel.y);
    const double dPred = (0.5 * o2::gpu::GPUTPCGeometry::NPads(sel.row) * padWidth - std::abs(yPred)) / padWidth;
    if (dPred < mapParams.clusterBiasRejectPads) {
      ++mClusterBiasCount[1];
      return false;
    }
    const float pad = cl.getPad();
    const int cls = cb::classOfFlags(cl.getFlags(), std::abs(pad - std::round(pad)) < 0.25f / 64.f);
    const int stack = cb::stackOfRow(sel.row);
    const double snp = std::clamp(trk.getSnp(), -0.999, 0.999);
    const double absSlope = std::abs(snp / std::sqrt((1. - snp) * (1. + snp)));
    const double bias = mClusterBias->inwardBias(mClusterBias->byRegion ? cb::regionOfRow(sel.row) : stack, cls, cb::driftLength(sel.z), absSlope,
                                                 std::max(0., dPred), mClusterBiasMult[stack]);
    sel.y += static_cast<float>((yPred >= 0. ? 1. : -1.) * bias * mapParams.clusterBiasScale);
    ++mClusterBiasCount[0];
    return true;
  };
  for (const auto& cc : leg == 0 ? cosm.clTPCBottom : cosm.clTPCTop) {
    if (!cc.isAttached() || cc.row > params.maxTPCPadRow || !isClusterAccepted(cc.cl, cc.row)) {
      continue;
    }
    auto& sel = selected.emplace_back(SelectedCluster{.sector = cc.sector, .row = cc.row, .qTot = static_cast<float>(cc.cl.getQtot())});
    mCorrMaps->Transform(cc.sector, cc.row, cc.cl.getPad(), cc.cl.getTime(), sel.x, sel.y, sel.z, tOffset);
    if (params.cosmTPCMaxZOverR > 0.f && std::abs(sel.z) > params.cosmTPCMaxZOverR * std::hypot(sel.x, sel.y)) {
      selected.pop_back(); // outside of the acceptance of the tracks the correction map was made of
      continue;
    }
    applyClusterShift(mClusterShift, sel);
    if (mClusterBias) {
      bool keep = true;
      if (mClusterBiasPass == 1) { // the fit for the prediction: no correction, the clusters measured near the sector edge (pad plane) left out
        const double padWidth = o2::gpu::GPUTPCGeometry::PadWidth(sel.row);
        float xIdeal = 0.f, yIdeal = 0.f, zIdeal = 0.f;
        mCorrMaps->TransformIdeal(sel.sector, sel.row, cc.cl.getPad(), cc.cl.getTime(), xIdeal, yIdeal, zIdeal, tOffset);
        keep = (0.5 * o2::gpu::GPUTPCGeometry::NPads(sel.row) * padWidth - std::abs(yIdeal)) / padWidth >= ClusterBiasPassOneEdgePads;
      } else {
        keep = correctClusterBias(cc.cl, sel);
      }
      if (!keep) {
        selected.pop_back();
      }
    }
  }
  std::sort(selected.begin(), selected.end(), [](const auto& a, const auto& b) { return a.yLab() < b.yLab(); });

  size_t iNext = 0;
  auto nextCluster = [&]() -> std::optional<SelectedCluster> {
    return iNext < selected.size() ? std::optional<SelectedCluster>(selected[iNext++]) : std::nullopt;
  };
  int npntCut = 0, nClUsed = 0;
  forEachSuperCluster(nextCluster, params.maxTPCRowsCombined, npntCut, [&](const SuperCluster& sc) {
    const auto pnt3D = sc.position();
    o2::track::TrackParCovD trkParam;
    if (!propagateCosmicRef(ref, sc.alpha(), pnt3D[0], trkParam)) {
      return true; // the supercluster is just not used
    }
    frames.push_back(superClusterFrame(mDetIdx, sc, pnt3D, superClusterCov(*mTPCParam, sc, pnt3D[2], trkParam.getSnp(), trkParam.getTgl())));
    nClUsed += sc.nClusters;
    return true;
  });
  return nClUsed;
}

int DetectorTPC::getMinCosmicMeasurements() const
{
  return Params::Instance().minTPCClusters;
}

void DetectorTPC::initClusterBias()
{
  mClusterBiasInit = true;
  const auto& mapParams = TPCMapParams::Instance();
  if (mapParams.clusterBiasTable.empty()) {
    return;
  }
  auto table = std::make_shared<o2::tpc::cluster_bias::ClusterBiasCorrection>();
  std::string error;
  if (!table->load(mapParams.clusterBiasTable, &error)) {
    LOGP(fatal, "Cannot load the cluster bias table {}: {}", mapParams.clusterBiasTable, error);
  }
  if (!mapParams.clusterBiasMult.empty()) {
    std::stringstream ss(mapParams.clusterBiasMult);
    std::string item;
    for (int s = 0; s < NStacks; ++s) {
      if (!std::getline(ss, item, ',')) {
        LOGP(fatal, "TPCMapParams.clusterBiasMult needs {} values (IROC, OROC1-3): {}", NStacks, mapParams.clusterBiasMult);
      }
      mClusterBiasMult[s] = std::stof(item);
    }
  } else if (table->hasOccupancy()) {
    LOGP(warn, "Cluster bias table {} has an occupancy term (and any data correction folded into it), but no TPCMapParams.clusterBiasMult: left out",
         mapParams.clusterBiasTable);
  }
  LOGP(info, "Cosmics: per-cluster bias correction {} ({} records, occupancy term {}, mult {} {} {} {}), x {}, clusters with d_pred < {} pads rejected",
       mapParams.clusterBiasTable, table->byRegion ? "pad-region" : "stack", mClusterBiasMult[0] >= 0.f && table->hasOccupancy() ? "on" : "off",
       mClusterBiasMult[0], mClusterBiasMult[1], mClusterBiasMult[2], mClusterBiasMult[3], mapParams.clusterBiasScale, mapParams.clusterBiasRejectPads);
  mClusterBias = table;
}

std::string DetectorTPC::getClusterBiasStats() const
{
  return std::format("corrected {}, rejected (d_pred below the cut) {}, kept without a prediction {}", mClusterBiasCount[0], mClusterBiasCount[1],
                     mClusterBiasCount[2]);
}

} // namespace o2::alignrs
