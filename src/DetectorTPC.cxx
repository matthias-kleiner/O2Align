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
#include <unordered_map>
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
#include "O2Align/ResidualReader.h"
#include "O2Align/DetectorTPC.h"
#include "O2Align/Params.h"
#include "O2Align/SensorTPC.h"
#include "O2Align/TPCClusterCorrections.h"
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
  bool road{false}; // cosmics: found in the road, not attached to the leg
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
    const float w = c.qTot > 0.f ? c.qTot : 1.f; // clusters stored w/o the charge: equal weights instead of a NaN position
    sumX += w * c.x;
    sumY += w * c.y;
    sumZ += w * c.z;
    sumRow += w * c.row;
    charge += w;
    state |= c.state;
    road |= c.road;
    ++nClusters;
  }
  std::array<float, 3> position() const { return {static_cast<float>(sumX / charge), static_cast<float>(sumY / charge), static_cast<float>(sumZ / charge)}; }
  uint8_t meanRow() const { return static_cast<uint8_t>(std::lround(sumRow / charge)); }
  double alpha() const { return o2::math_utils::detail::sector2Angle<double>(sector % (o2::tpc::constants::MAXSECTOR / 2)); }

  uint8_t sector{0}, firstRow{0};
  int16_t state{0};
  bool road{false};
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
  if (sc.road) {
    pnt.setTPCRoad();
  }
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
    if (!prop->propagateToAlphaX(trkParam, refLin, sc.alpha(), pnt3D[0], false, params.maxSnp, params.maxStep, 1, static_cast<o2::base::PropagatorD::MatCorrType>(params.corrType))) {
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

// The TPC clusters of a residual track (ResidualReader): IDEAL positions (the residual file's reference drift), corrected here with the
// static map set as mCorrMaps + the derivative + the sector-edge term (correctCluster), with the tables by the track in the lookup pass
// (setResidualReference: residualTables), shifted by the TPCMap cluster shift if any, merged to superclusters and fitted outward from the
// seed prepared by the inner detectors, as prepareTrack. Rows outside [minTPCPadRow, maxTPCPadRow], the stack-edge rows
// (discardEdgePadrows) and the clusters flagged by validateTrack are not used, nor those measured within max(discardEdgePadDepth,
// residualTPCEdgePads) pads from the sector edge or tpcCentrePads from the centre line (pad-plane y) unless the lookup pass handles them
// by the track (clusterBiasTable / centreTable).
bool DetectorTPC::prepareResidualTrack(const ResidualTrack& trk, Track& resTrack)
{
  const auto& params = Params::Instance();
  if (!mCorrMaps || !mTPCParam) {
    LOGP(fatal, "TPC correction maps ({}) and GPUParam ({}) must be set before processing TPC tracks", (void*)mCorrMaps, (void*)mTPCParam);
  }
  if (static_cast<int>(trk.tpc.size()) < params.minTPCClusters) {
    return false;
  }
  std::vector<const ResidualTPCCluster*> clusters;
  clusters.reserve(trk.tpc.size());
  for (const auto& c : trk.tpc) {
    clusters.push_back(&c);
  }
  std::stable_sort(clusters.begin(), clusters.end(), [](const auto* a, const auto* b) { return a->row < b->row; }); // outward
  const bool lookup = mResidualRef != nullptr && hasTrackTables();
  const bool edgeByTrack = lookup && mClusterBias != nullptr;
  const bool centreByTrack = lookup && mCentreTable != nullptr;
  const float edgeCut = std::max(static_cast<float>(params.discardEdgePadDepth), params.residualTPCEdgePads);
  if (mDriftNorm.on && !(mResRefVDrift > 0.f)) {
    LOGP(fatal, "TPCMapParams.driftNorm needs the residual file's drift reference (none set)");
  }
  o2::track::TrackParCovD refAt; // the reference propagated outward from cluster to cluster (the clusters are ordered by row)
  bool refOK = lookup;
  if (lookup) {
    refAt = *mResidualRef;
  }
  // a cut on the measured pad-plane y: the edge zone and the centre zone unless the lookup handles them
  auto measuredCut = [&](const ResidualTPCCluster& c, bool byTrack) {
    const double padWidth = o2::gpu::GPUTPCGeometry::PadWidth(c.row);
    if (!(byTrack && edgeByTrack) && edgeCut > 0.f && (0.5 * o2::gpu::GPUTPCGeometry::NPads(c.row) * padWidth - std::abs(c.y)) / padWidth < edgeCut) {
      return true;
    }
    return !(byTrack && centreByTrack) && params.tpcCentrePads > 0.f && std::abs(c.y) < params.tpcCentrePads * padWidth;
  };
  auto prop = o2::base::PropagatorD::Instance();
  auto trkParam = resTrack.track;
  o2::track::TrackParD trkRef, *refLin = nullptr;
  if (params.useStableRef) {
    refLin = &(trkRef = trkParam);
  }
  float chi2 = 0.f;
  const size_t nPointsIni = resTrack.info.size();
  size_t iNext = 0;
  auto nextCluster = [&]() -> std::optional<SelectedCluster> {
    while (iNext < clusters.size()) {
      const auto& c = *clusters[iNext++];
      if (c.row > params.maxTPCPadRow) {
        iNext = clusters.size();
        break;
      }
      if (c.rejected || c.row < params.minTPCPadRow || (params.discardEdgePadrows > 0 && getDistanceToStackEdge(c.row) < params.discardEdgePadrows)) {
        continue;
      }
      if (measuredCut(c, true)) {
        continue;
      }
      SelectedCluster sel{.sector = static_cast<uint8_t>(c.sector), .row = static_cast<uint8_t>(c.row)};
      sel.x = c.x;
      sel.y = c.y;
      sel.z = mDriftNorm.on ? driftNormZ(mDriftNorm, c.z, c.sector >= 18, mResRefVDrift, mResRefT0) : c.z;
      const float zNominal = sel.z;
      correctCluster(c.sector, c.row, c.x, sel.x, sel.y, sel.z);
      if (lookup) {
        if (refOK) {
          refOK = propagateCosmicRef(refAt, o2::math_utils::detail::sector2Angle<double>(c.sector % (o2::tpc::constants::MAXSECTOR / 2)), sel.x, refAt);
        }
        if (!refOK) { // no prediction: the measured cuts, no table
          ++mClusterBiasCount[2];
          if (measuredCut(c, false)) {
            continue;
          }
        } else {
          float dy = 0.f;
          if (!residualTables(trk, c, sel.y, zNominal, refAt, dy)) {
            continue;
          }
          sel.y += dy;
        }
      }
      applyClusterShift(mClusterShift, sel);
      return sel;
    }
    return std::nullopt;
  };
  int npntCut = params.minTPCClusters;
  const int npoints = forEachSuperCluster(nextCluster, params.maxTPCRowsCombined, npntCut, [&](const SuperCluster& sc) {
    const auto pnt3D = sc.position();
    if (!prop->propagateToAlphaX(trkParam, refLin, sc.alpha(), pnt3D[0], false, params.maxSnp, params.maxStep, 1, static_cast<o2::base::PropagatorD::MatCorrType>(params.corrType))) {
      return false;
    }
    const auto cov = superClusterCov(*mTPCParam, sc, pnt3D[2], trkParam.getSnp(), trkParam.getTgl());
    const std::array<double, 2> pos{pnt3D[1], pnt3D[2]};
    const std::array<double, 3> covD{cov[0], cov[1], cov[2]};
    chi2 += static_cast<float>(trkParam.getPredictedChi2Quiet(pos, covD));
    if (!trkParam.update(pos, covD)) {
      return false;
    }
    if (refLin) {
      refLin->setY(pos[0]);
      refLin->setZ(pos[1]);
    }
    resTrack.info.push_back(superClusterFrame(mDetIdx, sc, pnt3D, cov));
    return true;
  });
  if (npoints < npntCut) {
    resTrack.info.resize(nPointsIni);
    return false;
  }
  resTrack.track = trkParam;
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
      ++mCentreTableCount[1];
      return true;
    }
    // the table is defined on the PAD PLANE (where the charge landed): the prediction there = the cluster's ideal (uncorrected) y +
    // (predicted - measured), the map being smooth over the residual; in PbPb the map moves clusters by mm, near the edges by ~1 pad
    float xIdeal = 0.f, yIdeal = 0.f, zIdeal = 0.f;
    mCorrMaps->TransformIdeal(sel.sector, sel.row, cl.getPad(), cl.getTime(), xIdeal, yIdeal, zIdeal, tOffset);
    const double padWidth = o2::gpu::GPUTPCGeometry::PadWidth(sel.row);
    const double yPred = yIdeal + (trk.getY() - sel.y);
    const double dPred = (0.5 * o2::gpu::GPUTPCGeometry::NPads(sel.row) * padWidth - std::abs(yPred)) / padWidth;
    if (mClusterBias && dPred < mapParams.clusterBiasRejectPads) {
      ++mClusterBiasCount[1];
      return false;
    }
    const float pad = cl.getPad();
    const int cls = cb::classOfFlags(cl.getFlags(), std::abs(pad - std::round(pad)) < 0.25f / 64.f);
    if (mCentreTable) { // the spacer-cross table at the prediction (TPCMapParams.centreTable)
      const double shift = centreShift(sel.sector, sel.row, cls, yPred);
      if (shift != 0.) {
        sel.y -= static_cast<float>(shift);
        ++mCentreTableCount[0];
      }
    }
    if (!mClusterBias) {
      return true;
    }
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
    const bool road = !cc.isAttached();
    if (road && !(params.cosmTPCUseRoad > 0 && cc.isCorridor() && (params.cosmTPCUseRoad > 1 || !cc.isUsed()))) {
      continue;
    }
    if (cc.row > params.maxTPCPadRow || !isClusterAccepted(cc.cl, cc.row)) {
      continue;
    }
    if (params.tpcCentrePads > 0.f && !(mCentreTable && mClusterBiasPass != 1) &&
        std::abs(cc.cl.getPad() - 0.5f * (o2::gpu::GPUTPCGeometry::NPads(cc.row) - 1)) < params.tpcCentrePads) {
      continue; // the sector centre line (spacer cross) out, unless the centre table corrects it (pass 2)
    }
    auto& sel = selected.emplace_back(SelectedCluster{.sector = cc.sector, .row = cc.row, .qTot = static_cast<float>(cc.cl.getQtot())});
    sel.road = road;
    if (mDriftNorm.on || mDerivative || mEdgeTerm) { // the ideal transformation; on the maps' drift frame: with vD0 / t00 (set by the
                                                     // caller) and the window's vD slope; then the map + derivative + sector-edge term
      mCorrMaps->TransformIdeal(cc.sector, cc.row, cc.cl.getPad(), cc.cl.getTime(), sel.x, sel.y, sel.z, tOffset);
      if (mDriftNorm.on) {
        sel.z = driftNormZ(mDriftNorm, sel.z, cc.sector >= 18, mDriftNorm.vD0, mDriftNorm.t00);
      }
      correctCluster(cc.sector, cc.row, sel.x, sel.x, sel.y, sel.z);
    } else {
      mCorrMaps->Transform(cc.sector, cc.row, cc.cl.getPad(), cc.cl.getTime(), sel.x, sel.y, sel.z, tOffset);
    }
    if (params.cosmTPCMaxZOverR > 0.f && std::abs(sel.z) > params.cosmTPCMaxZOverR * std::hypot(sel.x, sel.y)) {
      selected.pop_back(); // outside of the acceptance of the tracks the correction map was made of
      continue;
    }
    applyClusterShift(mClusterShift, sel);
    if (mClusterBias || mCentreTable) {
      bool keep = true;
      if (mClusterBiasPass == 1) { // the fit for the prediction: no correction, the clusters measured near the sector edge (pad plane) left out
        if (mClusterBias) {
          const double padWidth = o2::gpu::GPUTPCGeometry::PadWidth(sel.row);
          float xIdeal = 0.f, yIdeal = 0.f, zIdeal = 0.f;
          mCorrMaps->TransformIdeal(sel.sector, sel.row, cc.cl.getPad(), cc.cl.getTime(), xIdeal, yIdeal, zIdeal, tOffset);
          keep = (0.5 * o2::gpu::GPUTPCGeometry::NPads(sel.row) * padWidth - std::abs(yIdeal)) / padWidth >= ClusterBiasPassOneEdgePads;
        }
      } else {
        keep = correctClusterBias(cc.cl, sel);
      }
      if (!keep) {
        selected.pop_back();
      }
    }
  }
  if (params.cosmTPCUseRoad > 0) {
    // road clusters: only on the pad rows w/o an attached cluster, per row the one closest to the reference track, inside the window
    constexpr int NSectorsPerSide = o2::tpc::constants::MAXSECTOR / 2;
    const auto rowKey = [](const SelectedCluster& c) { return static_cast<int>(c.sector) * 256 + c.row; };
    std::unordered_map<int, std::pair<int, float>> best; // row -> (index, distance^2); index -1: the row has an attached cluster
    for (int i = 0; i < static_cast<int>(selected.size()); ++i) {
      if (!selected[i].road) {
        best[rowKey(selected[i])] = {-1, 0.f};
      }
    }
    for (int i = 0; i < static_cast<int>(selected.size()); ++i) {
      const auto& sel = selected[i];
      if (!sel.road) {
        continue;
      }
      const auto it = best.find(rowKey(sel));
      if (it != best.end() && it->second.first < 0) {
        continue;
      }
      o2::track::TrackParCovD trk;
      if (!propagateCosmicRef(ref, o2::math_utils::detail::sector2Angle<double>(sel.sector % NSectorsPerSide), sel.x, trk)) {
        continue;
      }
      const float dy = sel.y - static_cast<float>(trk.getY()), dz = sel.z - static_cast<float>(trk.getZ());
      if (std::abs(dy) > params.cosmTPCRoadMaxDY || std::abs(dz) > params.cosmTPCRoadMaxDZ) {
        continue;
      }
      const float d2 = dy * dy + dz * dz;
      if (it == best.end() || d2 < it->second.second) {
        best[rowKey(sel)] = {i, d2};
      }
    }
    std::vector<SelectedCluster> kept;
    kept.reserve(selected.size());
    for (int i = 0; i < static_cast<int>(selected.size()); ++i) {
      if (!selected[i].road) {
        kept.push_back(selected[i]);
      } else if (const auto it = best.find(rowKey(selected[i])); it != best.end() && it->second.first == i) {
        kept.push_back(selected[i]);
      }
    }
    selected.swap(kept);
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
  if (params.verbose > 1) {
    const auto& all = leg == 0 ? cosm.clTPCBottom : cosm.clTPCTop;
    LOGP(info, "Cosmic TPC leg {}: {} clusters in the input, {} selected, {} used; time {} mus, 1st selected: sector {} row {} x {} y {} z {} qTot {}; ref: {}", leg, all.size(),
         selected.size(), nClUsed, timeMUS, selected.empty() ? -1 : selected[0].sector, selected.empty() ? -1 : selected[0].row, selected.empty() ? 0.f : selected[0].x,
         selected.empty() ? 0.f : selected[0].y, selected.empty() ? 0.f : selected[0].z, selected.empty() ? 0.f : selected[0].qTot, ref.asString());
  }
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
    LOGP(info, "Cluster bias table {} has an occupancy term: cosmics leave it out (no TPCMapParams.clusterBiasMult), residual-ext takes the track's multStack",
         mapParams.clusterBiasTable);
  }
  LOGP(info, "TPC per-cluster bias correction {} ({} records, occupancy term with clusterBiasMult {}, mult {} {} {} {}), x {}, clusters with d_pred < {} pads rejected",
       mapParams.clusterBiasTable, table->byRegion ? "pad-region" : "stack", mClusterBiasMult[0] >= 0.f && table->hasOccupancy() ? "on" : "off",
       mClusterBiasMult[0], mClusterBiasMult[1], mClusterBiasMult[2], mClusterBiasMult[3], mapParams.clusterBiasScale, mapParams.clusterBiasRejectPads);
  mClusterBias = table;
}

void DetectorTPC::correctCluster(int sector, int row, float rowX, float& x, float& y, float& z) const
{
  const float yNominal = y, zNominal = z;
  mCorrMaps->TransformXYZ(sector, row, x, y, z);
  if (mDerivative && mDerivativeScale != 0.f) { // at the nominal position, as the map (production: map + lumi scale x derivative)
    float dx = 0.f, dy = 0.f, dz = 0.f;
    mDerivative->getCorrectionLocal(sector, row, yNominal, zNominal, dx, dy, dz);
    x += mDerivativeScale * dx;
    y += mDerivativeScale * dy;
    z += mDerivativeScale * dz;
  }
  if (mEdgeTerm) { // on the corrected cluster
    mEdgeTerm->apply(sector, row, rowX, x, y, z, mEdgeScales);
  }
}

double DetectorTPC::centreShift(int sector, int row, int cls, double trackPadY) const
{
  const double padWidth = o2::gpu::GPUTPCGeometry::PadWidth(row);
  const double shift = mCentreTable->shift(sector, trackPadY < 0. ? 0 : 1, row, cls, std::abs(trackPadY) / padWidth);
  return (trackPadY < 0. ? -1. : 1.) * shift;
}

// As the map chain's DCA validation (gblfit collisionDcaSeries): the track's position on the pad plane = the measured pad-plane y minus the
// residual (corrected cluster - track); its distance d from the pad-plane edge on the cluster's side decides the drop (d below
// clusterBiasRejectPads) and the edge table's bias (moved outward on the cluster's side, occupancy = the track's multStack of the stack);
// the centre table is looked up at |track y| with the cluster class from its flags and pad phase.
bool DetectorTPC::residualTables(const ResidualTrack& trk, const ResidualTPCCluster& c, float yCorr, float zNominal, const o2::track::TrackParCovD& trkAt, float& dy)
{
  namespace cb = o2::tpc::cluster_bias;
  const auto& mapParams = TPCMapParams::Instance();
  const double padWidth = o2::gpu::GPUTPCGeometry::PadWidth(c.row);
  const double halfWidth = 0.5 * o2::gpu::GPUTPCGeometry::NPads(c.row) * padWidth;
  const double trackPadY = c.y - (yCorr - trkAt.getY());
  const double side = c.y < 0.f ? -1. : 1.;
  const double d = (halfWidth - side * trackPadY) / padWidth;
  dy = 0.f;
  if (mClusterBias && d < mapParams.clusterBiasRejectPads) {
    ++mClusterBiasCount[1];
    return false;
  }
  const double padMeasured = (halfWidth - std::abs(c.y)) / padWidth - 0.5;
  const int cls = cb::classOfFlags(c.flags, std::abs(padMeasured - std::round(padMeasured)) < 0.25 / 64.);
  if (mCentreTable) {
    const double shift = centreShift(c.sector, c.row, cls, trackPadY);
    if (shift != 0.) {
      dy -= static_cast<float>(shift);
      ++mCentreTableCount[0];
    }
  }
  if (mClusterBias) {
    const double snp = std::clamp(trkAt.getSnp(), -0.999, 0.999);
    const double absSlope = std::abs(snp / std::sqrt((1. - snp) * (1. + snp)));
    const int stack = cb::stackOfRow(c.row);
    const double mult = mClusterBias->hasOccupancy() ? static_cast<double>(trk.multStack[stack]) : -1.;
    const double bias = mClusterBias->inwardBias(mClusterBias->unitOfRow(c.row), cls, cb::driftLength(zNominal), absSlope, d, mult);
    if (bias != 0.) {
      dy += static_cast<float>(side * bias * mapParams.clusterBiasScale);
      ++mClusterBiasCount[0];
    }
  }
  return true;
}

std::string DetectorTPC::getCentreTableStats() const
{
  return std::format("shifted {}, kept without a prediction {}", mCentreTableCount[0], mCentreTableCount[1]);
}

std::string DetectorTPC::getClusterBiasStats() const
{
  return std::format("corrected {}, rejected (d_pred below the cut) {}, kept without a prediction {}", mClusterBiasCount[0], mClusterBiasCount[1],
                     mClusterBiasCount[2]);
}

} // namespace o2::alignrs
