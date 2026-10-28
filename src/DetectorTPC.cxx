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
#include "CommonConstants/LHCConstants.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
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
#include "O2Align/DetectorTPC.h"
#include "O2Align/Params.h"
#include "O2Align/SensorTPC.h"
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
};

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

  uint8_t sector{0}, firstRow{0};
  int16_t state{0};
  int nClusters{0};
  double sumX{0.}, sumY{0.}, sumZ{0.}, sumRow{0.}, charge{0.};
};
} // namespace

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

  constexpr float TAN10 = 0.17632698f; // tan of the half sector opening angle
  constexpr int NSectorsPerSide = o2::tpc::constants::MAXSECTOR / 2;
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
      if (row < params.minTPCPadRow) { // the following clusters still have a chance to be accepted
        continue;
      }
      if (params.discardEdgePadrows > 0 && getDistanceToStackEdge(row) < params.discardEdgePadrows) {
        continue;
      }
      auto padFromEdge = cl.getPad();
      int npads = o2::gpu::GPUTPCGeometry::NPads(row);
      if (padFromEdge > npads / 2) {
        padFromEdge = npads - 1 - padFromEdge;
      }
      if (padFromEdge < params.discardEdgePadDepth) {
        continue;
      }
      SelectedCluster sel{.sector = sector, .row = row, .qTot = static_cast<float>(cl.getQtot()), .state = shMap[&cl - clusterNativeAccess.clustersLinear]};
      mCorrMaps->Transform(sector, row, cl.getPad(), cl.getTime(), sel.x, sel.y, sel.z, tOffset);
      return sel;
    }
    return std::nullopt;
  };

  // the min number of points is relaxed by each merging of clusters into a supercluster
  int npntCut = params.minTPCClusters;
  int npoints = 0;
  auto pending = nextCluster();
  while (pending) {
    // merge the following clusters of the same sector within maxTPCRowsCombined rows from the 1st one
    SuperCluster sc(*pending);
    while ((pending = nextCluster()) && pending->sector == sc.sector && std::abs(pending->row - sc.firstRow) < params.maxTPCRowsCombined) {
      sc.add(*pending);
    }
    npntCut -= sc.nClusters - 1;
    const auto pnt3D = sc.position();
    const uint8_t meanRow = sc.meanRow();

    const double alpha = o2::math_utils::detail::sector2Angle<double>(sc.sector % NSectorsPerSide);
    if (!prop->propagateToAlphaX(trkParam, refLin, alpha, pnt3D[0], false, params.maxSnp, params.maxStep, 1, params.corrType)) {
      break;
    }
    std::array<float, 3> cov{0.f, 0.f, 0.f};
    // TODO: this disables the occupancy / charge components of the error estimation
    mTPCParam->GetClusterErrors2(sc.sector, meanRow, pnt3D[2], trkParam.getSnp(), trkParam.getTgl(), -1.f, 0.f, 0.f, cov[0], cov[2]);
    mTPCParam->UpdateClusterError2ByState(sc.state, cov[0], cov[2]);
    if (sc.nClusters > 1) { // conservative: the merged clusters are not independent measurements
      const float fact = 1.f / std::sqrt(static_cast<float>(sc.nClusters));
      cov[0] *= fact;
      cov[2] *= fact;
    }
    cov[0] += params.extraClsErrYTPC * params.extraClsErrYTPC;
    cov[2] += params.extraClsErrZTPC * params.extraClsErrZTPC;

    const std::array<double, 2> pos{pnt3D[1], pnt3D[2]};
    const std::array<double, 3> covD{cov[0], cov[1], cov[2]};
    chi2 += static_cast<float>(trkParam.getPredictedChi2Quiet(pos, covD));
    if (!trkParam.update(pos, covD)) {
      break;
    }
    if (refLin) { // displace the reference to the last updated cluster
      refLin->setY(pos[0]);
      refLin->setZ(pos[1]);
    }

    auto& pnt = resTrack.info.emplace_back();
    pnt.detIdx = mDetIdx;
    pnt.lr = static_cast<int8_t>(getStack(meanRow));
    pnt.label = Label(mDetIdx, sc.sector, true);
    pnt.x = pnt3D[0];
    pnt.alpha = static_cast<float>(alpha);
    pnt.cluster = o2::BaseCluster<float>(static_cast<int16_t>(sc.sector), pnt3D[0], pnt3D[1], pnt3D[2], cov[0], cov[2], cov[1]);
    ++npoints;
  }

  if (npoints < npntCut) {
    resTrack.info.resize(nPointsIni);
    return false;
  }
  resTrack.track = trkParam; // the seed is replaced only by a successful fit
  resTrack.kfFit.chi2 += chi2;
  return true;
}

bool DetectorTPC::acceptClusterPosition(int row, float pad)
{
  const auto& params = Params::Instance();
  if (row < params.minTPCPadRow || row > params.maxTPCPadRow) {
    return false;
  }
  if (params.discardEdgePadrows > 0 && getDistanceToStackEdge(row) < params.discardEdgePadrows) {
    return false;
  }
  auto padFromEdge = pad;
  const int npads = o2::gpu::GPUTPCGeometry::NPads(row);
  if (padFromEdge > npads / 2) {
    padFromEdge = npads - 1 - padFromEdge;
  }
  return padFromEdge >= params.discardEdgePadDepth;
}

bool DetectorTPC::addCosmicMaterialPoints(double xGloTo, double yGloTo, o2::track::TrackParCovD& trk, o2::track::TrackParD* refLin, Track& resTrack) const
{
  const auto& params = Params::Instance();
  auto prop = o2::base::PropagatorD::Instance();
  // the frame rotated to the direction of the track: snp = 0 there, the track moves towards +x
  double alphaDir = trk.getPhi();
  o2::math_utils::detail::bringToPMPi<double>(alphaDir);
  const double cosDir = std::cos(alphaDir);
  const double sinDir = std::sin(alphaDir);
  std::array<double, 3> xyz{};
  trk.getXYZGlo(xyz);
  const double xFrom = xyz[0] * cosDir + xyz[1] * sinDir;
  const double xEnd = xGloTo * cosDir + yGloTo * sinDir;
  const double step = params.cosmicsMaterialStep;
  for (double x = xFrom + step; step > 0. && x < xEnd - 0.5 * step; x += step) {
    if (!prop->propagateToAlphaX(trk, refLin, alphaDir, x, false, params.maxSnp, params.maxStep, 1, params.corrType)) {
      return false;
    }
    auto& pnt = resTrack.info.emplace_back();
    pnt.detIdx = mDetIdx;
    pnt.lr = 0;
    pnt.setMaterialOnly();
    pnt.x = static_cast<float>(x);
    pnt.alpha = static_cast<float>(alphaDir);
    pnt.cluster = o2::BaseCluster<float>(-1, x, trk.getY(), trk.getZ(), 0.f, 0.f, 0.f);
  }
  return true;
}

DetectorTPC::CosmicStatus DetectorTPC::prepareCosmicExt(const o2::dataformats::TrackCosmicsExtended& cosm, float timeMUS, float ptB0, Track& resTrack)
{
  const auto& params = Params::Instance();
  if (!mCorrMaps || !mTPCParam) {
    LOGP(fatal, "TPC correction maps ({}) and GPUParam ({}) must be set before processing TPC cosmics",
         (void*)mCorrMaps, (void*)mTPCParam);
  }
  constexpr int NSectorsPerSide = o2::tpc::constants::MAXSECTOR / 2;
  const float tOffset = timeMUS / (o2::constants::lhc::LHCBunchSpacingMUS * 8);

  // seed at the start of the motion: the outer parameters of the upper leg inverted to the direction of the motion, or the inner
  // ones of the lower leg if the upper one has no cluster
  o2::track::TrackParCov seed;
  if (!cosm.clTPCTop.empty()) {
    seed = cosm.tpcTop.getOuterParam();
    seed.invert();
  } else {
    seed = cosm.tpcBottom;
  }
  const double phiDir = seed.getPhi();
  const double cosLam = 1. / std::sqrt(1. + seed.getTgl() * seed.getTgl());
  const std::array<double, 3> dir{std::cos(phiDir) * cosLam, std::sin(phiDir) * cosLam, seed.getTgl() * cosLam};

  // the clusters of both legs ordered along the direction of the motion; each is measured in its sector frame or, where the muon
  // moves inwards there, in the sector frame rotated by pi, so that the motion always points to +x as the parameterisation requires
  struct OrderedCluster {
    SelectedCluster cl;
    double xGlo{0.};
    double yGlo{0.};
    double s{0.};        // position along the direction of the motion
    bool rotated{false}; // measured in the sector frame rotated by pi: x, y -> -x, -y
  };
  std::vector<OrderedCluster> clusters;
  for (const auto* leg : {&cosm.clTPCTop, &cosm.clTPCBottom}) {
    for (const auto& c : *leg) {
      if (!c.isAttached() || !acceptClusterPosition(c.row, c.cl.getPad())) { // the corridor clusters are not used
        continue;
      }
      // no shared-cluster map for the cosmics: the state is given by the cluster flags (split, edge, single) alone
      OrderedCluster oc{.cl = {.sector = c.sector, .row = c.row, .qTot = static_cast<float>(c.cl.getQtot()), .state = static_cast<int16_t>(c.cl.getFlags())}};
      mCorrMaps->Transform(c.sector, c.row, c.cl.getPad(), c.cl.getTime(), oc.cl.x, oc.cl.y, oc.cl.z, tOffset);
      if (params.verbose > 2) {
        LOGP(info, "COSMCL {} {} {:.4f} {:.4f} {:.5f} {:.5f} {:.5f}", c.sector, c.row, c.cl.getPad(), c.cl.getTime(), oc.cl.x, oc.cl.y, oc.cl.z);
      }
      const double alpha = o2::math_utils::detail::sector2Angle<double>(c.sector % NSectorsPerSide);
      const double cosA = std::cos(alpha);
      const double sinA = std::sin(alpha);
      oc.xGlo = oc.cl.x * cosA - oc.cl.y * sinA;
      oc.yGlo = oc.cl.x * sinA + oc.cl.y * cosA;
      oc.s = oc.xGlo * dir[0] + oc.yGlo * dir[1] + oc.cl.z * dir[2];
      oc.rotated = cosA * dir[0] + sinA * dir[1] < 0.;
      clusters.push_back(oc);
    }
  }
  if (clusters.size() < static_cast<size_t>(params.minTPCClusters)) {
    return CosmicStatus::TooFewClusters;
  }
  std::stable_sort(clusters.begin(), clusters.end(), [](const auto& a, const auto& b) { return a.s < b.s; });

  resTrack.track = convertTrack<double>(seed);
  resTrack.track.setPID(o2::track::PID::Muon); // energy loss and scattering of a muon, not of the PID guess of the TPC track
  if (ptB0 > 0.f) {
    resTrack.track.setQ2Pt(1. / ptB0);
  }
  resetTrackCovariance(resTrack.track);
  resTrack.kfFit.chi2 = 0.f;

  auto prop = o2::base::PropagatorD::Instance();
  auto trkParam = resTrack.track;
  o2::track::TrackParD trkRef, *refLin = nullptr;
  if (params.useStableRef) {
    refLin = &(trkRef = trkParam);
  }
  float chi2 = 0.f;
  size_t nAccepted = 0;
  double sPrev = clusters.front().s;
  int planePrev = -1; // pad-row plane (azimuthal sector, row, frame rotation) of the last accepted cluster
  for (const auto& oc : clusters) {
    const auto& sc = oc.cl;
    // a 2nd cluster on the same pad-row plane (e.g. shared by both legs, by a leg and its split piece, or on the other side of the
    // central electrode, where the sector of the same azimuth has the same planes) would be a zero-length GBL step, whose transport
    // jacobian is refused as an identity
    const int plane = ((sc.sector % NSectorsPerSide) * o2::tpc::constants::MAXGLOBALPADROW + sc.row) * 2 + oc.rotated;
    if (plane == planePrev) {
      ++mCosmicClustersSkipped;
      continue;
    }
    // a gap, e.g. between the legs: the material crossed in it is lumped at points along it, not all at the next cluster
    if (nAccepted && oc.s - sPrev > 2. * params.cosmicsMaterialStep && !addCosmicMaterialPoints(oc.xGlo, oc.yGlo, trkParam, refLin, resTrack)) {
      return CosmicStatus::PropagationFailed;
    }
    double alpha = o2::math_utils::detail::sector2Angle<double>(sc.sector % NSectorsPerSide);
    double x = sc.x;
    double y = sc.y;
    if (oc.rotated) {
      alpha += o2::constants::math::PI;
      o2::math_utils::detail::bringToPMPi<double>(alpha);
      x = -x;
      y = -y;
    }
    // tried on copies: a cluster the track crosses (almost) parallel to its pad row, typically at the inner end of a leg of a
    // cosmic passing outside the inner field cage, is skipped rather than the whole cosmic; the margin to maxSnp keeps the
    // numerical jacobians of the GBL step to it within the propagation limits
    auto trial = trkParam;
    o2::track::TrackParD trialRef;
    if (refLin) {
      trialRef = *refLin;
    }
    if (!prop->propagateToAlphaX(trial, refLin ? &trialRef : nullptr, alpha, x, false, params.maxSnp, params.maxStep, 1, params.corrType) ||
        std::abs(trial.getSnp()) > params.cosmicsMaxClusterSnp) {
      if (params.verbose > 1) {
        LOGP(info, "Cosmic cluster skipped: sector {} row {} s {} (rotated {}), npoints {}, track {}", sc.sector, sc.row, oc.s, oc.rotated,
             resTrack.info.size(), trkParam.asString());
      }
      ++mCosmicClustersSkipped;
      continue;
    }
    trkParam = trial;
    if (refLin) {
      *refLin = trialRef;
    }
    sPrev = oc.s;
    planePrev = plane;
    std::array<float, 3> cov{0.f, 0.f, 0.f};
    // TODO: as for the collisions, the occupancy / charge components of the error estimation are disabled
    mTPCParam->GetClusterErrors2(sc.sector, sc.row, sc.z, trkParam.getSnp(), trkParam.getTgl(), -1.f, 0.f, 0.f, cov[0], cov[2]);
    mTPCParam->UpdateClusterError2ByState(sc.state, cov[0], cov[2]);
    cov[0] += params.extraClsErrYTPC * params.extraClsErrYTPC;
    cov[2] += params.extraClsErrZTPC * params.extraClsErrZTPC;

    const std::array<double, 2> pos{y, sc.z};
    const std::array<double, 3> covD{cov[0], cov[1], cov[2]};
    chi2 += static_cast<float>(trkParam.getPredictedChi2Quiet(pos, covD));
    if (!trkParam.update(pos, covD)) {
      return CosmicStatus::UpdateFailed;
    }
    if (refLin) { // displace the reference to the last updated cluster
      refLin->setY(pos[0]);
      refLin->setZ(pos[1]);
    }

    auto& pnt = resTrack.info.emplace_back();
    pnt.detIdx = mDetIdx;
    pnt.lr = static_cast<int8_t>(getStack(sc.row));
    pnt.label = Label(mDetIdx, sc.sector, true);
    pnt.x = static_cast<float>(x);
    pnt.alpha = static_cast<float>(alpha);
    pnt.cluster = o2::BaseCluster<float>(static_cast<int16_t>(sc.sector), x, y, sc.z, cov[0], cov[2], cov[1]);
    ++nAccepted;
  }
  if (nAccepted < static_cast<size_t>(params.minTPCClusters)) {
    return CosmicStatus::TooFewClusters;
  }
  resTrack.track = trkParam;
  resTrack.kfFit.chi2 += chi2;
  return CosmicStatus::OK;
}

} // namespace o2::alignrs
