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

// TRD / TOF residuals of the `gbltrack` tree of debug_alg.root (see doc/DebugInfo.md), in an o2::HistoManager.
//
// trdTofResid(): per TRD chamber (sector, layer, stack) dY vs snp (nSnp bins in |snp| < maxSnp) and dZ vs tgl of the track
// (nTglBins bins between Zmin/X and Zmax/X of the chamber, X = nominal tracking X of the layer), per TOF sector and
// Z slice (nTOFSlices equal slices of the strip centre Z, i.e. groups of whole strips) the same, the tgl range being
// the Z limits of the slice divided by the strip radius. All the residuals are in cm. Only points with hasMeas() are used.
// The TH2F are fitted with FitSlicesY (gaus in [-resRange, resRange]), the mean and sigma TH1 are kept in the manager.
//
// o2::HistoManager IDs:   Sect*10000 + Lr*1000 + St*100 + iYZ*10  (+1: mean, +2: sigma of the fitted slices)
//   Sect 0:17, Lr 0:5 = TRD layer, 6 = TOF, St = TRD stack or TOF Z slice, iYZ 0 (Y) or 1 (Z).
//
// drawTRDTOFResid(): superimposes several such managers, the same way as `o2-check-resid-workfow --draw-external-only`
// does it with the ext_hm_list, ext_leg_list, ext_color_list, ext_marker_list keys of CheckResidConfig (see
// doc/reResidAll.sh). One pdf page per sector and layer: 4 rows (mean dY, sigma dY, mean dZ, sigma dZ), one column per
// TRD stack or TOF Z slice.
//
// Usage (O2 environment, libO2AlignRS and the O2Align headers findable, ROOT_INCLUDE_PATH):
//   root -b -q 'trdTofResid.C+("debug_alg.root", "resTRDTOF_Unb_hman.root", "o2sim_geometry-aligned.root", 1)'
//   root -b -q 'trdTofResid.C+("files.txt", "resTRDTOF_Pre_hman.root", "", 0)'    // text file: one root file per line
//   root -l 'trdTofResid.C+' -e 'drawTRDTOFResid("checkresid.ext_hm_list=a_hman.root,b_hman.root;checkresid.ext_leg_list=A,B;checkresid.ext_color_list=kBlack,kRed;checkresid.ext_marker_list=20,24;", "cmp.pdf")'

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TH1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fmt/format.h>
#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "CommonUtils/StringUtils.h"
#include "DetectorsBase/GeometryManager.h"
#include "GlobalTrackingStudy/CheckResidConfig.h"
#include "GlobalTrackingStudy/HistoManager.h"
#include "TOFBase/Geo.h"
#include "TRDBase/Geometry.h"
#include "O2Align/AlignmentDbg.h"

namespace trdtof
{

constexpr int NSect = 18, NTRDLr = 6, NTRDSt = 5, LrTOF = 6, DetTRD = 3, DetTOF = 4;
constexpr int NStripSM = o2::tof::Geo::NSTRIPXSECTOR;
constexpr std::array<const char*, 3> ResTypeName{"Pre", "Unb", "Fit"};

inline int histoID(int sect, int lr, int st, int iyz) { return sect * 10000 + lr * 1000 + st * 100 + iyz * 10; }

/// strips grouped in Z slices of the TOF sector, from the ideal geometry
struct TOFSlices {
  int nSlices = 0;
  std::vector<int> sliceOfStrip;             // by continuous strip number, -1 for missing strips
  std::vector<std::array<float, 3>> zLoHiR;  // per slice: Z min, Z max of the slice, mean radius of its strips
  std::vector<std::array<int, 2>> stripRange; // per slice and sector (index slice*NSect + sect): first, last strip number in the sector
};

/// chain of the `gbltrack` trees from a root file or from a text file with one root file per line (not recursive)
TChain* loadDbgChain(const char* inpData, const char* chName)
{
  auto* chain = new TChain(chName);
  TString inp = inpData;
  if (inp.EndsWith(".root")) {
    chain->AddFile(inpData);
  } else {
    std::ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n", inpData);
      delete chain;
      return nullptr;
    }
    TString flName;
    while (flName.ReadLine(inpf)) {
      flName = flName.Strip(TString::kBoth, ' ');
      if (flName.IsNull() || flName.BeginsWith("//") || flName.BeginsWith("#")) {
        continue;
      }
      flName = flName.Strip(TString::kBoth, ',');
      flName = flName.Strip(TString::kBoth, '"');
      printf("Adding %s\n", flName.Data());
      chain->AddFile(flName.Data());
    }
  }
  if (chain->GetEntries() < 1) {
    printf("Obtained chain is empty\n");
    delete chain;
    return nullptr;
  }
  printf("Opened %s chain with %lld entries\n", chName, chain->GetEntries());
  return chain;
}

TOFSlices makeTOFSlices(int nSlices)
{
  TOFSlices ts;
  ts.nSlices = nSlices;
  ts.sliceOfStrip.assign(NSect * NStripSM, -1);
  std::vector<float> zc(NSect * NStripSM, 0.f), rc(NSect * NStripSM, 0.f);
  float zMin = 1e9, zMax = -1e9;
  for (int isect = 0; isect < NSect; ++isect) {
    for (int istr = 0; istr < NStripSM; ++istr) {
      auto* pne = gGeoManager->GetAlignableEntry(Form("TOF/sm%02d/strip%02d", isect, istr + 1));
      if (!pne) {
        continue;
      }
      auto* pn = pne->GetPhysicalNode();
      if (!pn) {
        pn = gGeoManager->MakeAlignablePN(pne);
      }
      const double loc[3] = {0, 0, 0};
      double glo[3];
      pn->GetMatrix()->LocalToMaster(loc, glo);
      const int sid = isect * NStripSM + istr;
      zc[sid] = glo[2];
      rc[sid] = std::hypot(glo[0], glo[1]);
      zMin = std::min(zMin, zc[sid]);
      zMax = std::max(zMax, zc[sid]);
    }
  }
  const float halfPitch = 0.5f * (zMax - zMin) / (NStripSM - 1);
  zMin -= halfPitch;
  zMax += halfPitch;
  const float dz = (zMax - zMin) / nSlices;
  ts.zLoHiR.assign(nSlices, {0.f, 0.f, 0.f});
  std::vector<int> nStr(nSlices, 0);
  ts.stripRange.assign(nSlices * NSect, {NStripSM, -1});
  for (int is = 0; is < nSlices; ++is) {
    ts.zLoHiR[is][0] = zMin + is * dz;
    ts.zLoHiR[is][1] = zMin + (is + 1) * dz;
  }
  for (int sid = 0; sid < NSect * NStripSM; ++sid) {
    if (ts.sliceOfStrip[sid] < 0 && rc[sid] > 0.f) {
      const int is = std::min(nSlices - 1, std::max(0, int((zc[sid] - zMin) / dz)));
      ts.sliceOfStrip[sid] = is;
      ts.zLoHiR[is][2] += rc[sid];
      nStr[is]++;
      auto& rng = ts.stripRange[is * NSect + sid / NStripSM];
      rng[0] = std::min(rng[0], sid % NStripSM);
      rng[1] = std::max(rng[1], sid % NStripSM);
    }
  }
  for (int is = 0; is < nSlices; ++is) {
    ts.zLoHiR[is][2] = nStr[is] ? ts.zLoHiR[is][2] / nStr[is] : 380.f;
  }
  return ts;
}

/// books dY vs snp and dZ vs tgl for one chamber / slice, tgl in [tglMin, tglMax]
void bookPair(o2::HistoManager& hm, int id, const std::string& name, const std::string& title, int resType, int nSnp, float maxSnp,
              int nTgl, float tglMin, float tglMax, int nBinRes, float resRangeY, float resRangeZ)
{
  for (int iyz = 0; iyz < 2; ++iyz) {
    const std::string rname = fmt::format("d{}{}", iyz ? 'Z' : 'Y', ResTypeName[resType]);
    const std::string nm = fmt::format("{}_{}", name, rname);
    const std::string tit = fmt::format("{} {}", title, rname);
    const std::string ytit = fmt::format("{} (cm)", iyz ? "#DeltaZ" : "#DeltaY");
    TH2F* h2 = iyz == 0 ? new TH2F(nm.c_str(), fmt::format("{};snp;{}", tit, ytit).c_str(), nSnp, -maxSnp, maxSnp, nBinRes, -resRangeY, resRangeY)
                        : new TH2F(nm.c_str(), fmt::format("{};tg#lambda;{}", tit, ytit).c_str(), nTgl, tglMin, tglMax, nBinRes, -resRangeZ, resRangeZ);
    hm.addHisto(h2, id + iyz * 10);
  }
}

void bookHistos(o2::HistoManager& hm, const TOFSlices& ts, int resType, int nSnp, float maxSnp, int nTglBins, int nBinRes, float resRangeY, float resRangeZ)
{
  auto* trdGeo = o2::trd::Geometry::instance();
  for (int isect = 0; isect < NSect; ++isect) {
    for (int ilr = 0; ilr < NTRDLr; ++ilr) {
      for (int ist = 0; ist < NTRDSt; ++ist) {
        const auto symChamber = fmt::format("TRD/sm{:02d}/st{}/pl{}", isect, ist, ilr);
        if (!gGeoManager->GetAlignableEntry(symChamber.c_str())) { // missing chamber
          continue;
        }
        const int det = o2::trd::Geometry::getDetector(ilr, ist, isect);
        const float x = o2::trd::Geometry::getTime0(ilr);
        const float z0 = trdGeo->getRowEnd(ilr, ist), z1 = trdGeo->getRow0(ilr, ist);
        bookPair(hm, histoID(isect, ilr, ist, 0), fmt::format("trd_s{:02d}_l{}_st{}", isect, ilr, ist), fmt::format("{} ch{}", symChamber, det),
                 resType, nSnp, maxSnp, nTglBins, std::min(z0, z1) / x, std::max(z0, z1) / x, nBinRes, resRangeY, resRangeZ);
      }
    }
    for (int is = 0; is < ts.nSlices; ++is) {
      const auto& rng = ts.stripRange[is * NSect + isect];
      if (rng[1] < 0) {
        continue;
      }
      const auto& lim = ts.zLoHiR[is];
      bookPair(hm, histoID(isect, LrTOF, is, 0), fmt::format("tof_s{:02d}_z{}", isect, is), fmt::format("TOF/sm{:02d}/strip{}:{}", isect, rng[0], rng[1]),
               resType, nSnp, maxSnp, nTglBins, lim[0] / lim[2], lim[1] / lim[2], nBinRes, resRangeY, resRangeZ);
    }
  }
}

/// the ID of the TH2F receiving the point, -1 if the point is not TRD/TOF or has no booked histogram
int pointToHistoID(const o2::alignrs::DbgPoint& p, const TOFSlices& ts)
{
  if (p.det == DetTRD && p.sens >= 0) {
    return histoID(p.sens / (NTRDLr * NTRDSt), p.sens % NTRDLr, (p.sens % (NTRDLr * NTRDSt)) / NTRDLr, 0);
  }
  if (p.det == DetTOF && p.sens >= 0 && p.sens < (int)ts.sliceOfStrip.size() && ts.sliceOfStrip[p.sens] >= 0) {
    return histoID(p.sens / NStripSM, LrTOF, ts.sliceOfStrip[p.sens], 0);
  }
  return -1;
}

void fillHistos(TChain* chain, o2::HistoManager& hm, const TOFSlices& ts, int resType, bool onlyAccepted, long maxTracks)
{
  auto* dt = new o2::alignrs::DbgTrack;
  chain->SetBranchAddress("t", &dt);
  const long nEnt = maxTracks > 0 ? std::min<long>(maxTracks, chain->GetEntries()) : chain->GetEntries();
  for (long i = 0; i < nEnt; ++i) {
    chain->GetEntry(i);
    if (onlyAccepted && !dt->accepted) {
      continue;
    }
    const float tgl = dt->getTgl();
    for (const auto& p : dt->points) {
      if (!p.hasMeas()) {
        continue;
      }
      const int id = pointToHistoID(p, ts);
      if (id < 0) {
        continue;
      }
      float dy = p.dyPre, dz = p.dzPre;
      if (resType == 1) {
        if (p.sigYUnb <= 0.f) { // the unbiased residuals were not computed
          continue;
        }
        dy = p.dyUnb;
        dz = p.dzUnb;
      } else if (resType == 2) {
        dy = p.dyFit;
        dz = p.dzFit;
      }
      if (auto* hy = hm.getHisto2F(id)) {
        hy->Fill(p.snp, dy);
      }
      if (auto* hz = hm.getHisto2F(id + 10)) {
        hz->Fill(tgl, dz);
      }
    }
  }
  delete dt;
  chain->ResetBranchAddresses();
}

/// profile of the 2D histogram converted to TH1D: mean (error of the mean) or sigma (spread, error sigma/sqrt(2n)) per bin; empty bins have zero error
TH1* profileToTH1(TH2F* h2, bool sigma, const char* name)
{
  auto* prf = h2->ProfileX(Form("%s_pfx%d", h2->GetName(), (int)sigma), 1, -1, sigma ? "s" : "");
  auto* h1 = new TH1D(name, h2->GetTitle(), h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
  h1->SetDirectory(nullptr);
  for (int i = 1; i <= h1->GetNbinsX(); ++i) {
    const double n = prf->GetBinEntries(i);
    if (n < 2) { // the spread of a single entry is undefined
      continue;
    }
    h1->SetBinContent(i, sigma ? prf->GetBinError(i) : prf->GetBinContent(i));
    h1->SetBinError(i, sigma ? prf->GetBinError(i) / std::sqrt(2. * n) : prf->GetBinError(i));
  }
  delete prf;
  return h1;
}

/// mean/sigma of all the base 2D histograms of the manager, added with IDs +1/+2: FitSlicesY (gaus within the Y range of the
/// histogram) if it has at least minEntriesFit entries, otherwise the profile converted to TH1D (mean and spread)
void postProcessHistos(o2::HistoManager& hm, int minSliceStat, int minEntriesFit)
{
  printf("Fitting histos\n");
  TObjArray arr;
  arr.SetOwner(true);
  for (int isect = 0; isect < NSect; ++isect) {
    for (int ilr = 0; ilr <= LrTOF; ++ilr) {
      for (int ist = 0; ist < (ilr == LrTOF ? 10 : NTRDSt); ++ist) {
        for (int iyz = 0; iyz < 2; ++iyz) {
          const int id = histoID(isect, ilr, ist, iyz);
          auto* h2 = hm.getHisto2F(id);
          if (!h2 || h2->GetEntries() < 1) {
            continue;
          }
          TH1 *hmean = nullptr, *hsig = nullptr;
          if (h2->GetEntries() < minEntriesFit) {
            hmean = profileToTH1(h2, false, Form("%s_1", h2->GetName()));
            hsig = profileToTH1(h2, true, Form("%s_2", h2->GetName()));
          } else {
            const double rng = h2->GetYaxis()->GetXmax();
            TF1 gs("gs", "gaus", -rng, rng);
            h2->FitSlicesY(&gs, 0, -1, minSliceStat, "QNR", &arr);
            hmean = (TH1*)arr.RemoveAt(1);
            hsig = (TH1*)arr.RemoveAt(2);
            arr.Delete();
          }
          if (hmean) {
            hmean->SetTitle(Form("<%s>", h2->GetTitle()));
            hmean->SetXTitle(h2->GetXaxis()->GetTitle());
            hmean->SetYTitle(Form("<%s>", h2->GetYaxis()->GetTitle()));
            hm.addHisto(hmean, id + 1);
          }
          if (hsig) {
            hsig->SetTitle(Form("#sigma(%s)", h2->GetTitle()));
            hsig->SetXTitle(h2->GetXaxis()->GetTitle());
            hsig->SetYTitle(Form("#sigma(%s)", h2->GetYaxis()->GetTitle()));
            hm.addHisto(hsig, id + 2);
          }
        }
      }
    }
  }
}

} // namespace trdtof

/// Fill and fit the TRD/TOF residual histograms.
/// inpName   : debug_alg.root, or a text file with one root file per line
/// outName   : output file of the o2::HistoManager
/// geomName  : geometry file (alignable entries select the existing chambers/strips and give the TOF strip positions), "" = default
/// resType   : 0: dPre, 1: dUnb, 2: dFit
/// maxSnp, nSnp          : |snp| range and number of snp bins
/// nTglBins              : number of tgl bins (limits from the Z extent of the chamber / TOF slice divided by X)
/// nTOFSlices            : number of TOF Z slices
/// nBinRes, resRangeY/Z   : bins and half-ranges (cm) of the residual axes of dY and dZ, also the ranges of the gaussian fits
/// minEntriesFit         : histograms with fewer entries use the profile (mean and spread) instead of FitSlicesY
/// minSliceStat          : minimal number of non-empty bins in the slice to be fitted
/// onlyAccepted          : skip the tracks rejected by the chi2 cut of the GBL fit
/// maxTracks             : process at most that many tracks (<=0: all)
void trdTofResid(const char* inpName, const char* outName = "trdTofResid_hman.root", const char* geomName = "", int resType = 0,
                 float maxSnp = 0.4, int nSnp = 7, int nTglBins = 5, int nTOFSlices = 5, int nBinRes = 100, float resRangeY = 3,
                 float resRangeZ = 5, int minSliceStat = 3, int minEntriesFit = 100, bool onlyAccepted = true, long maxTracks = -1)
{
  using namespace trdtof;
  if (nTOFSlices < 1 || nTOFSlices > 10) {
    printf("nTOFSlices must be within 1:10\n");
    return;
  }
  if (resType < 0 || resType > 2) {
    printf("resType must be 0 (dPre), 1 (dUnb) or 2 (dFit)\n");
    return;
  }
  o2::base::GeometryManager::loadGeometry(geomName);
  auto* chain = loadDbgChain(inpName, "gbltrack");
  if (!chain) {
    return;
  }
  const auto ts = makeTOFSlices(nTOFSlices);
  o2::HistoManager hm("", outName);
  hm.SetName(gSystem->BaseName(outName));
  bookHistos(hm, ts, resType, nSnp, maxSnp, nTglBins, nBinRes, resRangeY, resRangeZ);
  fillHistos(chain, hm, ts, resType, onlyAccepted, maxTracks);
  delete chain;
  postProcessHistos(hm, minSliceStat, minEntriesFit);
  hm.write(outName);
  printf("Saved %d histograms to %s\n", hm.getNHistos(), outName);
}

//____________________________________________________________________________________________
// Drawing

namespace trdtof
{
// the same EColor name parsing as CheckResidSpec::drawHistos
std::optional<int> parseColor(const std::string& token)
{
  static constexpr std::array<std::pair<std::string_view, int>, 15> rootColors{{{"kWhite", kWhite}, {"kBlack", kBlack}, {"kGray", kGray}, {"kRed", kRed}, {"kGreen", kGreen}, {"kBlue", kBlue}, {"kYellow", kYellow}, {"kMagenta", kMagenta}, {"kCyan", kCyan}, {"kOrange", kOrange}, {"kSpring", kSpring}, {"kTeal", kTeal}, {"kAzure", kAzure}, {"kViolet", kViolet}, {"kPink", kPink}}};
  std::string normalized;
  for (unsigned char ch : token) {
    if (!std::isspace(ch)) {
      normalized.push_back(ch);
    }
  }
  for (const auto& [name, value] : rootColors) {
    if (normalized.compare(0, name.size(), name) != 0) {
      continue;
    }
    const auto increment = std::string_view(normalized).substr(name.size());
    if (increment.empty()) {
      return value;
    }
    if (increment.size() < 2 || (increment.front() != '+' && increment.front() != '-')) {
      return std::nullopt;
    }
    int offset = 0;
    const auto* begin = increment.data() + 1;
    const auto* end = increment.data() + increment.size();
    auto [parsedEnd, error] = std::from_chars(begin, end, offset);
    if (error != std::errc{} || parsedEnd != end) {
      return std::nullopt;
    }
    return increment.front() == '-' ? value - offset : value + offset;
  }
  return std::nullopt;
}

/// range of a row of pads (one sector and layer, all stacks/slices and managers) from the bins with a successful fit:
/// mean: +-(|median| + 4 MAD); sigma: [max(0, median - 4 MAD), median + 4 MAD]. Returns {0, 0} if there are no bins.
std::array<float, 2> rowRange(const std::vector<std::unique_ptr<o2::HistoManager>>& hmv, int isect, int ilr, int nCol, int iyz, bool isSigma)
{
  std::vector<double> vals;
  for (int ic = 0; ic < nCol; ++ic) {
    for (const auto& hm : hmv) {
      auto* h = hm ? hm->getHisto(histoID(isect, ilr, ic, iyz) + 1 + isSigma) : nullptr;
      if (!h) {
        continue;
      }
      for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
        if (h->GetBinError(ib) > 0.) {
          vals.push_back(h->GetBinContent(ib));
        }
      }
    }
  }
  if (vals.empty()) {
    return {0.f, 0.f};
  }
  auto median = [](std::vector<double>& v) {
    const size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end());
    double m = v[n];
    if (v.size() % 2 == 0) {
      m = 0.5 * (m + *std::max_element(v.begin(), v.begin() + n));
    }
    return m;
  };
  const double med = median(vals);
  for (auto& v : vals) {
    v = std::abs(v - med);
  }
  const double mad = median(vals);
  if (isSigma) {
    return {(float)std::max(0., med - 4. * mad), (float)(med + 4. * mad)};
  }
  const float lim = std::abs(med) + 4. * mad;
  return {-lim, lim};
}

/// draws, superimposed, the histograms with the given ID of all the managers on the current pad, skipping the bins w/o a successful
/// fit; yRange (lo < hi) is the vertical range, otherwise it is set from the content.
/// Returns the graph of each manager actually drawn (the pad owns them).
void drawOnPad(const std::vector<std::unique_ptr<o2::HistoManager>>& hmv, int id, bool isSigma, std::array<float, 2> yRange, const std::vector<int>& colors, const std::vector<int>& markers, std::vector<TGraphErrors*>& drawn)
{
  drawn.assign(hmv.size(), nullptr);
  float vmin = 1e9, vmax = -1e9;
  TH1* frame = nullptr;
  for (size_t i = 0; i < hmv.size(); ++i) {
    auto* h = hmv[i] ? hmv[i]->getHisto(id) : nullptr;
    if (!h) {
      continue;
    }
    auto* gr = new TGraphErrors();
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
      if (h->GetBinError(ib) > 0.) {
        const int n = gr->GetN();
        gr->SetPoint(n, h->GetBinCenter(ib), h->GetBinContent(ib));
        gr->SetPointError(n, 0., h->GetBinError(ib));
        vmin = std::min(vmin, (float)(h->GetBinContent(ib) - h->GetBinError(ib)));
        vmax = std::max(vmax, (float)(h->GetBinContent(ib) + h->GetBinError(ib)));
      }
    }
    gr->SetMarkerStyle(markers[i]);
    gr->SetMarkerSize(1.0);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(colors[i]);
    gr->SetLineColor(colors[i]);
    drawn[i] = gr;
    if (!frame) {
      frame = new TH1F(Form("%s_frame", h->GetName()), h->GetTitle(), h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
      frame->SetDirectory(nullptr);
      frame->SetStats(false);
      frame->SetXTitle(h->GetXaxis()->GetTitle());
      frame->SetYTitle(h->GetYaxis()->GetTitle());
    }
  }
  if (!frame) {
    return;
  }
  float lo = vmin, hi = vmax;
  if (yRange[1] > yRange[0]) {
    lo = yRange[0];
    hi = yRange[1];
  } else if (vmax > vmin) {
    const float mrg = 0.15f * (vmax - vmin);
    lo = isSigma ? 0.f : vmin - mrg;
    hi = vmax + mrg;
  } else {
    lo = isSigma ? 0.f : -1.f;
    hi = 1.f;
  }
  frame->SetMinimum(lo);
  frame->SetMaximum(hi);
  for (int ib = 1; ib <= frame->GetNbinsX(); ++ib) { // empty frame whose line is clipped away
    frame->SetBinContent(ib, lo - 10. * (hi - lo));
  }
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.03);
  for (auto* ax : {frame->GetXaxis(), frame->GetYaxis()}) {
    ax->SetTitleSize(0.065);
    ax->SetLabelSize(0.055);
  }
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetLineColorAlpha(kWhite, 0.);
  frame->Draw();
  for (auto* gr : drawn) {
    if (gr && gr->GetN()) {
      gr->Draw("P");
    }
  }
  gPad->SetGrid();
}
} // namespace trdtof

/// Superimpose the HistoManagers listed in the CheckResidConfig keys ext_hm_list, ext_leg_list, ext_color_list, ext_marker_list
/// (same format as in the --configKeyValues of o2-check-resid-workfow --draw-external-only, e.g. in doc/reResidAll.sh), written
/// by trdTofResid, into pdfName: one page per sector and TRD layer / TOF, 4 rows (mean dY, sigma dY, mean dZ, sigma dZ) and
/// one column per TRD stack / TOF Z slice.
/// cfgKeyValues : "checkresid.ext_hm_list=a.root,b.root;checkresid.ext_leg_list=...;..." (empty: the keys are taken as they are)
/// rerunPostProcess : refit/reconvert the 2D histograms of the loaded managers (minSliceStat, minEntriesFit as in trdTofResid)
/// maxMeanY, maxSigY, maxMeanZ, maxSigZ : fixed vertical ranges (cm), 0 = automatic: per row of a page (all stacks and managers)
///   median and MAD of the plotted bins, mean: +-(|median| + 4 MAD), sigma: [max(0, median - 4 MAD), median + 4 MAD]
void drawTRDTOFResid(const char* cfgKeyValues = "", const char* pdfName = "trdTofResid.pdf", bool rerunPostProcess = false, float maxMeanY = 0.f,
                     float maxSigY = 0.f, float maxMeanZ = 0.f, float maxSigZ = 0.f, int sectMin = 0, int sectMax = 17, int minSliceStat = 3,
                     int minEntriesFit = 100)
{
  using namespace trdtof;
  if (cfgKeyValues && cfgKeyValues[0]) {
    o2::conf::ConfigurableParam::updateFromString(cfgKeyValues);
  }
  const auto& params = o2::checkresid::CheckResidConfig::Instance();
  const auto hmNames = o2::utils::Str::tokenize(params.ext_hm_list, ",;");
  const auto vecLegends = o2::utils::Str::tokenize(params.ext_leg_list, ',');
  const auto colorNames = o2::utils::Str::tokenize(params.ext_color_list, ",;");
  const auto markerIDs = o2::utils::Str::tokenize(params.ext_marker_list, ",;");
  if (hmNames.empty()) {
    printf("checkresid.ext_hm_list is empty, nothing to draw\n");
    return;
  }
  const int nhm = hmNames.size();
  bool useLeg = vecLegends.size() == hmNames.size();
  bool useColors = colorNames.size() == hmNames.size() && !colorNames.empty();
  bool useMarkers = markerIDs.size() == hmNames.size() && !markerIDs.empty();
  if (!useLeg) {
    printf("Warning: %zu legend names provided for %d external histomanagers, will use file names as legends\n", vecLegends.size(), nhm);
  }
  if (!useColors) {
    printf("Warning: %zu color names provided for %d external histomanagers, will use automatic colors\n", colorNames.size(), nhm);
  }
  if (!useMarkers) {
    printf("Warning: %zu marker IDs provided for %d external histomanagers, will use automatic markers\n", markerIDs.size(), nhm);
  }
  std::vector<int> colors, markers;
  if (useColors) { // an unrecognized name falls back to the automatic color of that entry only
    for (size_t i = 0; i < colorNames.size(); ++i) {
      auto c = parseColor(colorNames[i]);
      if (!c) {
        printf("Warning: unrecognized EColor name '%s', will use the automatic color for entry %zu\n", colorNames[i].c_str(), i);
      }
      colors.push_back(c ? *c : -1);
    }
  }
  if (useMarkers) {
    for (const auto& mk : markerIDs) {
      int marker = 0;
      auto [parsedEnd, error] = std::from_chars(mk.data(), mk.data() + mk.size(), marker);
      if (error != std::errc{} || parsedEnd != mk.data() + mk.size()) {
        printf("Warning: unrecognized marker ID '%s', will use automatic markers\n", mk.c_str());
        useMarkers = false;
        break;
      }
      markers.push_back(marker);
    }
  }
  constexpr std::array<unsigned int, 3> hcol{kRed, kBlue, kGreen + 2};
  std::vector<std::unique_ptr<o2::HistoManager>> hmv;
  std::vector<int> hmColor(nhm), hmMarker(nhm);
  for (int i = 0; i < nhm; ++i) {
    printf("Loading external HistoManager %s\n", hmNames[i].c_str());
    auto& hm = hmv.emplace_back(std::make_unique<o2::HistoManager>("", hmNames[i], true));
    if (hm->GetLast() < 1) {
      printf("Failed to load histograms from %s\n", hmNames[i].c_str());
      hm.reset();
      continue;
    }
    if (rerunPostProcess) {
      postProcessHistos(*hm, minSliceStat, minEntriesFit);
    }
    hm->SetName(useLeg ? vecLegends[i].c_str() : hmNames[i].c_str());
    hmColor[i] = useColors && colors[i] >= 0 ? colors[i] : hcol[i % hcol.size()];
    hmMarker[i] = useMarkers ? markers[i] : 20 + i + (i % 2) * 4;
    hm->setMarkerStyle(hmMarker[i], 0.6);
    hm->setColor(hmColor[i]);
  }

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  // pages: sector and layer with at least one histogram, nCol = number of columns (TRD stacks or TOF Z slices) to draw
  struct Page {
    int sect, lr, nCol;
  };
  std::vector<Page> pages;
  for (int isect = sectMin; isect <= sectMax; ++isect) {
    for (int ilr = 0; ilr <= LrTOF; ++ilr) {
      int nCol = 0;
      for (int ist = 0; ist < (ilr == LrTOF ? 10 : NTRDSt); ++ist) {
        for (const auto& hm : hmv) {
          if (hm && hm->getHisto(histoID(isect, ilr, ist, 0) + 1)) {
            nCol = ist + 1;
          }
        }
      }
      if (nCol) {
        nCol = ilr == LrTOF ? nCol : NTRDSt; // all the stacks of a TRD layer
        pages.push_back({isect, ilr, nCol});
      }
    }
  }
  if (pages.empty()) {
    printf("No histograms to draw\n");
    return;
  }
  TCanvas canv("ctrdtof", "", 2000, 2600); // same aspect ratio as the ROOT pdf page (20x26), else the content is misplaced
  canv.SetCanvasSize(2000, 2600);              // the exact size, w/o it the pdf content is shifted
  const float rowMax[4] = {maxMeanY, maxSigY, maxMeanZ, maxSigZ};
  for (size_t ip = 0; ip < pages.size(); ++ip) {
    const auto [isect, ilr, nCol] = pages[ip];
    canv.Clear();
    canv.Divide(nCol, 4, 0.001, 0.001);
    bool legDrawn = false;
    for (int ir = 0; ir < 4; ++ir) {
      const int iyz = ir / 2, isSigma = ir % 2;
      std::array<float, 2> yRange = rowRange(hmv, isect, ilr, nCol, iyz, isSigma);
      if (rowMax[ir] > 0.f) { // explicit fixed range overrides the median/MAD one
        yRange = {isSigma ? 0.f : -rowMax[ir], rowMax[ir]};
      }
      for (int ic = 0; ic < nCol; ++ic) {
        canv.cd(ir * nCol + ic + 1);
        std::vector<TGraphErrors*> drawn;
        drawOnPad(hmv, histoID(isect, ilr, ic, iyz) + 1 + isSigma, isSigma, yRange, hmColor, hmMarker, drawn);
        if (legDrawn || std::none_of(drawn.begin(), drawn.end(), [](const TGraphErrors* g) { return g != nullptr; })) {
          continue;
        }
        auto* lg = new TLegend(0.12, 0.13, 0.9, 0.13 + std::min(0.5f, nhm * 0.2f / 3.f)); // owned by the pad
        lg->SetNColumns(std::max(1, int(std::round(nhm/2))));
        lg->SetFillStyle(0);
        lg->SetBorderSize(0);
        for (int i = 0; i < nhm; ++i) {
          if (hmv[i]) {
            auto* mk = new TGraph(); // legend marker
            mk->SetMarkerStyle(hmMarker[i]);
            mk->SetMarkerColor(hmColor[i]);
            mk->SetLineColor(hmColor[i]);
            lg->AddEntry(mk, hmv[i]->GetName(), "lp")->SetTextColor(hmColor[i]);
          }
        }
        lg->Draw();
        legDrawn = true;
      }
    }
    canv.cd(0);
    const char* flag = ip == 0 ? "(" : (ip + 1 == pages.size() ? ")" : "");
    canv.Print(pages.size() == 1 ? pdfName : Form("%s%s", pdfName, flag));
  }
  const int nPages = pages.size();
  printf("Wrote %d pages to %s\n", nPages, pdfName);
}

/// Called by ROOT when the macro is given as `root trdTofResid.C+ -e '...'` (the function named as the file is executed): only loads the macro.
void trdTofResid()
{
  printf("trdTofResid(inp, out, geom, resType, ...) fills and fits, drawTRDTOFResid(cfgKeyValues, pdf, ...) draws, see the header of trdTofResid.C\n");
}
