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

// Replay of Mille binaries with a pede solution applied: validation of an alignment w/o regenerating the debug trees.
// Each record is corrected by the solution (r -> r - G.p), its local (track) parameters are refitted as in the final
// pede iteration (Cauchy outlier down-weighting, down-weight fraction cut, chi2 cut at the 3-sigma quantile of
// chi2(ndf)), and for the accepted records the "pressure" on every global parameter is accumulated: the gradient
// term b = sum w g r_fit and the curvature with the local parameters profiled out, Hp = sum w g^2 - u^T A^-1 u.
// b/Hp is the shift the data still ask for, b^2/Hp the chi2 it would gain (each parameter alone). Free parameters of a
// converged fit have b ~ 0, unless they enter a constraint (their common mode is then balanced by the Lagrange
// multiplier); fixed ones tell what freeing them would bring. Groups of parameters (a regex on the volume name and
// a DOF name) are accumulated as one parameter whose derivative is the sum of the members', e.g. the common radial
// shift of all the chambers of a TRD layer.
// See doc/MilleReplay.md.
//
// Built and installed as o2-dev-mille-replay (stand-alone: no O2/ROOT/pede needed).
//   o2-dev-mille-replay -tree mp2tree.txt [-res millepede.res] [-con mp2con.txt] [-groups groups.txt] [-sel regex]
//                 [-dw 4] [-dwcut 0.15] [-cut 1] [-maxloc 1500] [-maxrec N] [-j 8] -o out file1.bin [file2.bin ...]
//   -tree   : mp2tree.txt of the RECORDING run (it lists the labels of all DOFs present in the records)
//   -res    : solution to apply (millepede.res); w/o it the records are replayed at the recorded reference
//   -con    : add the Measurement priors of this mp2con.txt to b and Hp of the parameters
//   -groups : lines "name DOF regex", e.g. "TRD_L5_radial TX TRD/sm[0-9]+/st[0-4]/pl5"
//   -sel    : accumulate the per-parameter pressure only for the volumes matching this regex (default: all)
// Output: out_files.txt (per input file: records, rejection and <chi2/ndf> for records with/without vertex labels),
//         out_pars.txt (per label), out_groups.txt (per group). The sums are additive over input files.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <mutex>
#include <regex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace
{

struct LabelInfo {
  std::string vol;
  std::string dof;
  int slot = -1; // calibration time slot, -1 for rigid-body DOFs
  bool isMV = false;
  bool selected = true;
  std::vector<int> groups;
};

struct Prior {
  double value = 0, sigma = 1;
  std::vector<std::pair<int, double>> terms;
};

struct Options {
  std::string tree, res, con, groups, sel, out;
  std::vector<std::string> inputs;
  int nDW = 4;
  double dwCut = 0.15, chi2Cut = 1.;
  int maxLoc = 1500, nThreads = 1;
  long maxRec = -1;
};

struct ParAcc {
  double b = 0, hu = 0, hp = 0;
  long nRow = 0, nRec = 0;
  void add(const ParAcc& o)
  {
    b += o.b;
    hu += o.hu;
    hp += o.hp;
    nRow += o.nRow;
    nRec += o.nRec;
  }
};

struct FileStat {
  std::string name;
  long nRec = 0, nSkip = 0;
  long n[2] = {0, 0}, nRej[2] = {0, 0}, nRejDW[2] = {0, 0};
  double chi2[2] = {0, 0}, ndf[2] = {0, 0};
};

struct Accumulator {
  std::unordered_map<int, ParAcc> pars;
  std::vector<ParAcc> groups;
  std::vector<FileStat> files;
};

struct Row {
  double meas = 0, sig = 0, w = 0, rfit = 0;
  std::vector<std::pair<int, double>> loc, glo;
};

// ---------------------------------------------------------------------------------------------- input files

// "NAME (kind) RB[ TX(l) ... ] CAL@k[ X(l) ... ]"
std::unordered_map<int, LabelInfo> readLabelTree(const std::string& fn)
{
  std::unordered_map<int, LabelInfo> labels;
  std::ifstream in(fn);
  if (!in) {
    fprintf(stderr, "cannot open %s\n", fn.c_str());
    exit(1);
  }
  const std::regex tok(R"((RB|CAL@(\d+))\[|(\w+)\((\d+)\))");
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::string name;
    ls >> name;
    if (name.empty()) {
      continue;
    }
    int slot = -1;
    for (auto it = std::sregex_iterator(line.begin(), line.end(), tok); it != std::sregex_iterator(); ++it) {
      const auto& m = *it;
      if (m[1].matched) {
        slot = m[2].matched ? std::stoi(m[2].str()) : -1;
        continue;
      }
      LabelInfo li;
      li.vol = name;
      li.dof = m[3].str();
      li.slot = slot;
      li.isMV = name == "PVT/meanVertex";
      labels[std::stoi(m[4].str())] = li;
    }
  }
  return labels;
}

// millepede.res: "label value presigma ..." after a header line
std::unordered_map<int, double> readSolution(const std::string& fn)
{
  std::unordered_map<int, double> sol;
  if (fn.empty()) {
    return sol;
  }
  std::ifstream in(fn);
  if (!in) {
    fprintf(stderr, "cannot open %s\n", fn.c_str());
    exit(1);
  }
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    long lab;
    double val;
    if (ls >> lab >> val && val != 0.) {
      sol[(int)lab] = val;
    }
  }
  return sol;
}

// mp2con.txt: "Measurement value sigma" followed by "label coefficient" lines
std::vector<Prior> readPriors(const std::string& fn)
{
  std::vector<Prior> priors;
  if (fn.empty()) {
    return priors;
  }
  std::ifstream in(fn);
  std::string line;
  bool inMeas = false;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::string first;
    if (!(ls >> first) || first[0] == '!') {
      continue;
    }
    if (first == "Measurement") {
      Prior p;
      ls >> p.value >> p.sigma;
      priors.push_back(p);
      inMeas = true;
    } else if (first == "Constraint" || first == "Wconstraint") {
      inMeas = false;
    } else if (inMeas) {
      double c;
      ls >> c;
      priors.back().terms.push_back({std::stoi(first), c});
    }
  }
  return priors;
}

std::vector<std::string> readGroups(const std::string& fn, std::unordered_map<int, LabelInfo>& labels)
{
  std::vector<std::string> names;
  if (fn.empty()) {
    return names;
  }
  std::ifstream in(fn);
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::string name, dof, rx;
    if (!(ls >> name >> dof >> rx) || name[0] == '#') {
      continue;
    }
    const std::regex re(rx);
    const int ig = names.size();
    names.push_back(name + " " + dof + " " + rx);
    int nm = 0;
    for (auto& [lab, li] : labels) {
      if (li.dof == dof && std::regex_match(li.vol, re)) {
        li.groups.push_back(ig);
        nm++;
      }
    }
    fprintf(stderr, "group %s: %d labels\n", names.back().c_str(), nm);
  }
  return names;
}

// ---------------------------------------------------------------------------------------------- linear algebra

// in-place Cholesky factor (lower, row-major) of the SPD matrix A
bool cholFactor(std::vector<double>& A, int n)
{
  for (int j = 0; j < n; j++) {
    double d = A[j * n + j];
    for (int k = 0; k < j; k++) {
      d -= A[j * n + k] * A[j * n + k];
    }
    if (d <= 1e-30) {
      return false;
    }
    d = std::sqrt(d);
    A[j * n + j] = d;
    for (int i = j + 1; i < n; i++) {
      double s = A[i * n + j];
      for (int k = 0; k < j; k++) {
        s -= A[i * n + k] * A[j * n + k];
      }
      A[i * n + j] = s / d;
    }
  }
  return true;
}

// y = L^-1 v (forward substitution)
void forwardSolve(const std::vector<double>& L, int n, std::vector<double>& v)
{
  for (int i = 0; i < n; i++) {
    double s = v[i];
    for (int k = 0; k < i; k++) {
      s -= L[i * n + k] * v[k];
    }
    v[i] = s / L[i * n + i];
  }
}

// x = L^-T y (backward substitution)
void backwardSolve(const std::vector<double>& L, int n, std::vector<double>& v)
{
  for (int i = n - 1; i >= 0; i--) {
    double s = v[i];
    for (int k = i + 1; k < n; k++) {
      s -= L[k * n + i] * v[k];
    }
    v[i] = s / L[i * n + i];
  }
}

// pede's down-weighting of a measurement with normalised residual z in a local-fit pass: none in pass 1, Huber in
// passes 2-3, Cauchy afterwards
constexpr double HuberC = 1.345, CauchyC = 2.3849;

double downWeight(int pass, double z)
{
  const double az = std::fabs(z);
  if (pass == 1) {
    return 1.;
  }
  if (pass <= 3) {
    return az > HuberC ? HuberC / az : 1.;
  }
  return 1. / (1. + (z / CauchyC) * (z / CauchyC));
}

// the matching chi2 contribution
double downWeightedChi2(int pass, double z)
{
  const double az = std::fabs(z);
  if (pass == 1) {
    return z * z;
  }
  if (pass <= 3) {
    return az > HuberC ? 2. * HuberC * (az - 0.5 * HuberC) : z * z;
  }
  return CauchyC * CauchyC * std::log(1. + (z / CauchyC) * (z / CauchyC));
}

// 3-sigma (p = 0.9973) quantile of chi2(ndf), Wilson-Hilferty approximation
double chi2Quantile3Sigma(int ndf)
{
  const double z = 2.7822, a = 2. / (9. * ndf);
  const double t = 1. - a + z * std::sqrt(a);
  return ndf * t * t * t;
}

// ---------------------------------------------------------------------------------------------- record processing

class RecordProcessor
{
 public:
  RecordProcessor(const Options& opt, const std::unordered_map<int, LabelInfo>& labels, const std::unordered_map<int, double>& sol, int nGroups)
    : mOpt(opt), mLabels(labels), mSol(sol), mNGroups(nGroups) {}

  void processFile(const std::string& fn, Accumulator& acc);

 private:
  bool readRecord(FILE* f);
  void applySolution();
  bool fitLocals(double& chi2, int& ndf, double& dwFrac);
  void accumulatePressure(Accumulator& acc);
  const LabelInfo* info(int lab) const
  {
    auto it = mLabels.find(lab);
    return it == mLabels.end() ? nullptr : &it->second;
  }

  const Options& mOpt;
  const std::unordered_map<int, LabelInfo>& mLabels;
  const std::unordered_map<int, double>& mSol;
  int mNGroups;
  // per-record buffers, reused (performance-critical loop)
  std::vector<Row> mRows;
  size_t mNRows = 0;
  std::vector<float> mBufF;
  std::vector<double> mBufD, mVal;
  std::vector<int> mBufI;
  std::vector<int> mLocIdx; // local label -> index, -1 if absent
  std::vector<int> mLocUsed;
  int mNLoc = 0;
  std::vector<double> mA, mX;
  std::unordered_map<int, int> mParSlot; // parameter key (label, or -1-group) -> slot in the per-record buffers
  std::vector<int> mSlotKey;
  std::vector<std::vector<double>> mU;
  std::vector<ParAcc> mSlotAcc;
};

bool RecordProcessor::readRecord(FILE* f)
{
  int n2;
  if (fread(&n2, 4, 1, f) != 1) {
    return false;
  }
  const bool dbl = n2 < 0;
  const int n = std::abs(n2) / 2;
  mVal.resize(n);
  if (dbl) {
    mBufD.resize(n);
    if (fread(mBufD.data(), 8, n, f) != (size_t)n) {
      return false;
    }
    std::copy(mBufD.begin(), mBufD.end(), mVal.begin());
  } else {
    mBufF.resize(n);
    if (fread(mBufF.data(), 4, n, f) != (size_t)n) {
      return false;
    }
    std::copy(mBufF.begin(), mBufF.end(), mVal.begin());
  }
  mBufI.resize(n);
  if (fread(mBufI.data(), 4, n, f) != (size_t)n) {
    return false;
  }
  // layout after the leading (0,0) word: meas (label 0), local derivatives, sigma (label 0), global derivatives
  mNRows = 0;
  int state = 0;
  for (int i = 1; i < n; i++) {
    const double v = mVal[i];
    const int k = mBufI[i];
    if (k == 0) {
      if (state == 0 || state == 2) {
        if (mNRows == mRows.size()) {
          mRows.emplace_back();
        }
        Row& r = mRows[mNRows++];
        r.meas = v;
        r.sig = 0;
        r.loc.clear();
        r.glo.clear();
        state = 1;
      } else {
        mRows[mNRows - 1].sig = v;
        state = 2;
      }
    } else if (state == 1) {
      mRows[mNRows - 1].loc.push_back({k, v});
    } else {
      mRows[mNRows - 1].glo.push_back({k, v});
    }
  }
  return true;
}

void RecordProcessor::applySolution()
{
  if (mSol.empty()) {
    return;
  }
  for (size_t ir = 0; ir < mNRows; ir++) {
    auto& r = mRows[ir];
    for (const auto& [lab, g] : r.glo) {
      auto it = mSol.find(lab);
      if (it != mSol.end()) {
        r.meas -= g * it->second;
      }
    }
  }
}

// local fit with Cauchy down-weighting as in pede; on success mA holds the Cholesky factor of the last iteration
bool RecordProcessor::fitLocals(double& chi2, int& ndf, double& dwFrac)
{
  for (int l : mLocUsed) {
    mLocIdx[l] = -1;
  }
  mLocUsed.clear();
  int nMeas = 0;
  for (size_t ir = 0; ir < mNRows; ir++) {
    auto& r = mRows[ir];
    r.w = r.sig > 0 ? 1. / (r.sig * r.sig) : 0.;
    nMeas += r.sig > 0;
    for (const auto& [l, d] : r.loc) {
      if (l >= (int)mLocIdx.size()) {
        mLocIdx.resize(l + 1, -1);
      }
      if (mLocIdx[l] < 0) {
        mLocIdx[l] = mLocUsed.size();
        mLocUsed.push_back(l);
      }
    }
  }
  mNLoc = mLocUsed.size();
  const int nl = mNLoc;
  if (nl == 0 || nl > mOpt.maxLoc || nMeas <= nl) {
    return false;
  }
  // as pede (pede.f90, loop of the local fits): pass 1 w/o down-weighting, passes 2-3 Huber, then Cauchy, the weights
  // of a pass from the residuals of the previous one; the chi2 compared to the cut is the down-weighted one of the last
  // pass, the down-weight fraction 1 - <sqrt(factor)>
  const int nPass = std::max(1, mOpt.nDW);
  for (int pass = 1; pass <= nPass; pass++) {
    if (pass > 1) {
      for (size_t ir = 0; ir < mNRows; ir++) {
        auto& r = mRows[ir];
        if (r.sig > 0) {
          r.w = downWeight(pass, r.rfit / r.sig) / (r.sig * r.sig);
        }
      }
    }
    mA.assign((size_t)nl * nl, 0.);
    mX.assign(nl, 0.);
    for (size_t ir = 0; ir < mNRows; ir++) {
      const auto& r = mRows[ir];
      if (r.w <= 0) {
        continue;
      }
      for (const auto& [la, da] : r.loc) {
        const int i = mLocIdx[la];
        mX[i] += r.w * da * r.meas;
        for (const auto& [lb, db] : r.loc) {
          const int j = mLocIdx[lb];
          if (j <= i) {
            mA[(size_t)i * nl + j] += r.w * da * db;
          }
        }
      }
    }
    if (!cholFactor(mA, nl)) {
      return false;
    }
    forwardSolve(mA, nl, mX);
    backwardSolve(mA, nl, mX);
    chi2 = 0;
    double sumSqrtDW = 0;
    for (size_t ir = 0; ir < mNRows; ir++) {
      auto& r = mRows[ir];
      if (r.sig <= 0) {
        continue;
      }
      double p = 0;
      for (const auto& [l, d] : r.loc) {
        p += d * mX[mLocIdx[l]];
      }
      r.rfit = r.meas - p;
      const double z = r.rfit / r.sig;
      chi2 += downWeightedChi2(pass, z);
      sumSqrtDW += std::sqrt(downWeight(pass, z));
    }
    dwFrac = (nMeas - sumSqrtDW) / nMeas;
  }
  ndf = nMeas - nl;
  return true;
}

// b, Hu and the profiled Hp of every selected label and group of an accepted record
void RecordProcessor::accumulatePressure(Accumulator& acc)
{
  const int nl = mNLoc;
  for (auto& [key, slot] : mParSlot) {
    (void)key;
    mSlotAcc[slot] = ParAcc{};
  }
  mParSlot.clear();
  mSlotKey.clear();
  auto slotOf = [&](int key) {
    auto [it, isNew] = mParSlot.emplace(key, (int)mSlotKey.size());
    if (isNew) {
      mSlotKey.push_back(key);
      if (mU.size() < mSlotKey.size()) {
        mU.emplace_back();
        mSlotAcc.emplace_back();
      }
      mU[it->second].assign(nl, 0.);
      mSlotAcc[it->second] = ParAcc{};
    }
    return it->second;
  };
  std::vector<std::pair<int, double>> rowPars; // (slot, derivative) of this row, groups summed
  for (size_t ir = 0; ir < mNRows; ir++) {
    const auto& r = mRows[ir];
    if (r.w <= 0 || r.glo.empty()) {
      continue;
    }
    rowPars.clear();
    for (const auto& [lab, g] : r.glo) {
      const LabelInfo* li = info(lab);
      if (!li) {
        continue;
      }
      if (li->selected) {
        rowPars.push_back({slotOf(lab), g});
      }
      for (int ig : li->groups) {
        const int s = slotOf(-1 - ig);
        auto it = std::find_if(rowPars.begin(), rowPars.end(), [s](const auto& p) { return p.first == s; });
        if (it == rowPars.end()) {
          rowPars.push_back({s, g});
        } else {
          it->second += g;
        }
      }
    }
    for (const auto& [s, g] : rowPars) {
      auto& pa = mSlotAcc[s];
      pa.b += r.w * g * r.rfit;
      pa.hu += r.w * g * g;
      pa.nRow++;
      auto& u = mU[s];
      for (const auto& [l, d] : r.loc) {
        u[mLocIdx[l]] += r.w * g * d;
      }
    }
  }
  for (size_t s = 0; s < mSlotKey.size(); s++) {
    auto& u = mU[s];
    forwardSolve(mA, nl, u); // |L^-1 u|^2 = u^T A^-1 u
    double q = 0;
    for (int i = 0; i < nl; i++) {
      q += u[i] * u[i];
    }
    auto& pa = mSlotAcc[s];
    pa.hp = pa.hu - q;
    pa.nRec = 1;
    const int key = mSlotKey[s];
    if (key >= 0) {
      acc.pars[key].add(pa);
    } else {
      acc.groups[-1 - key].add(pa);
    }
  }
}

void RecordProcessor::processFile(const std::string& fn, Accumulator& acc)
{
  FileStat st;
  st.name = fn;
  FILE* f = fopen(fn.c_str(), "rb");
  if (!f) {
    fprintf(stderr, "cannot open %s\n", fn.c_str());
    acc.files.push_back(st);
    return;
  }
  while ((mOpt.maxRec < 0 || st.nRec < mOpt.maxRec) && readRecord(f)) {
    st.nRec++;
    bool hasPV = false;
    for (size_t ir = 0; ir < mNRows && !hasPV; ir++) {
      for (const auto& [lab, g] : mRows[ir].glo) {
        const LabelInfo* li = info(lab);
        if (li && li->isMV) {
          hasPV = true;
          break;
        }
      }
    }
    const int cat = hasPV;
    applySolution();
    double chi2 = 0, dwFrac = 0;
    int ndf = 0;
    if (!fitLocals(chi2, ndf, dwFrac)) {
      st.nSkip++;
      continue;
    }
    st.n[cat]++;
    if (dwFrac > mOpt.dwCut) {
      st.nRejDW[cat]++;
      continue;
    }
    if (chi2 > mOpt.chi2Cut * chi2Quantile3Sigma(ndf)) {
      st.nRej[cat]++;
      continue;
    }
    st.chi2[cat] += chi2;
    st.ndf[cat] += ndf;
    accumulatePressure(acc);
  }
  fclose(f);
  acc.files.push_back(st);
}

// ---------------------------------------------------------------------------------------------- output

void addPriors(const std::vector<Prior>& priors, const std::unordered_map<int, double>& sol, Accumulator& acc)
{
  for (const auto& p : priors) {
    double pred = 0;
    for (const auto& [lab, c] : p.terms) {
      auto it = sol.find(lab);
      pred += c * (it == sol.end() ? 0. : it->second);
    }
    const double w = 1. / (p.sigma * p.sigma);
    for (const auto& [lab, c] : p.terms) {
      auto& pa = acc.pars[lab];
      pa.b += w * c * (p.value - pred);
      pa.hu += w * c * c;
      pa.hp += w * c * c;
    }
  }
}

void writeOutput(const Options& opt, const Accumulator& acc, const std::unordered_map<int, LabelInfo>& labels,
                 const std::unordered_map<int, double>& sol, const std::vector<std::string>& groupNames)
{
  FILE* ff = fopen((opt.out + "_files.txt").c_str(), "w");
  fprintf(ff, "# file nRec nSkip | noPV: n rejChi2 rejDW <chi2/ndf>acc | PV: n rejChi2 rejDW <chi2/ndf>acc\n");
  for (const auto& s : acc.files) {
    fprintf(ff, "%s %ld %ld | %ld %ld %ld %.4f | %ld %ld %ld %.4f\n", s.name.c_str(), s.nRec, s.nSkip, s.n[0], s.nRej[0], s.nRejDW[0],
            s.ndf[0] > 0 ? s.chi2[0] / s.ndf[0] : 0., s.n[1], s.nRej[1], s.nRejDW[1], s.ndf[1] > 0 ? s.chi2[1] / s.ndf[1] : 0.);
  }
  fclose(ff);
  auto line = [](FILE* fo, const ParAcc& p) {
    const double shift = p.hp > 0 ? p.b / p.hp : 0., dchi2 = p.hp > 0 ? p.b * p.b / p.hp : 0.;
    fprintf(fo, " %ld %ld %.6e %.6e %.6e %+.6e %.4e\n", p.nRec, p.nRow, p.b, p.hu, p.hp, shift, dchi2);
  };
  FILE* fp = fopen((opt.out + "_pars.txt").c_str(), "w");
  fprintf(fp, "# label volume dof slot value nRec nRow b Hu Hp shift=b/Hp dchi2=b^2/Hp   (cm, rad)\n");
  std::map<int, const ParAcc*> sorted;
  for (const auto& [lab, p] : acc.pars) {
    sorted[lab] = &p;
  }
  for (const auto& [lab, p] : sorted) {
    auto it = labels.find(lab);
    auto is = sol.find(lab);
    fprintf(fp, "%d %s %s %d %+.6e", lab, it == labels.end() ? "?" : it->second.vol.c_str(), it == labels.end() ? "?" : it->second.dof.c_str(),
            it == labels.end() ? -1 : it->second.slot, is == sol.end() ? 0. : is->second);
    line(fp, *p);
  }
  fclose(fp);
  FILE* fg = fopen((opt.out + "_groups.txt").c_str(), "w");
  fprintf(fg, "# group(name dof regex) nRec nRow b Hu Hp shift=b/Hp dchi2=b^2/Hp\n");
  for (size_t i = 0; i < groupNames.size(); i++) {
    fprintf(fg, "%s", groupNames[i].c_str());
    line(fg, acc.groups[i]);
  }
  fclose(fg);
}

Options parseArgs(int argc, char** argv)
{
  Options o;
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        fprintf(stderr, "missing value for %s\n", a.c_str());
        exit(1);
      }
      return argv[++i];
    };
    if (a == "-tree") {
      o.tree = next();
    } else if (a == "-res") {
      o.res = next();
    } else if (a == "-con") {
      o.con = next();
    } else if (a == "-groups") {
      o.groups = next();
    } else if (a == "-sel") {
      o.sel = next();
    } else if (a == "-dw") {
      o.nDW = std::stoi(next());
    } else if (a == "-dwcut") {
      o.dwCut = std::stod(next());
    } else if (a == "-cut") {
      o.chi2Cut = std::stod(next());
    } else if (a == "-maxloc") {
      o.maxLoc = std::stoi(next());
    } else if (a == "-maxrec") {
      o.maxRec = std::stol(next());
    } else if (a == "-j") {
      o.nThreads = std::max(1, std::stoi(next()));
    } else if (a == "-o") {
      o.out = next();
    } else {
      o.inputs.push_back(a);
    }
  }
  if (o.tree.empty() || o.out.empty() || o.inputs.empty()) {
    fprintf(stderr, "usage: o2-dev-mille-replay -tree mp2tree.txt [-res millepede.res] [-con mp2con.txt] [-groups groups.txt] [-sel regex]\n"
                    "       [-dw 4] [-dwcut 0.15] [-cut 1] [-maxloc 1500] [-maxrec N] [-j 8] -o outPrefix file1.bin [file2.bin ...]\n");
    exit(1);
  }
  return o;
}

} // namespace

int main(int argc, char** argv)
{
  const Options opt = parseArgs(argc, argv);
  auto labels = readLabelTree(opt.tree);
  if (!opt.sel.empty()) {
    const std::regex re(opt.sel);
    for (auto& [lab, li] : labels) {
      li.selected = std::regex_match(li.vol, re);
    }
  }
  const auto groupNames = readGroups(opt.groups, labels);
  const auto sol = readSolution(opt.res);
  const auto priors = readPriors(opt.con);
  fprintf(stderr, "%zu labels, %zu solution values, %zu priors, %zu groups, %zu input files, %d threads\n", labels.size(), sol.size(),
          priors.size(), groupNames.size(), opt.inputs.size(), opt.nThreads);

  std::vector<Accumulator> accs(opt.nThreads);
  std::atomic<size_t> nextFile{0};
  std::mutex logMutex;
  auto worker = [&](int it) {
    RecordProcessor proc(opt, labels, sol, groupNames.size());
    accs[it].groups.resize(groupNames.size());
    for (size_t i = nextFile++; i < opt.inputs.size(); i = nextFile++) {
      proc.processFile(opt.inputs[i], accs[it]);
      std::lock_guard<std::mutex> lock(logMutex);
      const auto& s = accs[it].files.back();
      fprintf(stderr, "[%zu/%zu] %s: %ld records\n", i + 1, opt.inputs.size(), s.name.c_str(), s.nRec);
    }
  };
  std::vector<std::thread> threads;
  for (int it = 0; it < opt.nThreads; it++) {
    threads.emplace_back(worker, it);
  }
  for (auto& t : threads) {
    t.join();
  }
  Accumulator total;
  total.groups.resize(groupNames.size());
  for (auto& a : accs) {
    for (const auto& [lab, p] : a.pars) {
      total.pars[lab].add(p);
    }
    for (size_t i = 0; i < a.groups.size(); i++) {
      total.groups[i].add(a.groups[i]);
    }
    total.files.insert(total.files.end(), a.files.begin(), a.files.end());
  }
  std::sort(total.files.begin(), total.files.end(), [](const FileStat& a, const FileStat& b) { return a.name < b.name; });
  addPriors(priors, sol, total);
  writeOutput(opt, total, labels, sol, groupNames);
  return 0;
}
