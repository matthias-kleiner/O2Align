// Draws the weak-mode diagnostics written by the workflow with --output MilleRes,WeakModesDiag
// (AlignParams.weakModesOutBase, default weakModes.root). One summary page and one page of the
// per-sensor displacement profiles of every mode per detector branch go to <out>.pdf.
//
//   root -l -b -q 'drawWeakModes.C("weakModes.root", "weakModes.pdf")'
//
// Black points: fitted displacement of the sensor centres, red: the part of it the fitted mode
// (plus, for the radial/phi components, the global shift/rotation) explains. um everywhere.
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TTree.h>
#include <map>
#include <memory>
#include <string>
#include <vector>
#endif

namespace
{
// uMode[3*m+k] is the Cartesian contribution of mode column m. The column order follows
// Volume::NGlobalEq, so it must never be hard-coded here: modeCols() reads the names from the
// "modes" tree and colOf() resolves them, which keeps the macro valid when a mode is added.
std::map<std::string, int> gCol;

void modeCols(TTree* modes)
{
  gCol.clear();
  std::string mode, *pm = &mode;
  int b = -1;
  modes->SetBranchAddress("mode", &pm);
  modes->SetBranchAddress("branch", &b);
  for (Long64_t i = 0; i < modes->GetEntries(); ++i) {
    modes->GetEntry(i);
    if (gCol.find(mode) == gCol.end()) {
      gCol[mode] = static_cast<int>(gCol.size());
    }
  }
  modes->ResetBranchAddresses();
}

// -1 if this file has no such mode (e.g. written before the mode existed)
int colOf(const char* name)
{
  const auto it = gCol.find(name);
  return it == gCol.end() ? -1 : it->second;
}

std::string radialOf(int m) { return m < 0 ? "0" : Form("(uMode[%d]*x+uMode[%d]*y)/r", 3 * m, 3 * m + 1); }
std::string phiOf(int m) { return m < 0 ? "0" : Form("(-uMode[%d]*y+uMode[%d]*x)/r", 3 * m, 3 * m + 1); }
std::string zOf(int m) { return m < 0 ? "0" : Form("uMode[%d]", 3 * m + 2); }
std::string radialOf(const char* n) { return radialOf(colOf(n)); }
std::string phiOf(const char* n) { return phiOf(colOf(n)); }
std::string zOf(const char* n) { return zOf(colOf(n)); }

// sum of one component (0: radial, 1: r*dphi, 2: z) over every equation of a family present in the file
std::string familySum(const char* family, int comp)
{
  std::string s = "0";
  for (const auto& [name, col] : gCol) {
    if (name.rfind(family, 0) == 0) {
      s += "+" + (comp == 0 ? radialOf(col) : (comp == 1 ? phiOf(col) : zOf(col)));
    }
  }
  return "(" + s + ")";
}
std::string familySum(const char* family, bool radial) { return familySum(family, radial ? 0 : 1); }

std::unique_ptr<TGraph> graphOf(TTree* t, const std::string& yx, const std::string& cut, int color, int marker, double size)
{
  const auto n = t->Draw(yx.c_str(), cut.c_str(), "goff");
  auto g = std::make_unique<TGraph>(n > 0 ? n : 0, t->GetV2(), t->GetV1());
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(size);
  return g;
}

// data and model points of one panel: yData and yModel in um vs x
void panel(TTree* t, const std::string& cut, const char* title, const char* xTitle, const char* xExpr, const std::string& yData, const std::string& yModel)
{
  auto* mg = new TMultiGraph();
  mg->SetTitle(Form("%s;%s;#Delta [#mum]", title, xTitle));
  mg->Add(graphOf(t, Form("1e4*(%s):%s", yData.c_str(), xExpr), cut, kBlack, 20, 0.35).release(), "p");
  mg->Add(graphOf(t, Form("1e4*(%s):%s", yModel.c_str(), xExpr), cut, kRed, 20, 0.35).release(), "p");
  mg->Draw("ap");
}

void summaryPage(TTree* modes, int ib, const std::string& name, const std::string& pdf)
{
  const int n = modes->Draw("metricB:metricBErr:explained:metricA", Form("branch==%d", ib), "goff");
  std::vector<double> vB(modes->GetV1(), modes->GetV1() + n), vErr(modes->GetV2(), modes->GetV2() + n),
    vExp(modes->GetV3(), modes->GetV3() + n), vA(modes->GetV4(), modes->GetV4() + n);
  std::vector<std::string> labels;
  {
    std::string mode, *pm = &mode;
    int b = -1;
    modes->SetBranchAddress("mode", &pm);
    modes->SetBranchAddress("branch", &b);
    for (Long64_t i = 0; i < modes->GetEntries(); ++i) {
      modes->GetEntry(i);
      if (b == ib) {
        labels.push_back(mode);
      }
    }
    modes->ResetBranchAddresses();
  }
  auto* c = new TCanvas(Form("cs%d", ib), "summary", 1200, 800);
  c->Divide(2, 1);
  auto* hB = new TH1F(Form("hB%d", ib), Form("%s: mode amplitude;;#mum", name.c_str()), n, 0, n);
  auto* hA = new TH1F(Form("hA%d", ib), "", n, 0, n);
  auto* hE = new TH1F(Form("hE%d", ib), Form("%s: fraction of displacement power explained;;", name.c_str()), n, 0, n);
  for (int i = 0; i < n; ++i) {
    hB->SetBinContent(i + 1, vB[i]);
    hB->SetBinError(i + 1, vErr[i]);
    hA->SetBinContent(i + 1, vA[i]);
    hE->SetBinContent(i + 1, vExp[i]);
    hB->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hE->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  }
  c->cd(1);
  hB->SetStats(0);
  hB->SetFillColorAlpha(kAzure + 1, 0.5);
  hB->Draw("bar e");
  hA->SetLineColor(kRed);
  hA->SetLineWidth(2);
  hA->Draw("hist same");
  auto* leg = new TLegend(0.15, 0.78, 0.55, 0.88);
  leg->AddEntry(hB, "B: fitted displacement field", "f");
  leg->AddEntry(hA, "A: constraint-equation projection", "l");
  leg->Draw();
  c->cd(2);
  hE->SetStats(0);
  hE->SetFillColor(kGreen + 2);
  hE->Draw("bar");
  c->Print(pdf.c_str());
}

void profilePage(TTree* sens, int ib, const std::string& name, const std::string& pdf)
{
  const std::string cut = Form("branch==%d", ib);
  auto* c = new TCanvas(Form("cp%d", ib), "profiles", 1500, 900);
  c->Divide(3, 2);
  const std::string rigR = "(dRigid[0]*x+dRigid[1]*y)/r", rigP = "(-dRigid[0]*y+dRigid[1]*x)/r";
  c->cd(1);
  panel(sens, cut, Form("%s: radialScaling (#Deltar vs r)", name.c_str()), "r [cm]", "r", "dr", rigR + "+" + radialOf("radialScaling"));
  c->cd(2);
  panel(sens, cut, "zScaling (#Deltaz vs z)", "z [cm]", "z", "dz", "dRigid[2]+" + zOf("zScaling"));
  c->cd(3);
  panel(sens, cut, "curl (r#Delta#phi vs r)", "r [cm]", "r", "rdphi", rigP + "+" + phiOf("curl"));
  c->cd(4);
  panel(sens, cut, "twist (r#Delta#phi vs z)", "z [cm]", "z", "rdphi", rigP + "+" + phiOf("twist"));
  c->cd(5);
  panel(sens, cut, "elliptical (#Deltar vs #phi)", "#phi [rad]", "phi", "dr",
        rigR + "+" + radialOf("elliptical_cos2phi") + "+" + radialOf("elliptical_sin2phi"));
  c->cd(6);
  panel(sens, cut, "bowing (#Deltar vs z)", "z [cm]", "z", "dr", rigR + "+" + radialOf("bowing"));
  c->Print(pdf.c_str());
}

// The azimuthal-harmonic families, l = 1..GlobalHarmonicLMax: the modes measured to dominate the ITS
// corrections once the radial scaling and the curl are constrained away. Left column: the whole family
// summed; right column: each order on its own, so the dominant l is visible.
void harmonicPage(TTree* sens, int ib, const std::string& name, const std::string& pdf)
{
  const std::string cut = Form("branch==%d", ib);
  const std::string rigR = "(dRigid[0]*x+dRigid[1]*y)/r", rigP = "(-dRigid[0]*y+dRigid[1]*x)/r";
  auto* c = new TCanvas(Form("ch%d", ib), "harmonics", 1500, 900);
  c->Divide(2, 2);
  c->cd(1);
  panel(sens, cut, Form("%s: radialHarmonic, all l (#Deltar vs #phi)", name.c_str()), "#phi [rad]", "phi", "dr",
        rigR + "+" + familySum("radialHarmonic", true));
  c->cd(2);
  panel(sens, cut, "tangentialHarmonic, all l (r#Delta#phi vs #phi)", "#phi [rad]", "phi", "rdphi",
        rigP + "+" + familySum("tangentialHarmonic", false));
  // per-order amplitude of the two families, read straight from the modes tree by the summary page
  c->cd(3);
  panel(sens, cut, "radialHarmonic l=2,3 (#Deltar vs #phi)", "#phi [rad]", "phi", "dr",
        rigR + "+" + radialOf("radialHarmonic2_cos2phi") + "+" + radialOf("radialHarmonic2_sin2phi") + "+" +
          radialOf("radialHarmonic3_cos3phi") + "+" + radialOf("radialHarmonic3_sin3phi"));
  c->cd(4);
  panel(sens, cut, "tangentialHarmonic l=2 (r#Delta#phi vs #phi)", "#phi [rad]", "phi", "rdphi",
        rigP + "+" + phiOf("tangentialHarmonic2_cos2phi") + "+" + phiOf("tangentialHarmonic2_sin2phi"));
  c->Print(pdf.c_str());
}

// The polar families, fields r^p*(z/r)^k: plotted vs t = z/r, the variable they depend on, one panel per
// displacement component with the whole family summed
void polarPage(TTree* sens, int ib, const std::string& name, const std::string& pdf)
{
  if (colOf("radialPolar_r0t1") < 0) { // written before the polar families existed
    return;
  }
  const std::string cut = Form("branch==%d", ib);
  const std::string rigR = "(dRigid[0]*x+dRigid[1]*y)/r", rigP = "(-dRigid[0]*y+dRigid[1]*x)/r";
  auto* c = new TCanvas(Form("cpol%d", ib), "polar", 1500, 900);
  c->Divide(3, 1);
  c->cd(1);
  panel(sens, cut, Form("%s: radialPolar (#Deltar vs z/r)", name.c_str()), "z/r", "z/r", "dr", rigR + "+" + familySum("radialPolar", 0));
  c->cd(2);
  panel(sens, cut, "tangentialPolar (r#Delta#phi vs z/r)", "z/r", "z/r", "rdphi", rigP + "+" + familySum("tangentialPolar", 1));
  c->cd(3);
  panel(sens, cut, "longitudinalPolar (#Deltaz vs z/r)", "z/r", "z/r", "dz", "dRigid[2]+" + familySum("longitudinalPolar", 2));
  c->Print(pdf.c_str());
}

// #Deltar map in (z, #phi): any structure the modes do not describe shows up here
void mapPage(TTree* sens, int ib, const std::string& name, const std::string& pdf)
{
  auto* c = new TCanvas(Form("cm%d", ib), "map", 1200, 800);
  c->Divide(2, 1);
  c->cd(1);
  auto* p = new TProfile2D(Form("pm%d", ib), Form("%s: #Deltar [#mum];z [cm];#phi [rad]", name.c_str()), 40, -40, 40, 36, -3.1416, 3.1416);
  sens->Draw(Form("phi:z:1e4*dr>>pm%d", ib), Form("branch==%d", ib), "colz");
  p->SetStats(0);
  c->cd(2);
  auto* q = new TProfile2D(Form("pq%d", ib), Form("%s: r#Delta#phi [#mum];z [cm];#phi [rad]", name.c_str()), 40, -40, 40, 36, -3.1416, 3.1416);
  sens->Draw(Form("phi:z:1e4*rdphi>>pq%d", ib), Form("branch==%d", ib), "colz");
  q->SetStats(0);
  c->Print(pdf.c_str());
}
} // namespace

void drawWeakModes(const char* inFile = "weakModes.root", const char* outPdf = "weakModes.pdf")
{
  std::unique_ptr<TFile> f(TFile::Open(inFile));
  if (!f || f->IsZombie()) {
    printf("Cannot open %s\n", inFile);
    return;
  }
  auto* sens = f->Get<TTree>("sensors");
  auto* modes = f->Get<TTree>("modes");
  if (!sens || !modes) {
    printf("%s lacks the sensors/modes trees\n", inFile);
    return;
  }
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  modeCols(modes); // the uMode column layout comes from the file, never hard-coded
  const std::string pdf = outPdf;
  const int nBranch = static_cast<int>(modes->GetMaximum("branch")) + 1;
  bool first = true;
  for (int ib = 0; ib < nBranch; ++ib) {
    std::string name = "branch" + std::to_string(ib);
    {
      std::string det, *pd = &det;
      int b = -1;
      modes->SetBranchAddress("detector", &pd);
      modes->SetBranchAddress("branch", &b);
      for (Long64_t i = 0; i < modes->GetEntries() && b != ib; ++i) {
        modes->GetEntry(i);
      }
      name = det;
      modes->ResetBranchAddresses();
    }
    if (first) {
      TCanvas dummy("d", "d");
      dummy.Print((pdf + "[").c_str());
      first = false;
    }
    summaryPage(modes, ib, name, pdf);
    profilePage(sens, ib, name, pdf);
    harmonicPage(sens, ib, name, pdf);
    polarPage(sens, ib, name, pdf);
    mapPage(sens, ib, name, pdf);
  }
  if (!first) {
    TCanvas dummy("d", "d");
    dummy.Print((pdf + "]").c_str());
  }
}
