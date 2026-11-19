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

// Lab reference followed by the axes of the frames in which O2Align defines the rigid-body DOFs, parent -> child.
// For a sensor that is the local frame its Sensor* class builds in defineMatrixL2G (a TRD chamber: the TGeo chamber
// frame rotated to X radial, Y rphi, Z along the beam; an ITS chip: the embedded sensor layer), otherwise the TGeo
// matrix of the physical node.
// Call as drawLocalFrames(<geomName>, std::vector<std::string> volNames, <azimuth = 35>, <elevation = 25>)
// to draw the frames as seen in the lab for the whole hierarchy of the volume. The Sensor* classes come from
// libO2AlignRS, which must be loaded before the macro is compiled:
//   root -l -e 'gSystem->Load("libO2AlignRS")' -e '.L drawLocalFrames.C+' -e 'drawLocalFrames(...)'
// E.g drawLocalFrames("o2sim_geometry-aligned.root",
//                     {"ITS/ITSULayer6/ITSUHalfBarrel0/ITSUStave6/ITSUHalfStave0/ITSUModule0/ITSUChip4",
//                      "ITS/ITSULayer6/ITSUHalfBarrel0/ITSUStave0/ITSUHalfStave0/ITSUModule0/ITSUChip4",
//                      "ITS/ITSULayer0/ITSUHalfBarrel0/ITSUStave3/ITSUChip4",
//                      "TRD/sm02/st3/pl0","TOF/sm03/strip35"})
// draws the doc/local_frames_ITS.pdf; doc/local_frames_TRD.pdf has TRD/sm00/st1/pl1, TRD/sm02/st4/pl1, TRD/sm05/st1/pl1,
// TRD/sm13/st4/pl1, doc/local_frames_TOF.pdf TOF/sm00/strip01, strip25, strip45, strip66, strip90.
//
// At zero angles: lab Z right, lab Y up, lab X into the screen.
// Azimuth rotates about lab Y from this view; elevation tilts the camera up.
// Both angles are in degrees.
// Display origins are separated for readability; physical translations are not drawn.

#include <TArrow.h>
#include <TCanvas.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TString.h>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <Framework/Logger.h>
#include <DetectorsBase/GeometryManager.h>
#include <ITSBase/GeometryTGeo.h>
#include "O2Align/Label.h"
#include "O2Align/SensorITS.h"
#include "O2Align/SensorTOF.h"
#include "O2Align/SensorTRD.h"

namespace
{
/// local-to-lab matrix built by the sensor class S for the alignable sym (sensorId: the ITS chip index)
template <class S>
TGeoHMatrix sensorL2G(const std::string& sym, int sensorId = -1)
{
  S sensor(sym.c_str(), o2::alignrs::Label(0u, 0u, true));
  if (sensorId >= 0) {
    sensor.setSensorId(sensorId);
  }
  sensor.defineMatrixL2G();
  return sensor.getMatrixL2G();
}

/// ITS chip index of the chip alignable sym, as assigned in DetectorITS::buildHierarchy; -1 if not found
int itsChipIndex(const std::string& sym)
{
  auto* geom = o2::its::GeometryTGeo::Instance();
  int lay, hba, sta, ssta, mod, chip;
  for (int ich = 0; ich < geom->getNumberOfChips(); ich++) {
    geom->getChipId(ich, lay, hba, sta, ssta, mod, chip);
    const int nch = mod < 0 ? geom->getNumberOfChipsPerStave(lay) : geom->getNumberOfChipsPerModule(lay);
    if (sym == geom->composeSymNameChip(lay, hba, sta, ssta, mod, chip % nch)) {
      return ich;
    }
  }
  return -1;
}

/// local-to-lab matrix of the frame of the rigid-body DOFs of the alignable sym, whose physical-node matrix is pnMatrix
TGeoHMatrix dofFrameL2G(const std::string& sym, const TGeoHMatrix& pnMatrix)
{
  static const std::regex trdChamber(R"(TRD/sm\d+/st\d/pl\d)"), tofStrip(R"(TOF/sm\d+/strip\d+)"), itsChip(R"(ITS/.*ITSUChip\d+)");
  if (std::regex_match(sym, trdChamber)) {
    return sensorL2G<o2::alignrs::SensorTRD>(sym);
  }
  if (std::regex_match(sym, tofStrip)) {
    return sensorL2G<o2::alignrs::SensorTOF>(sym);
  }
  if (std::regex_match(sym, itsChip)) {
    const int ich = itsChipIndex(sym);
    if (ich >= 0) {
      return sensorL2G<o2::alignrs::SensorITS>(sym, ich);
    }
    std::cerr << "No ITS chip index for " << sym << ", drawing its physical-node frame\n";
  }
  return pnMatrix;
}
} // namespace

void drawLocalFrames(const std::string& volumeName, double azimuth = 35., double elevation = 25.)
{
  auto* target = gPad;
  std::string volume{volumeName};
  while (!volume.empty() && volume.front() == '/')
    volume.erase(0, 1);
  while (!volume.empty() && volume.back() == '/')
    volume.pop_back();
  if (!gGeoManager->GetAlignableEntry(volume.c_str())) {
    std::cerr << "No alignable entry: " << volume << '\n';
    return;
  }
  struct Frame {
    std::string name;
    std::string label;
    TGeoHMatrix l2g;     // frame of the rigid-body DOFs
    TGeoHMatrix l2gTGeo; // TGeo physical-node frame
    bool showTGeo;
  };
  std::vector<Frame> frames;
  std::size_t labelStart = 0;
  // Symbolic hierarchy, not the intervening mechanical nodes in the TGeo path.
  for (std::size_t end = volume.find('/');; end = volume.find('/', end + 1)) {
    const auto prefix = volume.substr(0, end);
    auto* entry = gGeoManager->GetAlignableEntry(prefix.c_str());
    if (entry) {
      auto* node = entry->GetPhysicalNode();
      if (!node)
        node = gGeoManager->MakePhysicalNode(entry->GetTitle());
      if (!node || !node->GetMatrix()) {
        std::cerr << "Cannot obtain local-to-lab matrix: " << prefix << '\n';
        return;
      }
      frames.push_back({prefix, prefix.substr(labelStart), dofFrameL2G(prefix, *node->GetMatrix()), *node->GetMatrix(), true});
      labelStart = prefix.size() + 1;
    } else {
      std::cerr << "Skipping symbolic prefix without an alignable entry: " << prefix << '\n';
    }
    if (end == std::string::npos)
      break;
  }

  // Identity local-to-lab transformation supplies the lab reference triplet.
  frames.insert(frames.begin(), Frame{"Lab frame", "Lab frame", TGeoHMatrix{}, TGeoHMatrix{}, false});

  const double az = azimuth * TMath::DegToRad();
  const double el = elevation * TMath::DegToRad();
  // Screen right/up; camera position direction is
  // (-cos(el)*cos(az), sin(el), cos(el)*sin(az)). No perspective distortion.
  const std::array<double, 3> right{{std::sin(az), 0., std::cos(az)}};
  const std::array<double, 3> up{{std::sin(el) * std::cos(az), std::cos(el), -std::sin(el) * std::sin(az)}};
  const int colors[3] = {kRed + 1, kGreen + 2, kBlue + 1};
  const char* labels[3] = {"X", "Y", "Z"};
  static unsigned serial = 0;
  ++serial;
  const double width = target->GetCanvas()->GetWw() * target->GetAbsWNDC();
  const double height = target->GetCanvas()->GetWh() * target->GetAbsHNDC();
  auto* heading = new TLatex;
  heading->SetNDC();
  heading->SetTextFont(42);
  heading->SetTextSize(0.035);
  heading->DrawLatex(0.015, 0.95, "DOFs rigid-body (and local TGeo - dashed) axes in Lab orientation: parent #rightarrow child");
  heading->SetTextSize(0.026);
  heading->DrawLatex(0.015, 0.89, Form("Orthographic lab view: azimuth %.1f#circ, elevation %.1f#circ; origins separated", azimuth, elevation));

  for (std::size_t i = 0; i < frames.size(); ++i) {
    target->cd();
    auto* pad = new TPad(Form("frame_%u_%zu", serial, i), frames[i].name.c_str(),
                         double(i) / frames.size(), 0., double(i + 1) / frames.size(), 0.83);
    pad->Draw();
    pad->cd();
    // Pad dimensions and range give equal horizontal/vertical pixel scales.
    const double aspect = width / (frames.size() * height * 0.83);
    const double halfX = 1.65 * (aspect > 1. ? aspect : 1.);
    pad->Range(-halfX, -halfX / aspect, halfX, halfX / aspect);
    auto* text = new TLatex;
    text->SetTextFont(42);
    text->SetTextAlign(22);
    text->SetTextSize(0.065);
    // Include non-alignable components since the preceding alignable ancestor.
    text->SetNDC();
    text->SetTextSize(0.055);
    text->DrawLatex(0.5, 0.88, frames[i].label.c_str());
    text->SetNDC(false);
    text->SetTextSize(0.065);
    // projection of the local axis of the matrix to the screen; prints its lab components
    auto drawTriplet = [&](const TGeoHMatrix& l2g, double ox, double oy, bool tgeo) {
      std::cout << frames[i].name << (tgeo ? "\n  TGeo" : "\n  DOF") << " local axes in lab (X,Y,Z):\n";
      for (int axis = 0; axis < 3; ++axis) {
        double local[3] = {0., 0., 0.}, lab[3];
        local[axis] = 1.;
        l2g.LocalToMasterVect(local, lab);
        double x = 0., y = 0.;
        for (int k = 0; k < 3; ++k) {
          x += right[k] * lab[k];
          y += up[k] * lab[k];
        }
        std::cout << "  " << labels[axis] << ": (" << lab[0] << ", " << lab[1] << ", " << lab[2] << ")\n";
        const TString name = tgeo ? Form("(%s)", labels[axis]) : labels[axis];
        text->SetTextColor(colors[axis]);
        text->SetTextSize(tgeo ? 0.05 : 0.065);
        if (std::hypot(x, y) < 1.e-8) {
          const double towardsCamera = -std::cos(el) * std::cos(az) * lab[0] +
                                       std::sin(el) * lab[1] + std::cos(el) * std::sin(az) * lab[2];
          text->DrawLatex(ox, oy - 0.22, Form("%s %s", name.Data(), towardsCamera > 0. ? "#odot" : "#otimes"));
          continue;
        }
        auto* arrow = new TArrow(ox, oy, ox + x, oy + y, tgeo ? 0.004 : 0.006, "|>");
        arrow->SetLineColor(colors[axis]);
        arrow->SetFillColor(colors[axis]);
        arrow->SetLineWidth(tgeo ? 1 : 3);
        arrow->SetLineStyle(tgeo ? 2 : 1);
        arrow->Draw();
        text->DrawLatex(ox + 1.18 * x, oy + 1.18 * y, name.Data());
      }
    };
    drawTriplet(frames[i].l2g, 0., 0., false);
    if (frames[i].showTGeo) { // TGeo triplet next to the DOF one
      drawTriplet(frames[i].l2gTGeo, 0.35 * halfX, -0.3 * halfX / aspect, true);
    }
  }
  target->cd();
  target->Modified();
  target->Update();
}

TCanvas* cnv = nullptr;

void drawLocalFrames(const std::string& geomFile, std::vector<std::string> vols, double azimuth = 35., double elevation = 25.)
{
  if (!o2::base::GeometryManager::isGeometryLoaded()) {
    o2::base::GeometryManager::loadGeometry(geomFile);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGP(fatal, "failed to load geometry {}", geomFile);
    }
  }
  int nv = vols.size();
  if (!nv) {
    LOGP(warn, "vector of volumes to draw is empty");
    return;
  }
  cnv = new TCanvas("frames", "frames", 800, 1000);
  cnv->Divide(1, nv);
  int cnvCnt = 0;
  for (auto& vn : vols) {
    cnv->cd(++cnvCnt);
    drawLocalFrames(vn, azimuth, elevation);
  }
}
