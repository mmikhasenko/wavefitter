// Copyright [2017] Misha Mikhasenko
// the program plots isobars used for the COMPASS PW expansion

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"

#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "waves.h"

#include "constants.h"
#include "mstructures.h"
#include "deflib.h"

int main(int ac, char **av) {

  {
  // Create isobars
  MIsobar rho(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.); rho.setIntU();
  MIsobar  f2(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.); f2.setIntU();
  MIsobar rho3(1.69, 0.16, PI_MASS, PI_MASS, 1, 5.); rho3.setIntU();
  MIsobarPiPiS pipiS; pipiS.setIntU();
  MIsobar *iso[] = {&pipiS, &rho, &f2, &rho3};
  MIsobar f980(0.99, 0.04, PI_MASS, PI_MASS, 0); f980.setIntU();
  MIsobar f1500(1.504, 0.11, PI_MASS, PI_MASS, 0); f1500.setIntU();
  MIsobar *iso_scalars[] = {&pipiS, &f980, &f1500};

//  const uint Npoints = 1000;
//  std::pair<double, double> ranges(make_pair(0.3, 2.3));
//  std::vector<double> x(Npoints), y(Npoints);
//  for (uint i = 0; i < Npoints; i++)
//    x[i] = ranges.first + (ranges.second - ranges.first)/(Npoints-1)*i;
  TCanvas c1("c1");
  TMultiGraph *m =
    combine(
          SET3(
               draw([&](double m)->double{return abs(iso_scalars[0]->ToneVertex(m*m));}, 0.3, 2.3, 250),
               SetTitle("(#pi#pi)_{S}"),
               SetFillColor(0),
               SetLineColor(kBlack)),
          SET3(
               draw([&](double m)->double{return abs(iso_scalars[1]->ToneVertex(m*m));}, 0.3, 2.3, 250),
               SetTitle("f_{0}(980)"),
               SetFillColor(0),
               SetLineColor(kMagenta) ),
          SET3(
               draw([&](double m)->double{return abs(iso_scalars[2]->ToneVertex(m*m));}, 0.3, 2.3, 250),
               SetTitle("f_{0}(1500)"),
               SetFillColor(0),
               SetLineColor(kAzure) ) );
  m->Draw("al");
  m->SetTitle("Scalar resonances");
  m->GetXaxis()->SetTitle("M_{#pi#pi}");
  m->GetYaxis()->SetTitle("absolute value |T_{l}|");
  c1.BuildLegend(0.6, 0.75, 0.9, 0.9);
  c1.Print("/tmp/test_compass_isobars.pdf(");

  // Second plot
  m = combine(
          SET3(
               draw([&](double m)->double{return abs(iso[1]->ToneVertex(m*m));}, 0.3, 2.3, 250),
               SetTitle("#rho(770)"),
               SetFillColor(0),
               SetLineColor(kBlack)),
          SET3(
               draw([&](double m)->double{return abs(iso[2]->ToneVertex(m*m));}, 0.3, 2.3, 250),
               SetTitle("f_{2}(1270)"),
               SetFillColor(0),
               SetLineColor(kMagenta) ),
          SET3(
               draw([&](double m)->double{return abs(iso[3]->ToneVertex(m*m));}, 0.3, 2.3, 250),
               SetTitle("#rho_{3}(1690)"),
               SetFillColor(0),
               SetLineColor(kAzure) ) );
  m->Draw("al");
  m->SetTitle("Non-scalar resonances");
  m->GetXaxis()->SetTitle("M_{#pi#pi}");
  m->GetYaxis()->SetTitle("absolute value |T_{l}|");
  c1.BuildLegend(0.6, 0.75, 0.9, 0.9);
  c1.Print("/tmp/test_compass_isobars.pdf)");
  }

  {
  MIsobar rho(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5);
  MIsobar f2(F2_MASS, F2_WIDTH, PI_MASS, PI_MASS, 2, 5);
  std::function<cd(double)> isobarT[3];
  isobarT[0] = [&](double s) -> cd { return waves::GKPY::T(s); };
  isobarT[1] = [&](double s) -> cd { return rho.T(s); };
  isobarT[2] = [&](double s) -> cd { return f2.T(s); };

  TCanvas c1("c1");
TMultiGraph *m;
  for (uint i = 0; i < 3; i++) {
  m =
    combine(
          SET3(
               draw([&](double m)->double{return real(isobarT[i](m*m));}, 0.3, 2.3, 250),
                 SetTitle(TString::Format("Re t_{%d}", i)),
               SetFillColor(0),
               SetLineColor(kBlack)),
          SET4(
               draw([&](double m)->double{return imag(isobarT[i](m*m));}, 0.3, 2.3, 250),
               SetTitle(TString::Format("Im t_{%d}", i)),
               SetFillColor(0),
               SetLineStyle(kDashed),
               SetLineColor(kBlack) ) );
  m->Draw("al");
  m->SetTitle("Elastic #pi#pi scattering amplitude");
  m->GetXaxis()->SetTitle("M_{#pi#pi}");
  c1.BuildLegend(0.65, 0.7, 0.85, 0.85);
  if (i == 0) c1.Print("/tmp/test_deck_isobars.pdf(");
  if (i == 1) c1.Print("/tmp/test_deck_isobars.pdf");
  if (i == 2) c1.Print("/tmp/test_deck_isobars.pdf)");

  }
}

  return 0;
}
