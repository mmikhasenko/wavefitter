// Copyright [2016] Mikhail Mikhasenko

#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MIsobarChannel.h"
#include "MIsobarPiPiS.h"
#include "MTwoBodyChannel.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {
  
  MIsobarPiPiS  pipiS_iso;
  MIsobar rho_iso(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
  MIsobar  f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.);

  TCanvas c1("c1");
  TMultiGraph *m;  
  /* RHO */
  m = combine(
          SET3(
               draw([&](double m)->double{return real(pipiS_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part"),
               SetFillColor(0),
               SetLineColor(kBlue)),
          SET3(
               draw([&](double m)->double{return imag(pipiS_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part"),
               SetFillColor(0),
               SetLineColor(kOrange) ),
          SET3(
               draw([&](double m)->double{return abs(pipiS_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Abs part"),
               SetFillColor(0),
               SetLineColor(kGreen) ) );
  m->SetTitle("Amplitude pipiS(s)");
  m->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_isobar_shale.pdf(");
  c1.Clear();

  /* RHO */
  m = combine(
          SET3(
               draw([&](double m)->double{return real(rho_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part"),
               SetFillColor(0),
               SetLineColor(kBlue)),
          SET3(
               draw([&](double m)->double{return imag(rho_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part"),
               SetFillColor(0),
               SetLineColor(kOrange) ),
          SET3(
               draw([&](double m)->double{return abs(rho_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Abs part"),
               SetFillColor(0),
               SetLineColor(kGreen) ) );
  m->SetTitle("Amplitude rho(s)");
  m->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_isobar_shale.pdf");
  c1.Clear();

  /* F2 */
  m = combine(
          SET3(
               draw([&](double m)->double{return real(f2_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part"),
               SetFillColor(0),
               SetLineColor(kBlue)),
          SET3(
               draw([&](double m)->double{return imag(f2_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part"),
               SetFillColor(0),
               SetLineColor(kOrange) ),
          SET3(
               draw([&](double m)->double{return abs(f2_iso.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Abs part"),
               SetFillColor(0),
               SetLineColor(kGreen) ) );
  m->SetTitle("Amplitude f2(s)");
  m->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_isobar_shale.pdf)");
  c1.Clear();

  return 0;
}

