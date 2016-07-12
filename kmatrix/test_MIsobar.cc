// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MIsobarChannel.h"
#include "MTwoBodyChannel.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {
  // stable
  MIsobar rho(0.77, 0.15, 0.14, 0.14, 1, 5.);
  MIsobarChannel rho_q(rho, 0.14, 2, 5.);
  MTwoBodyChannel rho_s(rho.GetM(), 0.14, 2, 5.);
  rho_q.makeLookupTable(0.0, 6.0, 100);

  TCanvas c1("c1");
  combine(
          SET1(
               draw([&](double e)->double{return rho_s.rho(e*e);}, 0.1, 2.5),
               SetLineColor(kBlue) ),
          SET1(
               draw([&](double e)->double{return rho_q.rho(e*e);}, 0.1, 2.5),
               SetLineColor(kGreen) )
          )
    ->Draw("al");
  c1.SaveAs("/tmp/a.pdf");

  // disperse integral
  rho_q.makeDisperseLookupTable(0.0, 10.0, 100);
  rho_s.makeDisperseLookupTable(0.0, 10.0, 100);
  combine(
          SET1(
               draw([&](double e)->double{return real(rho_q.rholtilde(e*e));}, 0.001, 2.5, 100),
               SetLineColor(kBlack) ),
          SET1(
               draw([&](double e)->double{return imag(rho_q.rholtilde(e*e));}, 0.001, 2.5, 100),
               SetLineColor(kRed) ),
          SET2(
               draw([&](double e)->double{return real(rho_s.rholtilde(e*e));}, 0.001, 2.5, 100),
               SetLineColor(kBlack), SetLineStyle(kDashed) ),
          SET2(
               draw([&](double e)->double{return imag(rho_s.rholtilde(e*e));}, 0.001, 2.5, 100),
               SetLineColor(kRed), SetLineStyle(kDashed) )
               )
    ->Draw("al");
  c1.SaveAs("/tmp/b.pdf");

  return 0;
}
