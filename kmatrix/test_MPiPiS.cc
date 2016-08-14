// Copyright [2016] Mikhail Mikhasenko

#include <complex>

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH2D.h"
// #include "Math/SpecFuncMathMore.h"

#include "MIsobarPiPiS.h"
// #include "MAscoli.h"
// #include "MmatrixK.h"
// #include "MProductionPhysics.h"
// #include "MParKeeper.h"
#include "mstructures.h"
#include "deflib.h"
#include "constants.h"

int main () {
  
  MIsobarPiPiS pipiS;
  std::cout << pipiS.GetM() << "\n";

  TCanvas c1("c1");
  combine(
          SET3(
               draw([&](double m)->double{return real(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part"),
               SetFillColor(0),
               SetLineColor(kBlack)),
          SET3(
               draw([&](double m)->double{return imag(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part"),
               SetFillColor(0),
               SetLineColor(kRed) ),
          SET3(
               draw([&](double m)->double{return abs(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Abs part"),
               SetFillColor(0),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_PiPiS.pdf(");

  // Phase vs M-matrix
  combine(
          SET3(
               draw([&](double m)->double{return imag(
                                                      2./RHO(m*m, POW2(PI_MASS), POW2(PI_MASS))*
                                                      1./(1./tan(arg(pipiS.T(m*m))) - cd(0., 1) )
                                                      );}, 2*PI_MASS+1.e-3, 0.4),
               SetTitle("Imag part of amplitude from phase"),
               SetFillColor(0),
               SetLineColor(kRed) ),
          SET3(
               draw([&](double m)->double{return real(
                                                      2./RHO(m*m, POW2(PI_MASS), POW2(PI_MASS))*
                                                      1./(1./tan(arg(pipiS.T(m*m))) - cd(0., 1) )
                                                      );}, 2*PI_MASS+1.e-3, 0.4),
               SetTitle("Real part of amplitude from phase"),
               SetFillColor(0),
               SetLineColor(kBlack) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_PiPiS.pdf");

  // Phase vs M-matrix
  combine(
          SET3(
               draw([&](double m)->double{return real(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part"),
               SetFillColor(0),
               SetLineColor(kBlack)),
          SET3(
               draw([&](double m)->double{return imag(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part"),
               SetFillColor(0),
               SetLineColor(kRed) ),
          SET4(
               draw([&](double m)->double{return imag(
                                                      2./RHO(m*m, POW2(PI_MASS), POW2(PI_MASS))*
                                                      1./(1./tan(arg(pipiS.T(m*m))) - cd(0., 1) )
                                                      );}, 2*PI_MASS+1.e-3, 2.2),
               SetTitle("Imag part of amplitude from phase"),
               SetFillColor(0),
               SetLineStyle(2),
               SetLineColor(kRed) ),
          SET4(
               draw([&](double m)->double{return real(
                                                      2./RHO(m*m, POW2(PI_MASS), POW2(PI_MASS))*
                                                      1./(1./tan(arg(pipiS.T(m*m))) - cd(0., 1) )
                                                      );}, 2*PI_MASS+1.e-3, 2.2),
               SetTitle("Real part of amplitude from phase"),
               SetFillColor(0),
               SetLineStyle(2),
               SetLineColor(kBlack) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_PiPiS.pdf(");

  // 
  SET3(
       draw([&](double m)->double{return 180./M_PI*arg(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
       SetTitle("Phase shift"),
       GetXaxis()->SetTitle("GeV"),
       SetLineColor(kBlack) )->Draw("al");
  c1.SaveAs("/tmp/test_PiPiS.pdf");

  /* Plot: real part T(Real s) and real part of T(s+i0) */
  combine(
          SET3(
               draw([&](double m)->double{return real(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetFillColor(0),
               SetTitle("Real part [real s]"),
               SetLineColor(kBlack)),
          SET4(
               draw([&](double m)->double{return real(pipiS.T(cd(m*m,1e-3)));}, 2*PI_MASS, 2.2),
               SetTitle("Real part [s + i0]"),
               SetLineColor(kBlack),
               SetFillColor(0),
               SetLineStyle(2) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_PiPiS.pdf");

  /* Plot: |T|^2*RHO vs U(s+i0) */
  combine(
          SET3(
               draw([&](double m)->double{return norm(pipiS.T(m*m)) *
                     1./(8*M_PI)*sqrt(LAMBDA(m*m, POW2(PI_MASS), POW2(PI_MASS)))/(m*m);
                 }, 2*PI_MASS, 2.2),
               SetTitle("Real part [real s]"),
               SetFillColor(0),
               SetLineColor(kBlack)),
          SET4(
               draw([&](double m)->double{return real(pipiS.U(cd(m*m, 1e-3)));}, 2*PI_MASS, 2.2),
               SetTitle("Real part [s + i0]"),
               SetFillColor(0),
               SetLineColor(kBlack),
               SetLineStyle(2) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_PiPiS.pdf");

  TH2D h2("h2","Complex sheet of U(s);Re@s;Im@s", 100, 0.1, POW2(2.2), 100 ,-0.4, 0.4);
  for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
    for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
      cd s(h2.GetXaxis()->GetBinCenter(i+1), h2.GetYaxis()->GetBinCenter(j+1));
      double U_II = abs(pipiS.U(s));
      h2.SetBinContent(i+1, j+1, U_II);
    }
  h2.SetStats(kFALSE);
  // h2.GetZaxis()->SetRangeUser(0,30.);
  h2.Draw("colz");
  c1.SaveAs("/tmp/test_PiPiS.pdf");

  /* compare to BW */
  MIsobar pipiS_bw(0.5, 0.5,  PI_MASS, PI_MASS, 0);

  combine(
          SET3(
               draw([&](double m)->double{return real(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part AMP"),
               SetFillColor(0),
               SetLineColor(kBlack)),
          SET3(
               draw([&](double m)->double{return imag(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part AMP"),
               SetFillColor(0),
               SetLineColor(kRed) ),
          SET3(
               draw([&](double m)->double{return abs(pipiS.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Abs part AMP"),
               SetFillColor(0),
               SetLineColor(kGreen) ),
          SET4(
               draw([&](double m)->double{return real(pipiS_bw.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part BW"),
               SetFillColor(0),
               SetLineStyle(2),
               SetLineColor(kBlack)),
          SET4(
               draw([&](double m)->double{return imag(pipiS_bw.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part BW"),
               SetFillColor(0),
               SetLineStyle(2),
               SetLineColor(kRed) ),
          SET4(
               draw([&](double m)->double{return abs(pipiS_bw.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Abs part BW"),
               SetFillColor(0),
               SetLineStyle(2),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_PiPiS.pdf");

  for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
    for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
      cd s(h2.GetXaxis()->GetBinCenter(i+1), h2.GetYaxis()->GetBinCenter(j+1));
      double U_II = abs(pipiS_bw.U(s));
      h2.SetBinContent(i+1, j+1, U_II);
    }
  h2.SetStats(kFALSE);
  // h2.GetZaxis()->SetRangeUser(0,30.);
  h2.Draw("colz");
  c1.SaveAs("/tmp/test_PiPiS.pdf)");
  
  return 0;
}
