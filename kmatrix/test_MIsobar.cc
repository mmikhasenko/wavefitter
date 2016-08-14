// Copyright [2016] Mikhail Mikhasenko

#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MIsobarChannel.h"
#include "MTwoBodyChannel.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {
  // stable
  MIsobar rho(0.77, 0.15, PI_MASS, PI_MASS, 1, 5.);
  std::cout << rho.GetM() << "\n";

  TCanvas c1("c1");
  combine(
          SET2(
               draw([&](double m)->double{return real(rho.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part"),
               SetLineColor(kBlack)),
          SET2(
               draw([&](double m)->double{return imag(rho.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part"),
               SetLineColor(kRed) ),
          SET2(
               draw([&](double m)->double{return abs(rho.T(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Abs part"),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.SaveAs("/tmp/test_MIsobar.pdf(");

  SET2(
       draw([&](double m)->double{return 180./M_PI*arg(rho.T(m*m));}, 2*PI_MASS, 2.2),
       SetTitle("Phase shift"),
       SetLineColor(kBlack) )->Draw("al");
  c1.SaveAs("/tmp/test_MIsobar.pdf");

   TH2D h2("h2","Second sheet of U(s)", 100, 0.1, 1, 100 ,-0.5, 0.5);
   for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
     for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
       cd s(h2.GetXaxis()->GetBinCenter(i+1), h2.GetYaxis()->GetBinCenter(j+1));
       double U_II = abs(rho.U(s));
       h2.SetBinContent(i+1, j+1, U_II);
     }
 
   h2.GetZaxis()->SetRangeUser(0,30.); h2.Draw("lego2");
   c1.SaveAs("/tmp/test_MIsobar.pdf)");
    
  return 0;
}


  
//  TCanvas c1("c1");
//  combine(
//          SET1(
//               draw([&](double e)->double{return rho_s.rho(e*e);}, 0.1, 2.5),
//               SetLineColor(kBlue) ),
//          SET1(
//               draw([&](double e)->double{return rho_q.rho(e*e);}, 0.1, 2.5),
//               SetLineColor(kGreen) )
//          )
//    ->Draw("al");
//  c1.SaveAs("/tmp/test_MIsobar.pdf(");
//
//  // disperse integral
//  rho_q.makeDisperseLookupTable(0.0, 10.0, 100);
//  rho_s.makeDisperseLookupTable(0.0, 10.0, 100);
//  combine(
//          SET1(
//               draw([&](double e)->double{return real(rho_q.rholtilde(e*e));}, 0.001, 2.5, 100),
//               SetLineColor(kBlack) ),
//          SET1(
//               draw([&](double e)->double{return imag(rho_q.rholtilde(e*e));}, 0.001, 2.5, 100),
//               SetLineColor(kRed) ),
//          SET2(
//               draw([&](double e)->double{return real(rho_s.rholtilde(e*e));}, 0.001, 2.5, 100),
//               SetLineColor(kBlack), SetLineStyle(kDashed) ),
//          SET2(
//               draw([&](double e)->double{return imag(rho_s.rholtilde(e*e));}, 0.001, 2.5, 100),
//               SetLineColor(kRed), SetLineStyle(kDashed) )
//               )
//    ->Draw("al");
//  c1.SaveAs("/tmp/test_MIsobar.pdf)");
