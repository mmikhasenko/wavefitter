// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "MIsobarPiPiS.h"
#include "MIsobarChannel.h"
#include "MTwoBodyChannel.h"
#include "mstructures.h"

#include "MmatrixK.h"
#include "MParKeeper.h"

int main(int argc, char *argv[]) {

  // stable
  MIsobarPiPiS pipiS;
  // MIsobar pipiS(1.1, 0.7, 0.14, 0.14, 1, 5.);
  MIsobarChannel ch_q(pipiS, 0.14, 0, 5.);
  MTwoBodyChannel ch_s(1., 0.14, 0, 5.);
  ch_q.makeLookupTable(0.0, 3.9, 100);
  TCanvas c1("c1");
  combine(
          SET1(
               draw([&](double e)->double{return ch_s.rho(e*e);}, 0.1, 4.5),
               SetLineColor(kBlue) ),
          SET1(
               draw([&](double e)->double{return ch_q.rho(e*e);}, 0.1, 4.5),
               SetLineColor(kGreen) )
          )->Draw("al");
  c1.SaveAs("/tmp/test_MIsobarChannelPiPis.pdf(");

  // disperse integral
  ch_q.makeDisperseLookupTable(0.0, 20.0, 100);
  ch_s.makeDisperseLookupTable(0.0, 20.0, 100);
  
  TH2D h2("h2","Complex sheet of  #tilde{#rho}(s);Re@s;Im@s", 20, 0.1, POW2(2.2), 20 ,-0.4, 0.4);
  for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
    for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
      cd s(h2.GetXaxis()->GetBinCenter(i+1), h2.GetYaxis()->GetBinCenter(j+1));
      cd rholtilde = pipiS.U(s);// ch_q.rho(s);
      h2.SetBinContent(i+1, j+1, real(rholtilde));
    }
  h2.SetStats(kFALSE);
  h2.Draw("lego2");
  c1.SaveAs("/tmp/test_MIsobarChannelPiPis.pdf");

  MmatrixK km({&ch_q}, 0);
  km.addPole("m0", "g");
  MParKeeper::gI()->set("m0", 1.322);
  MParKeeper::gI()->set("g0", 5.5);
  
  for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
    for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
      cd s(h2.GetXaxis()->GetBinCenter(i+1), h2.GetYaxis()->GetBinCenter(j+1));
      cd value = km.getK(s)(0,0);
      h2.SetBinContent(i+1, j+1, abs(value));
   }
  h2.SetTitle("Re@Kmatrix");
  h2.Draw("colz");
  c1.SaveAs("/tmp/test_MIsobarChannelPiPis.pdf");

  for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
    for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
      cd s(h2.GetXaxis()->GetBinCenter(i+1), h2.GetYaxis()->GetBinCenter(j+1));
      cd value = km.getFSdenominator(s)(0,0);
      h2.SetBinContent(i+1, j+1, imag(value));
    }
  h2.SetTitle("Im[I - i #rho K]_{I}");
  h2.Draw("colz");
  c1.SaveAs("/tmp/test_MIsobarChannelPiPis.pdf");

  for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
    for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
      cd s(h2.GetXaxis()->GetBinCenter(i+1), h2.GetYaxis()->GetBinCenter(j+1));
      cd value = km.getSSdenominator(s)(0,0);
      h2.SetBinContent(i+1, j+1, imag(value));
    }
  h2.SetTitle("Im[I - i #rho K]_{II}");
  h2.Draw("colz");
  c1.SaveAs("/tmp/test_MIsobarChannelPiPis.pdf");

  for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
    for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
      cd s(h2.GetXaxis()->GetBinCenter(i+1), h2.GetYaxis()->GetBinCenter(j+1));
      cd value = (imag(s) > 0) ? km.getFSdenominator(s)(0,0) :  km.getSSdenominator(s)(0,0);
      h2.SetBinContent(i+1, j+1, imag(value));
    }
  h2.SetTitle("Im[I - i #rho K]_{II}");
  h2.Draw("colz");
  c1.SaveAs("/tmp/test_MIsobarChannelPiPis.pdf)");

  return 0;
}
