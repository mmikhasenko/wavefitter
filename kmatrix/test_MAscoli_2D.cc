// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLegend.h"
#include "TH2D.h"
#include "Math/SpecFuncMathMore.h"

#include "MIsobar.h"
#include "MAscoli.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.h"
#include "mintegrate.h"
#include "deflib.h"

int main(int argc, char *argv[]) {
  // define all variables
  double mAsq = POW2(PI_MASS);  // beam
  double mBsq = POW2(PROT_MASS);  // target
  double mDsq = POW2(PROT_MASS);  // recoil
  double mtRsq = POW2(PI_MASS);  // exchange
  double m1sq = POW2(PI_MASS);  // bachelor
  double stot = 2*190*PROT_MASS;
  double t = -0.01;
  
  // Isobars
  MIsobar sigma(0.6, 0.15, 0.14, 0.14);
  MIsobar rho(0.77, 0.15, 0.14, 0.14, 1, 5.);
  MIsobar  f2(1.23, 0.2,  0.14, 0.14, 2, 5.);
  // pick current 
  // const MIsobar & ciso = rho;
  const MIsobar & ciso = sigma;

  // Deck channels
  uint lamS = 0;
  double R = 5.;

  TCanvas c1("c1");
  TH2D eVSz("heVSz", "#sqrt{s} vs z for Ascoli Deck;#sqrt{s}(GeV);z",
            100, ciso.GetM()+PI_MASS, 2.5,
            50, -1., 1.);
  for (int i = 0; i < eVSz.GetXaxis()->GetNbins(); i++)
    for (int j = 0; j < eVSz.GetYaxis()->GetNbins(); j++) {
      double e = eVSz.GetXaxis()->GetBinCenter(i+1);
      double z = eVSz.GetYaxis()->GetBinCenter(j+1);
      double f = 2*M_PI*MAscoli::getReducedDeck(z, M_PI/2.,
                                  POW2(ciso.GetM()), ciso.GetL(), lamS, R,
                                  e*e, t,
                                  mtRsq,
                                  stot,
                                  mAsq, mBsq, mDsq,
                                  m1sq);
      eVSz.SetBinContent(i+1, j+1, f);
    }
  eVSz.Draw("colz");
  c1.SaveAs("/tmp/test_MAscoli_2D.pdf(");
  eVSz.Draw("lego");
  c1.SetPhi(-30);
  c1.SaveAs("/tmp/test_MAscoli_2D.pdf)");
  c1.SaveAs("/tmp/test_MAscoli_2D.png");
  //c1.SetPhi(80);
  
  return 0;
}
