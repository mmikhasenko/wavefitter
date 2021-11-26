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

  // uint J = 2;

  double mAsq = POW2(PI_MASS);
  double mBsq = POW2(PROT_MASS);
  double mDsq = POW2(PROT_MASS);
  double mtRsq = POW2(PI_MASS);
  double stot = 2*190*PROT_MASS;
  double t = -0.01;

  // Isobars
  MIsobar rho(0.77, 0.15, 0.14, 0.14, 1, 5.);
  MIsobar  f2(1.23, 0.2,  0.14, 0.14, 2, 5.);
  // constexpr uint Nch = 4;
  // const MIsobar *iso[Nch] = {&rho, &rho, &f2, &f2};
  const MIsobar & ciso = f2;
  // double mR = ciso.GetM();
  // uint cL = 3;  // partial wave

  // Deck channels
  uint lamS = 0;
  double R = 5.;

  TCanvas c1("c1");
  /*  Angular distribution */
  SET3(
       draw([&](double z)->double{
           return MAscoli::upperPart(z, POW2(ciso.GetM()), ciso.GetL(), lamS, R, POW2(2.2),
                                     t, mtRsq, mAsq, POW2(PI_MASS)) *
             2*M_PI *
             MAscoli::sPionProton(z, M_PI/2, POW2(ciso.GetM()), POW2(2.2),
                                  t, stot, mAsq, mBsq, mDsq, POW2(PI_MASS)) *
             sqrt(2*ciso.GetL()+1);
         }, -1, 1, 300),
       SetLineColor(kOrange),
       SetLineWidth(1),
       SetFillStyle(0))->Draw("al");
  c1.SaveAs("/tmp/test_MAscoli.pdf(");

  /* Mass distribution */
  SET3(
       draw([&](double w)->double{
           return MAscoli::upperPart(1., POW2(ciso.GetM()), ciso.GetL(), lamS, R, w*w,
                                     t, mtRsq, mAsq, POW2(PI_MASS)) *
             2*M_PI *
             MAscoli::sPionProton(1., M_PI/2, POW2(ciso.GetM()), w*w,
                                  t, stot, mAsq, mBsq, mDsq, POW2(PI_MASS)) *
             sqrt(2*ciso.GetL()+1);
         }, ciso.GetM()+PI_MASS+0.01, 3.),
       SetLineColor(kRed),
       SetLineWidth(2),
       SetFillStyle(0))->Draw("alp");
  c1.SaveAs("/tmp/test_MAscoli.pdf)");
  
  return 0;
}
