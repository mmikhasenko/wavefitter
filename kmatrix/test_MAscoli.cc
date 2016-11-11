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
  MAscoli bch(mAsq, mBsq, POW2(2.2), mDsq, mtRsq, POW2(ciso.GetM()), stot, t,
              ciso.GetL(), lamS, R);

  TCanvas c1("c1");
  /*  Angular distribution */
  combine( 
          SET3(
               draw([&](double z)->double{
                   return bch.getValue(POW2(2.2), z);
                 }, -1, 1, 300),
               SetLineColor(kBlack),
               SetLineWidth(1),
               SetFillStyle(0)),
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
               SetFillStyle(0)))->Draw("al");
  c1.SaveAs("/tmp/test_MAscoli.pdf(");

  /* Mass distribution */
  combine(
          SET3(
               draw([&](double w)->double{
                   return bch.getValue(w*w, 1);
                 }, ciso.GetM()+PI_MASS+0.01, 3.),
               SetLineColor(kBlack),
               SetLineWidth(2),
               SetFillStyle(0)),
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
               SetFillStyle(0)))->Draw("alp");
  c1.SaveAs("/tmp/test_MAscoli.pdf)");



//   /*******************************************************************************************/
//   /*******************************************************************************************/
//   /*******************************************************************************************/
//   SET4(
//        draw([&](double w)->double{
//            double mR = rho.GetM();
//            uint S = rho.GetL();
//            uint L = 1;
//            double m23 = 2*PI_MASS+(mR-2*PI_MASS)*(1.-exp(-1./mR*(w-3*PI_MASS)));
//            return bch.getProjectedDeck(mAsq, mBsq, w*w, mDsq, mtRsq, m23*m23, s, t,
//                                         J, 0, L, S, R);
//          }, 3*PI_MASS, 3.),
//        SetTitle("Smeared Ascoli-Jones \"Deck\" for #rho#pi P-wave.;M_{3#pi}"),
//        SetLineColor(kRed),
//        SetLineWidth(2),
//        SetFillStyle(0))->Draw("alp");
//   c1.SaveAs("/tmp/test_MAscoli.pdf");
//   SET4(
//        draw([&](double w)->double{
//            double mR = rho.GetM();
//            uint S = rho.GetL();
//            uint L = 3;
//            double m23 = 2*PI_MASS+(mR-2*PI_MASS)*(1.-exp(-1./mR*(w-3*PI_MASS)));
//            return bch.getProjectedDeck(mAsq, mBsq, w*w, mDsq, mtRsq, m23*m23, s, t,
//                                         J, 0, L, S, R);
//          }, 3*PI_MASS, 3.),
//        SetTitle("Smeared Ascoli-Jones \"Deck\" for #rho#pi F-wave.;M_{3#pi}"),
//        SetLineColor(kRed),
//        SetLineWidth(2),
//        SetFillStyle(0))->Draw("alp");
//   c1.SaveAs("/tmp/test_MAscoli.pdf");
//   SET4(
//        draw([&](double w)->double{
//            double mR = f2.GetM();
//            uint S = f2.GetL();
//            uint L = 0;
//            double m23 = 2*PI_MASS+(mR-2*PI_MASS)*(1.-exp(-1./mR*(w-3*PI_MASS)));
//            return bch.getProjectedDeck(mAsq, mBsq, w*w, mDsq, mtRsq, m23*m23, s, t,
//                                         J, 0, L, S, R);
//          }, 3*PI_MASS, 3.),
//        SetTitle("Smeared Ascoli-Jones \"Deck\" for f_{2}#pi S-wave.;M_{3#pi}"),
//        SetLineColor(kRed),
//        SetLineWidth(2),
//        SetFillStyle(0))->Draw("alp");
//   c1.SaveAs("/tmp/test_MAscoli.pdf");
//   SET4(
//        draw([&](double w)->double{
//            double mR = f2.GetM();
//            uint S = f2.GetL();
//            uint L = 2;
//            double m23 = 2*PI_MASS+(mR-2*PI_MASS)*(1.-exp(-1./mR*(w-3*PI_MASS)));
//            return bch.getProjectedDeck(mAsq, mBsq, w*w, mDsq, mtRsq, m23*m23, s, t,
//                                         J, 0, L, S, R);
//          }, 3*PI_MASS, 3.),
//        SetTitle("Smeared Ascoli-Jones \"Deck\" for f_{2}#pi D-wave.;M_{3#pi}"),
//        SetLineColor(kRed),
//        SetLineWidth(2),
//        SetFillStyle(0))->Draw("alp");
//   c1.SaveAs("/tmp/test_MAscoli.pdf");
// 
//   //   SET4(
// //        draw([&, ciso](double w)->double{
// //            std::function<double(double)> integrand = [&, ciso](double s1)->double {
// //              double pD = bch.getProjection(w*w, s1, J, L);
// //              return ciso.U(s1)*POW2(pD)*1./(2*M_PI) * 1./(8*M_PI)*sqrt(LAMBDA(w*w, s1, POW2(PI_MASS)))/w*w;
// //            };
// //            std::function<double(double)> drho = [&, ciso](double s1)->double {
// //              return ciso.U(s1)*1./(2*M_PI) * 1./(8*M_PI)*sqrt(LAMBDA(w*w, s1, POW2(PI_MASS)))/(w*w);
// //            };
// //            double intD = integrate(integrand, ciso.sth(), POW2(w-PI_MASS));
// //            double intRho = integrate(drho, ciso.sth(), POW2(w-PI_MASS));
// //            std::cout << "intRho done\n";
// //            return intD;  // sqrt(intD/intRho);
// //          }, 3*PI_MASS+0.01, 3., 20),
// //        SetTitle("Integrated Ascoli-Jones \"Deck\".;M_{3#pi}"),
// //        SetLineColor(kRed),
// //        SetLineWidth(2),
// //        SetFillStyle(0))->Draw("alp");
// //   c1.SaveAs("/tmp/test_MAscoli.pdf");
// 
//   TH2D h2("2dDeck", "Ascoli-Jones \"Deck\", #rho#pi F-wave;M_{3#pi};M_{2#pi}",
//           20, 3*PI_MASS, 3.,
//           20, 2*PI_MASS, 2.5);
//   for (int i=0; i < h2.GetXaxis()->GetNbins(); i++)
//     for (int j=0; j < h2.GetYaxis()->GetNbins(); j++) {
//       double w =  h2.GetXaxis()->GetBinCenter(i+1);
//       double m23 =  h2.GetYaxis()->GetBinCenter(j+1);
//       if (w < m23+PI_MASS) continue;
//       h2.SetBinContent(i+1, j+1, bch.getProjection(w*w, m23*m23, J, cL) );
//     }
//   h2.Draw("colz");
//   c1.SaveAs("/tmp/test_MAscoli.pdf)");
//   c1.SaveAs("/tmp/test_MAscoli.root");
//   
//       //  c1.Clear();
// //  c1.DivideSquare(Nch);
// //  for (uint i = 0; i < Nch; i++) {
// //    c1.cd(i+1);
// //    SET1(
// //         draw([&, i](double s)->double{ return bch[i]->getProjection(s); },
// //              bch[i]->sth(), POW2(2.5)),
// //         SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}",
// //                                  bch[i]->L(), bch[i]->S(), bch[i]->J(), bch[i]->lamP())))->Draw("al");
// //  }
// //  c1.SaveAs("/tmp/e.test_MDeck.pdf");
// //
// //  constexpr uint Npoints = 20;
// //  constexpr double from = POW2(3*PI_MASS)+0.1;
// //  constexpr double to = POW2(7.5);
// //
// //  bch[2]->makeLookupTable(*iso[2], PI_MASS, from, to, Npoints);
  
  return 0;
}
