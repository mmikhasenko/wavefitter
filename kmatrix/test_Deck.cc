// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/SpecFuncMathMore.h"

#include "MIsobar.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.hh"
#include "mintegrate.hh"
#include "deflib.h"

#include "dFunction.hpp"


uint i_m1sq;
uint i_m2sq;
uint i_m3sq;
uint i_m4sq;
uint i_mtsq;

cd getDeck(double s, double z,
           uint Sp, int lamP,
           uint S , int lamS);

double getProjectedDeck(double s,
                 uint J, uint L,
                 uint Sp, int lamP,
                 uint S);

int main(int argc, char *argv[]) {
  // define all variables

  i_m1sq = MParKeeper::gI()->add("m1sq", POW2(PI_MASS));
  i_m2sq = MParKeeper::gI()->add("m2sq", -0.01);
  i_m3sq = MParKeeper::gI()->add("m3sq", POW2(RHO_MASS));
  i_m4sq = MParKeeper::gI()->add("m4sq", POW2(PI_MASS));
  i_mtsq = MParKeeper::gI()->add("mtsq", POW2(PI_MASS));

  double m3sq = MParKeeper::gI()->get(i_m3sq);
  double m4sq = MParKeeper::gI()->get(i_m4sq);

  uint _Sp = 1;
  uint _lamP = 0;
  uint _S = 1;
  uint _lamS = 0;
  uint _J = 2;
  uint _L = 3;

  TCanvas c1("c1");
  combine(SET1(
               draw([&](double z)->double{
                   return real(getDeck(2.1, z, _Sp, _lamP, _S, _lamS));
                 }, -1, 1), SetLineColor(kBlack)),
          SET1(
               draw([&](double z)->double{
                   return imag(getDeck(2.1, z, _Sp, _lamP, _S, _lamS));
                 }, -1, 1), SetLineColor(kRed) ) )->Draw("al");
  c1.SaveAs("/tmp/a.test_Deck.pdf");


  /* part with integrals */
  std::cout << "int = " << cintegrate([&](double z)->cd {
      return getDeck(2.1, z, _Sp, _lamP, _S, _lamS);
    }, -1, 1) << "\n";

  std::function<cd(double)> integB = [&](double s)->cd {
    return cintegrate([&](double z)->cd {
        return getDeck(s, z, _Sp, _lamP, _S, _lamS);
      }, -1, 1);
  };

  /* part with integrals */
  std::cout << "cross check int = " << integB(2.1) << "\n";

  // plot integral value
  combine(
          SET2(
               draw([&](double s)->double{
                   return real(cintegrate([&](double z)->cd {
                         return getDeck(s, z, _Sp, _lamP, _S, -1);
                       }, -1, 1));},
                 POW2(sqrt(m3sq)+sqrt(m4sq)), 2.2),
               SetTitle("Integrated Deck b_{S=1,-1}^{M=0}"),
               SetLineColor(kGreen)),
          SET2(
               draw([&](double s)->double{
                   return real(cintegrate([&](double z)->cd {
                         return getDeck(s, z, _Sp, _lamP, _S, 0);
                       }, -1, 1));},
                 POW2(sqrt(m3sq)+sqrt(m4sq)), 2.2),
               SetTitle("Integrated Deck b_{S=1,-1}^{M=0}"),
               SetLineColor(kOrange)),
          SET2(
               draw([&](double s)->double{
                   return real(cintegrate([&](double z)->cd {
                         return getDeck(s, z, _Sp, _lamP, _S, +1);
                       }, -1, 1));},
                 POW2(sqrt(m3sq)+sqrt(m4sq)), 2.2),
               SetTitle("Integrated Deck b_{S=1,+1}^{M=0}"),
               SetLineColor(kBlue)) )->Draw("al");
  c1.SaveAs("/tmp/b.test_Deck.pdf");

  // plot projection value
  SET1(
       draw([&](double s)->double{
           return getProjectedDeck(s, _J, _L, _Sp, _lamP, _S);
         },
         POW2(sqrt(m3sq)+sqrt(m4sq))+1e-4, 2.2),
       SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP)) )->Draw("al");
  c1.SaveAs("/tmp/c.test_Deck.pdf");

  // plot projection value
  SET1(
       draw([&](double s)->double{
           return getProjectedDeck(s, _J, _L, _Sp, _lamP, _S) /
             pow(LAMBDA(s, m3sq, m4sq)/(4*s), _L/2.);
         },
         POW2(sqrt(m3sq)+sqrt(m4sq))+1e-5, 2.2),
       SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP)) )->Draw("al");
  c1.SaveAs("/tmp/d.test_Deck.pdf");

  // complete plot fo projections
  /************************************************************/
  {
    const double Hlim = POW2(2.5);
    c1.Clear();
    c1.DivideSquare(4);
    // rho pi P-wave
    int _S = 1; int _L = 1;
    c1.cd(1);
    SET1(
         draw([&](double s)->double{
             return getProjectedDeck(s, _J, _L, _Sp, _lamP, _S);
           },
           POW2(sqrt(m3sq)+sqrt(m4sq))+1e-4, Hlim),
         SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP)) )->Draw("al");
    // rho pi F-wave
    _S = 1; _L = 3;
    c1.cd(2);
    SET1(
         draw([&](double s)->double{
             return getProjectedDeck(s, _J, _L, _Sp, _lamP, _S);
           },
           POW2(sqrt(m3sq)+sqrt(m4sq))+1e-4, Hlim),
         SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP)) )->Draw("al");

    // f2 pi S-wave
    MParKeeper::gI()->set(i_m3sq, POW2(F2_MASS));
    m3sq = MParKeeper::gI()->get(i_m3sq);
    _S = 2; _L = 0;
    c1.cd(3);
    SET1(
         draw([&](double s)->double{
             return getProjectedDeck(s, _J, _L, _Sp, _lamP, _S);
           },
           POW2(sqrt(m3sq)+sqrt(m4sq))+1e-4, Hlim),
         SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP)) )->Draw("al");
    // f2 pi D-wave
    _S = 2; _L = 2;
    c1.cd(4);
    SET1(
         draw([&](double s)->double{
             return getProjectedDeck(s, _J, _L, _Sp, _lamP, _S);
           },
           POW2(sqrt(m3sq)+sqrt(m4sq))+1e-4, Hlim),
         SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP)) )->Draw("al");
//    // f2 pi H-wave
//    _S = 2; _L = 4;
//    c1.cd(5);
//    SET1(
//         draw([&](double s)->double{
//             return getProjectedDeck(s, _J, _L, _Sp, _lamP, _S);
//           },
//           POW2(sqrt(m3sq)+sqrt(m4sq))+1e-4, Hlim),
//         SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP)) )->Draw("al");

    c1.SaveAs("/tmp/e.test_Deck.pdf");
  }


  // stable
  MIsobar rho(RHO_MASS, 0.15, PI_MASS, PI_MASS, 1, 5.);
  MIsobar  f2(F2_MASS , 0.2,  PI_MASS, PI_MASS, 2, 5.);
  // complete plot fo projections
  /************************************************************/
  {
    double s = 2.2;
//    s = 0.22;
//    std::cout << integrate(integrand, rho.sth(), POW2(sqrt(s)-PI_MASS)) << "\n";

    c1.Clear();
    c1.DivideSquare(4);
    // rho pi P-wave
    int _S = 1; int _L = 1;
    std::function<double(double)> integrandRHO = [&](double s12)->double {
      MParKeeper::gI()->set(i_m3sq, s12);
      double pD = getProjectedDeck(s, _J, _L, _Sp, _lamP, _S);
      return rho.U(s12)*POW2(pD)*1./(2*M_PI);
    };
    c1.cd(1);
    TGraph *g1 =
      SET2(
           draw([&](double sr)->double{
               s = sr;
               return integrate(integrandRHO, rho.sth(), POW2(sqrt(sr)-PI_MASS));
             },
             POW2(3*PI_MASS)+0.1, POW2(2.5), 40),
           SetLineColor(kBlack),
           SetTitle(TString::Format("Intensity of Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP) ) );
    gr1->Draw("al");
    std::cout << "rho pi P-wave" << std::endl;
    for (int i = 0; i < gr1->GetN(); i++) std::cout << gr1->GetX()[i] << gr1->GetY()[i] << std::endl;
    // f2 pi S-wave
    _S = 2; _L = 0;
    std::function<double(double)> integrandF2 = [&](double s12)->double {
      MParKeeper::gI()->set(i_m3sq, s12);
      double pD = getProjectedDeck(s, _J, _L, _Sp, _lamP, _S);
      return f2.U(s12)*POW2(pD)*1./(2*M_PI);
    };
    c1.cd(3);
    TGraph *g3 =
      SET2(
           draw([&](double sr)->double{
               s = sr;
               return integrate(integrandF2, rho.sth(), POW2(sqrt(sr)-PI_MASS));
             },
             POW2(3*PI_MASS)+0.1, POW2(2.5), 40),
           SetLineColor(kBlack),
           SetTitle(TString::Format("Intensity of Deck projections b_{%d%d}^{%d%d}", _L, _S, _J, _lamP) ) );
    g3->Draw("al");
    std::cout << "f2 pi S-wave" << std::endl;
    for (int i = 0; i < gr3->GetN(); i++) std::cout << gr3->GetX()[i] << gr3->GetY()[i] << std::endl;
    c1.SaveAs("/tmp/g.test_Deck.pdf");
  }

  return 0;
}


cd getDeck(double s, double z,
           uint Sp, int lamP,
           uint S , int lamS) {
  double m1sq = MParKeeper::gI()->get(i_m1sq);
  double m2sq = MParKeeper::gI()->get(i_m2sq);
  double m3sq = MParKeeper::gI()->get(i_m3sq);
  double m4sq = MParKeeper::gI()->get(i_m4sq);
  double mtsq = MParKeeper::gI()->get(i_mtsq);

  // calculate dependent variables
  double t = m1sq + m3sq - ((s + m3sq - m4sq)*(s + m1sq - m2sq))/(2.*s) +
    sqrt(LAMBDA(s, m1sq, m2sq)*LAMBDA(s, m3sq, m4sq))/(2.*s)*z;
  cd CosChi2 = ((s + m2sq - m1sq)*(t + m2sq - m4sq) - 2.*m2sq*(m2sq + m3sq - m1sq - m4sq)) /
    sqrt(LAMBDA(s, m1sq, m2sq)*LAMBDA(t, m2sq, m4sq));
  cd CosChi3 = ((s + m3sq - m4sq)*(t + m3sq - m1sq) - 2.*m3sq*(m2sq + m3sq - m1sq - m4sq)) /
    sqrt(LAMBDA(s, m3sq, m4sq)*LAMBDA(t, m1sq, m3sq));
  cd Chi2 = acos(CosChi2);
  cd Chi3 = acos(CosChi3);
  double ClebschP = pow(2., Sp/2.)*TMath::Factorial(Sp)/sqrt(TMath::Factorial(2*Sp));
  double ClebschS = pow(2., S/2.)*TMath::Factorial(S)/sqrt(TMath::Factorial(2*S));

  // damping
  const double R = 5.;

  // (* main fucnction *)
  cd f = ClebschP *
    pow(cd(LAMBDA(t, m1sq, m3sq))/(LAMBDA(t, m1sq, m3sq) - 4.*t/POW2(R)), S/2.) *
    1./(mtsq - t) *
    rpwa::dFunction<cd>(2*Sp, 0, 2*lamP, Chi2) *
    rpwa::dFunction<cd>(2*S,  0, 2*lamS, Chi3) *
    pow(cd(LAMBDA(t, m2sq, m4sq))/(LAMBDA(t, m2sq, m4sq) - 4.*t/POW2(R)), Sp/2.) /
    pow(cd(m2sq), lamP/2.) * /* to suppress complexity*/
    ClebschS;
  return f;
}


double getProjectedDeck(double s,
                        uint J, uint L,
                        uint Sp, int lamP,
                        uint S) {
  // integrate
  int minSJ = (S < J) ? S : J;
  double int_val =
    integrate([&, s, J, L, Sp, lamP, S, minSJ](double z)->double{
      cd val = 0;
      for (int lamS = -minSJ; lamS <= minSJ; lamS++)
        val +=
          ROOT::Math::wigner_3j(2*L, 2*S, 2*J, 2*0, 2*lamS, -2*lamS) *
          getDeck(s, z, Sp, lamP, S, -lamS) *
          rpwa::dFunction<double> (2*J, 2*lamP, 2*lamS, acos(z));
      return real(val);
      }, -1., 1.);
  int_val *= sqrt(2*J+1) * /*because of Clebsch and 3j*/
    sqrt(M_PI*(2*L+1)); /*because of normalisation*/

  return int_val;
}
