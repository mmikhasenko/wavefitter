// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TAxis.h"
#include "Math/SpecFuncMathMore.h"

#include "MIsobar.h"
#include "MDeck.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.h"
#include "mintegrate.h"
#include "deflib.h"

#include "dFunction.hpp"

std::vector<std::pair<double, double> >
makeLookupTable(const MDeck &b, const MIsobar &iso, double from, double to, uint Npoints);


int main(int argc, char *argv[]) {
  // define all variables

  uint Sp = 1;
  uint lamP = 0;
  uint J = 2;

  double m1sq = POW2(PI_MASS);
  double m2sq = -0.01;
  double m4sq = POW2(PI_MASS);
  double mtsq = POW2(PI_MASS);

  double R = 1.;
  // Isobars
  MIsobar rho(0.77, 0.15, 0.14, 0.14, 1, 5.);
  MIsobar  f2(1.23, 0.2,  0.14, 0.14, 2, 5.);
  constexpr uint Nch = 4;
  const MIsobar *iso[Nch] = {&rho, &rho, &f2, &f2};

  // Deck channels
  MDeck *bch[Nch];
  bch[0] = new MDeck(m1sq, m2sq, POW2(iso[0]->GetM()), m4sq, mtsq, J, 1, Sp, lamP, 1, 0, R);
  bch[1] = new MDeck(m1sq, m2sq, POW2(iso[1]->GetM()), m4sq, mtsq, J, 3, Sp, lamP, 1, 0, R);
  bch[2] = new MDeck(m1sq, m2sq, POW2(iso[2]->GetM()), m4sq, mtsq, J, 0, Sp, lamP, 2, 0, R);
  bch[3] = new MDeck(m1sq, m2sq, POW2(iso[3]->GetM()), m4sq, mtsq, J, 2, Sp, lamP, 2, 0, R);
  constexpr int colors[Nch] = {kBlue, kOrange, kGreen, kMagenta};

  TCanvas c1("c1");
  /*  Angular distribution */
  TMultiGraph m;
  for (uint i = 0; i < Nch; i += 2) {
    m.Add(
          SET4(
               draw([&, i](double z)->double{
                   return real(bch[i]->getValue(2.1, z));
                 }, -1, 1),
               SetTitle(TString::Format("Deck %d", i)),
               SetLineColor(colors[i]),
               SetLineWidth(2),
               SetFillStyle(0)) );
  }
  m.SetTitle("cos#theta dependence for s = 2.1 GeV^{2}");
  m.Draw("al");
  c1.BuildLegend(0.11, 0.11, 0.4, 0.25);
  c1.SaveAs("/tmp/a.test_MDeck.pdf");

  /* Mass distribution */
  c1.Clear();
  c1.DivideSquare(Nch);
  for (uint i = 0; i < Nch; i++) {
    c1.cd(i+1);
    SET1(
         draw([&, i](double s)->double{ return bch[i]->getProjection(s); },
              bch[i]->sth(), POW2(2.5)),
         SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}",
                                  bch[i]->L(), bch[i]->S(), bch[i]->J(), bch[i]->lamP())))->Draw("al");
  }
  c1.SaveAs("/tmp/e.test_MDeck.pdf");

  constexpr uint Npoints = 20;
  constexpr double from = POW2(3*PI_MASS)+0.1;
  constexpr double to = POW2(7.5);

  bch[2]->makeLookupTable(*iso[2], PI_MASS, from, to, Npoints);

//  c1.Clear();
//  draw([&](double s)->double{
//      return bch[2]->getPrecalculated(s);
//    }, from, to)->Draw("al");
//  c1.SaveAs("/tmp/g.test_MDeck.pdf");

  c1.Clear();
  c1.DivideSquare(Nch);
  for (uint i=0; i < Nch; i++) {
    c1.cd(i+1);
    SET2(draw([&, i](double s)->double{
          return fabs(bch[i]->getProjection(s))*sqrt(LAMBDA(s, POW2(iso[i]->GetM()), m4sq)/(4*s));
        }, bch[i]->sth(), to),
      SetTitle(TString::Format("J^{PC} = %u^{-+} S=%u %c-wave damping with left pole M = %2.2f GeV",
                               bch[i]->J(), iso[i]->GetL(),
                               (std::vector<char>{'S', 'P', 'D', 'F'})[bch[i]->L()],
                               1./R)),
      GetXaxis()->SetTitle("s (GeV)"))->Draw("al");
  }
  c1.SaveAs("/tmp/g.test_MDeck.pdf");

  return 0;
}

// std::vector<std::pair<double, double> > ltable = makeLookupTable(*bch[2], *iso[2],
//                                                                  from, to, Npoints);
std::vector<std::pair<double, double> >
makeLookupTable(const MDeck &b, const MIsobar &iso, double from, double to, uint Npoints) {
  // first create integrand
  double s;
  std::function<double(double)> integrand = [&](double s3)->double {
    double pD = b.getProjection(s, s3);
    return iso.U(s3)*POW2(pD)*1./(2*M_PI);
  };
  // then create and fill table
  std::vector<std::pair<double, double> > ltable;
  ltable.resize(Npoints);
  for (uint i = 0; i < Npoints; i++) {
    s = from + (to-from)/(Npoints-1)*i;
    ltable[i] = std::make_pair(s,
                               integrate(integrand, iso.sth(), POW2(sqrt(s)-PI_MASS)));
    std::cout << "--> " << s << ", " << ltable[i].second << std::endl;
  }
  return ltable;
}
