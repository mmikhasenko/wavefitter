// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TAxis.h"
#include "Math/SpecFuncMathMore.h"

#include "MIsobar.h"
#include "MHelicityDeck.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.h"
#include "mintegrate.h"
#include "deflib.h"

std::vector<std::pair<double, double> >
makeLookupTable(const MHelicityDeck &b, const MIsobar &iso, double from, double to, uint Npoints);


int main(int argc, char *argv[]) {
  // define all variables
}
// [this Deck I do not need]   uint Sp = 1;
// [this Deck I do not need]   uint lamP = 0;
// [this Deck I do not need]   uint J = 2;
// [this Deck I do not need] 
// [this Deck I do not need]   double m1sq = POW2(PI_MASS);
// [this Deck I do not need]   double m2sq = -0.01;
// [this Deck I do not need]   double m4sq = POW2(PI_MASS);
// [this Deck I do not need]   double mtsq = POW2(PI_MASS);
// [this Deck I do not need] 
// [this Deck I do not need]   double R = 1.;
// [this Deck I do not need]   // Isobars
// [this Deck I do not need]   MIsobar rho(0.77, 0.15, 0.14, 0.14, 1, 5.);
// [this Deck I do not need]   MIsobar  f2(1.23, 0.2,  0.14, 0.14, 2, 5.);
// [this Deck I do not need]   constexpr uint Nch = 4;
// [this Deck I do not need]   const MIsobar *iso[Nch] = {&rho, &rho, &f2, &f2};
// [this Deck I do not need] 
// [this Deck I do not need]   // Deck channels
// [this Deck I do not need]   MDeck *bch[Nch];
// [this Deck I do not need]   bch[0] = new MDeck(m1sq, m2sq, POW2(iso[0]->GetM()), m4sq, mtsq, J, 1, Sp, lamP, 1, 0, R);
// [this Deck I do not need]   bch[1] = new MDeck(m1sq, m2sq, POW2(iso[1]->GetM()), m4sq, mtsq, J, 3, Sp, lamP, 1, 0, R);
// [this Deck I do not need]   bch[2] = new MDeck(m1sq, m2sq, POW2(iso[2]->GetM()), m4sq, mtsq, J, 0, Sp, lamP, 2, 0, R);
// [this Deck I do not need]   bch[3] = new MDeck(m1sq, m2sq, POW2(iso[3]->GetM()), m4sq, mtsq, J, 2, Sp, lamP, 2, 0, R);
// [this Deck I do not need]   constexpr int colors[Nch] = {kBlue, kOrange, kGreen, kMagenta};
// [this Deck I do not need] 
// [this Deck I do not need]   TCanvas c1("c1");
// [this Deck I do not need]   /*  Angular distribution */
// [this Deck I do not need]   TMultiGraph m;
// [this Deck I do not need]   for (uint i = 0; i < Nch; i += 2) {
// [this Deck I do not need]     m.Add(
// [this Deck I do not need]           SET4(
// [this Deck I do not need]                draw([&, i](double z)->double{
// [this Deck I do not need]                    return real(bch[i]->getValue(2.1, z));
// [this Deck I do not need]                  }, -1, 1),
// [this Deck I do not need]                SetTitle(TString::Format("Deck %d", i)),
// [this Deck I do not need]                SetLineColor(colors[i]),
// [this Deck I do not need]                SetLineWidth(2),
// [this Deck I do not need]                SetFillStyle(0)) );
// [this Deck I do not need]   }
// [this Deck I do not need]   m.SetTitle("cos#theta dependence for s = 2.1 GeV^{2}");
// [this Deck I do not need]   m.Draw("al");
// [this Deck I do not need]   c1.BuildLegend(0.11, 0.11, 0.4, 0.25);
// [this Deck I do not need]   c1.SaveAs("/tmp/a.test_MDeck.pdf");
// [this Deck I do not need] 
// [this Deck I do not need]   /* Mass distribution */
// [this Deck I do not need]   c1.Clear();
// [this Deck I do not need]   c1.DivideSquare(Nch);
// [this Deck I do not need]   for (uint i = 0; i < Nch; i++) {
// [this Deck I do not need]     c1.cd(i+1);
// [this Deck I do not need]     SET1(
// [this Deck I do not need]          draw([&, i](double s)->double{ return bch[i]->getProjection(s); },
// [this Deck I do not need]               bch[i]->sth(), POW2(2.5)),
// [this Deck I do not need]          SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}",
// [this Deck I do not need]                                   bch[i]->L(), bch[i]->S(), bch[i]->J(), bch[i]->lamP())))->Draw("al");
// [this Deck I do not need]   }
// [this Deck I do not need]   c1.SaveAs("/tmp/e.test_MDeck.pdf");
// [this Deck I do not need] 
// [this Deck I do not need]   constexpr uint Npoints = 20;
// [this Deck I do not need]   constexpr double from = POW2(3*PI_MASS)+0.1;
// [this Deck I do not need]   constexpr double to = POW2(7.5);
// [this Deck I do not need] 
// [this Deck I do not need]   bch[2]->makeLookupTable(*iso[2], PI_MASS, from, to, Npoints);
// [this Deck I do not need] 
// [this Deck I do not need] //  c1.Clear();
// [this Deck I do not need] //  draw([&](double s)->double{
// [this Deck I do not need] //      return bch[2]->getPrecalculated(s);
// [this Deck I do not need] //    }, from, to)->Draw("al");
// [this Deck I do not need] //  c1.SaveAs("/tmp/g.test_MDeck.pdf");
// [this Deck I do not need] 
// [this Deck I do not need]   c1.Clear();
// [this Deck I do not need]   c1.DivideSquare(Nch);
// [this Deck I do not need]   for (uint i=0; i < Nch; i++) {
// [this Deck I do not need]     c1.cd(i+1);
// [this Deck I do not need]     SET2(draw([&, i](double s)->double{
// [this Deck I do not need]           return fabs(bch[i]->getProjection(s))*sqrt(LAMBDA(s, POW2(iso[i]->GetM()), m4sq)/(4*s));
// [this Deck I do not need]         }, bch[i]->sth(), to),
// [this Deck I do not need]       SetTitle(TString::Format("J^{PC} = %u^{-+} S=%u %c-wave damping with left pole M = %2.2f GeV",
// [this Deck I do not need]                                bch[i]->J(), iso[i]->GetL(),
// [this Deck I do not need]                                (std::vector<char>{'S', 'P', 'D', 'F'})[bch[i]->L()],
// [this Deck I do not need]                                1./R)),
// [this Deck I do not need]       GetXaxis()->SetTitle("s (GeV)"))->Draw("al");
// [this Deck I do not need]   }
// [this Deck I do not need]   c1.SaveAs("/tmp/g.test_MDeck.pdf");
// [this Deck I do not need] 
// [this Deck I do not need]   return 0;
// [this Deck I do not need] }
// [this Deck I do not need] 
// [this Deck I do not need] // std::vector<std::pair<double, double> > ltable = makeLookupTable(*bch[2], *iso[2],
// [this Deck I do not need] //                                                                  from, to, Npoints);
// [this Deck I do not need] std::vector<std::pair<double, double> >
// [this Deck I do not need] makeLookupTable(const MDeck &b, const MIsobar &iso, double from, double to, uint Npoints) {
// [this Deck I do not need]   // first create integrand
// [this Deck I do not need]   double s;
// [this Deck I do not need]   std::function<double(double)> integrand = [&](double s3)->double {
// [this Deck I do not need]     double pD = b.getProjection(s, s3);
// [this Deck I do not need]     return iso.U(s3)*POW2(pD)*1./(2*M_PI);
// [this Deck I do not need]   };
// [this Deck I do not need]   // then create and fill table
// [this Deck I do not need]   std::vector<std::pair<double, double> > ltable;
// [this Deck I do not need]   ltable.resize(Npoints);
// [this Deck I do not need]   for (uint i = 0; i < Npoints; i++) {
// [this Deck I do not need]     s = from + (to-from)/(Npoints-1)*i;
// [this Deck I do not need]     ltable[i] = std::make_pair(s,
// [this Deck I do not need]                                integrate(integrand, iso.sth(), POW2(sqrt(s)-PI_MASS)));
// [this Deck I do not need]     std::cout << "--> " << s << ", " << ltable[i].second << std::endl;
// [this Deck I do not need]   }
// [this Deck I do not need]   return ltable;
// [this Deck I do not need] }
