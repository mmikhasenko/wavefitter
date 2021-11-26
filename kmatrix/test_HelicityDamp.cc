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

  uint Sp = 1;
  uint lamP = 0;
  uint J = 2;
  uint L = 3;
  uint S = 1;

  double m1sq = POW2(PI_MASS);
  double m2sq = -0.01;
  double m4sq = POW2(PI_MASS);
  double mtsq = POW2(PI_MASS);

  // Isobars
  MIsobar rho(0.77, 0.15, 0.14, 0.14, S, 5.);
  MIsobar *iso = &rho;

  // Deck channels
  constexpr uint Nch = 4;
  double R[Nch] = {0.0001, 1., 5., 50.};
  MHelicityDeck *bch[Nch];
  bch[0] = new MHelicityDeck(m1sq, m2sq, POW2(iso->GetM()), m4sq, mtsq, J, L, Sp, lamP, iso->GetL(), 0, R[0]);
  bch[1] = new MHelicityDeck(m1sq, m2sq, POW2(iso->GetM()), m4sq, mtsq, J, L, Sp, lamP, iso->GetL(), 0, R[1]);
  bch[2] = new MHelicityDeck(m1sq, m2sq, POW2(iso->GetM()), m4sq, mtsq, J, L, Sp, lamP, iso->GetL(), 0, R[2]);
  bch[3] = new MHelicityDeck(m1sq, m2sq, POW2(iso->GetM()), m4sq, mtsq, J, L, Sp, lamP, iso->GetL(), 0, R[3]);
  // colors
  constexpr int colors[Nch] = {kBlue, kOrange, kGreen, kMagenta};

  TCanvas c1("c1");
  /*  Angular distribution */
  TMultiGraph ma;
  for (uint i = 0; i < Nch; i++) {
    ma.Add(
          SET4(
               draw([&, i](double z)->double{
                   return real(bch[i]->getValue(2.1, z));
                 }, -1, 1),
               SetTitle(TString::Format("Deck with R = %2.1f", R[i])),
               SetLineColor(colors[i]),
               SetLineWidth(2),
               SetFillStyle(0)) );
  }
  ma.SetTitle("cos(#theta) dependence for Deck for s = 2.1 GeV^{2}");
  ma.Draw("al");
  ma.GetXaxis()->SetTitle("cos(#theta)");
  c1.BuildLegend(0.11, 0.11, 0.4, 0.25);
  c1.SaveAs("/tmp/a.test_Damp.pdf");

  /* Mass distribution */
  c1.Clear();
  TMultiGraph mb;
  const double sNorm = POW2(2.2);
  double norm[Nch]; for (uint i = 0; i < Nch; i++) { norm[i] = bch[i]->getProjection(sNorm); std::cout << "norm[i] = " << norm[i] << "\n"; }
  for (uint i = 0; i < Nch; i++) {
    mb.Add(
          SET4(
               draw([&, i](double s)->double{ return bch[i]->getProjection(s)/norm[i]; },
                    bch[i]->sth(), POW2(2.5), 100),
               SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d} with R = %2.1f",
                                        bch[i]->L(), bch[i]->S(), bch[i]->J(), bch[i]->lamP(), R[i])),
               SetLineColor(colors[i]),
               SetLineWidth(2),
               SetFillStyle(0)));
  }
  mb.SetTitle("Deck projection(s) for different damp-parameter");
  mb.Draw("al");
  mb.GetXaxis()->SetTitle("s (GeV^{2})");
  c1.BuildLegend(0.61, 0.11, 0.89, 0.25);
  c1.SaveAs("/tmp/b.test_Damp.pdf");

  /* Change with s3 */
  const uint Ns3Points = 50;
  double s3values[Ns3Points];
  for (uint j = 0; j < Ns3Points; j++) {
    s3values[j] = 0.1 + 1.4/(Ns3Points-1) * j;
    double s3 = s3values[j];

    c1.Clear();
    TMultiGraph mc;

    const double sNorm = POW2(2.2);
    double norm[Nch];
    for (uint i = 0; i < Nch; i++) { norm[i] = fabs(bch[i]->getProjection(sNorm, s3)); std::cout << "norm[i] = " << norm[i] << "\n"; }
    for (uint i = 0; i < Nch; i++) {
      mc.Add(
             SET4(
                  draw([&, i](double s)->double{ return bch[i]->getProjection(s, s3) / norm[i]; },
                       POW2(sqrt(s3)+PI_MASS), POW2(2.5), 100),
                  SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d} with R = %2.1f, M3 = %2.2f",
                                           bch[i]->L(), bch[i]->S(), bch[i]->J(), bch[i]->lamP(), R[i], sqrt(s3))),
                  SetLineColor(colors[i]),
                  SetLineWidth(2),
                  SetFillStyle(0)));
    }
    mc.SetTitle("Deck projection for different damp-parameter");
    mc.Draw("al");
    c1.BuildLegend(0.51, 0.11, 0.89, 0.25);
    c1.SaveAs(TString::Format("/tmp/c%02d.test_Damp.pdf", j));
  }

  /* Change with s3, remove threshold */
  for (uint j = 0; j < Ns3Points; j++) {
    double s3 = s3values[j];

    c1.Clear();
    TMultiGraph mc;

    const double sNorm = POW2(2.2);
    double norm[Nch];
    for (uint i = 0; i < Nch; i++) {
      norm[i] = fabs(bch[i]->getProjection(sNorm, s3) / pow(LAMBDA(sNorm, s3, m4sq)/(4.*sNorm), S/2.));
                     std::cout << "norm[i] = " << norm[i] << "\n";
                     }
    for (uint i = 0; i < Nch; i++) {
      mc.Add(
             SET4(
                  draw([&, i](double s)->double{ return bch[i]->getProjection(s, s3) / pow(LAMBDA(s, s3, m4sq)/(4.*s), L/2.) / norm[i]; },
                       POW2(sqrt(s3)+PI_MASS)+1e-2, POW2(2.5), 100),
                  SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d} with R = %2.1f, M_{3} = %2.2f",
                                           bch[i]->L(), bch[i]->S(), bch[i]->J(), bch[i]->lamP(), R[i], sqrt(s3))),
                  SetLineColor(colors[i]),
                  SetLineWidth(2),
                  SetFillStyle(0)));
    }
    mc.SetTitle("Deck projection for different damp-parameter");
    mc.Draw("al");
    c1.BuildLegend(0.51, 0.11, 0.89, 0.25);
    c1.SaveAs(TString::Format("/tmp/d%02d.test_Damp.pdf", j));
  }

  /* s3 dependence */
  TMultiGraph me;
  const double s0 = POW2(1.2);
  for (uint i = 0; i < Nch; i++) {
    me.Add(
           SET4(
                draw([&, i](double s3)->double{ return bch[i]->getProjection(s0, s3) / pow(LAMBDA(s0, s3, m4sq)/(4.*s0), L/2.); },
                     POW2(2*PI_MASS)+1e-2, POW2(sqrt(s0)-PI_MASS)-1e-2, 100),
                SetTitle(TString::Format("Deck projections b_{%d%d}^{%d%d}(%2.2f, m_{3}^{2}) with R = %2.1f",
                                         bch[i]->L(), bch[i]->S(), bch[i]->J(), bch[i]->lamP(), s0, R[i])),
                SetLineColor(colors[i]),
                SetLineWidth(2),
                SetFillStyle(0)));
  }
  me.SetTitle("Deck projection for different damp-parameter");
  me.Draw("al");
  me.GetXaxis()->SetTitle("m_{3}^{2} (GeV^{2})");
  c1.BuildLegend(0.51, 0.11, 0.89, 0.25);
  c1.SaveAs("/tmp/e.test_Damp.pdf");

  return 0;
}
