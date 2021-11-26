// Copyright Misha Mikhasenko [12.2016], test
// test_ProjectedDeck.4:
//   In this example I export  many Deck projections

#include "MAscoli.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"

#include <utility>
#include <string>

#include "constants.h"
#include "mintegrate.h"
#include "mstructures.h"

#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "MIsobarChannel.h"

typedef struct {
  uint J;
  uint M;
  uint S;
  uint L;
  std::string title;
  std::string name;
} wave;
wave make_wave(uint _J, uint _M, uint _S, uint _L,
               const std::string &_title = "", const std::string &_name = "") {
  wave n; n.J = _J; n.M = _M; n.S = _S; n.L = _L; n.title = _title; n.name = _name;
  return n;
};

int main() {

  std::vector<wave> waves;
  waves.push_back(make_wave(1,0,1,0, "1^{++}0^{+}", "h2"));
  waves.push_back(make_wave(2,0,1,1, "2^{-+}0^{+}", "h33"));
  waves.push_back(make_wave(2,0,2,0, "2^{-+}0^{+}", "h26"));
  waves.push_back(make_wave(3,0,1,2, "3^{++}0^{+}", "h45"));
  waves.push_back(make_wave(3,0,2,1, "3^{++}0^{+}", "h49"));
  // waves.push_back(make_wave(3,0,2,1, "3^{++}0^{+}", "53"));
  waves.push_back(make_wave(4,0,1,3, "4^{-+}0^{+}", "h56"));

  const uint Niso = 3;
  double IIso[Niso] = {1./3, 1./sqrt(2), 1./3};
  double R = 5.0;
  double E_BEAM_LAB = 190;
  double mAsq = POW2(PI_MASS);
  double mBsq = POW2(PROT_MASS);
  double mDsq = POW2(PROT_MASS);
  double m1sq = POW2(PI_MASS);
  double mtRsq = POW2(PI_MASS);
  double tP = -0.1;

  // isobar
  MIsobar * iso[Niso] = {
    new MIsobarPiPiS,
    new MIsobar(RHO_MASS, RHO_WIDTH,  PI_MASS, PI_MASS, 0, 5),
    new MIsobar(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 1, 5)
  };
  for (uint i=0; i < Niso; i++) iso[i]->setIntU();
  
  // Deck table
  uint Ninterp = 100;
  std::vector<std::pair<double, double> > ltable(Ninterp);
  std::pair<double, double> ranges(0.6,2.5);
  TCanvas c1("c1","projections",0,0,1000,1000);
  c1.DivideSquare(waves.size());
  TFile fout("/tmp/anal.proj_0.100000-0.112853.root", "RECREATE");
  for (uint w = 0; w < waves.size(); w++) {
    double mIso = iso[waves[w].S]->GetM(); 
    for (uint t=0; t < Ninterp; t++) {
      double e = ranges.first + (ranges.second-ranges.first)/(Ninterp-1)*t;  /* triky */
      double s = e*e;
      double m23 = mIso;//2*PI_MASS+(mIso-2*PI_MASS)*(1.-exp(-1./mIso*(e-3*PI_MASS)));
      double val = MAscoli::getProjectedReducedDeck(waves[w].J, waves[w].M, waves[w].L,
                                                    POW2(m23), waves[w].S, R,
                                                    s, tP,
                                                    mtRsq,
                                                    2*PROT_MASS*E_BEAM_LAB,
                                                    mAsq, mBsq, mDsq,
                                                    m1sq);
      val *= IIso[waves[w].S] * iso[waves[w].S]->IntU();
      // val *= pow(s,-3./2);
      ltable[t] = std::make_pair(e, val);
      // std::cout << ltable[t].first << ", " << ltable[t].second << "\n";
    }
    TGraph *gr = new TGraph(Ninterp);
    gr->SetTitle(TString::Format("Deck projection, %s %s#pi %c-wave; M_{3#pi}, GeV",
                                 waves[w].title.c_str(),
                                 (std::vector<std::string>{"(#pi#pi)_S", "#rho", "f_{2}"})[waves[w].S].c_str(),
                                 (std::vector<char>{'S', 'P', 'D', 'F', 'G', 'H'})[waves[w].L]));
    gr->SetName(waves[w].name.c_str());
    for (uint t=0; t < Ninterp; t++) {
      gr->GetX()[t] = ltable[t].first;
      double pbu = ((ltable[t].first < PI_MASS+mIso) ? 0 :
                    sqrt(LAMBDA(POW2(ltable[t].first),
                                POW2(PI_MASS),
                                POW2(mIso)))/(2*ltable[t].first));

      gr->GetY()[t] = POW2(ltable[t].second) *
        /* Phase space */
         1./(8*M_PI)* 2*pbu/ltable[t].first  *
        /* Blatt-Weisskopf is not needed */
        /* Clebsch */
        IIso[waves[w].S]*sqrt(2);
    }
    gr->Write();
    c1.cd(w+1); gr->Draw("al");
  }
  c1.SaveAs("/tmp/waves.pdf");
  fout.Close();
  
  return 1;
}
