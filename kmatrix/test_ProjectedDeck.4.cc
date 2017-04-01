// Copyright Misha Mikhasenko [12.2016], test
// test_ProjectedDeck.4:
//   In this example I export  many Deck projections

#include "MAscoli.h"
#include "TCanvas.h"
#include "TGraph.h"

#include <utility>
#include <string>

#include "constants.h"
#include "mintegrate.h"
#include "mstructures.h"

#include "MIsobar.h"
#include "MIsobarChannel.h"

typedef struct {
  uint J;
  uint M;
  uint S;
  uint L;
  std::string title;
} wave;
wave make_wave(uint _J, uint _M, uint _S, uint _L, const std::string &_title = "") {
  wave n; n.J = _J; n.M = _M; n.S = _S; n.L = _L; n.title = _title;
  return n;
};

int main() {

  std::vector<wave> waves;
  /* 1++ */
  waves.push_back(make_wave(1,0,1,0, "1^{++}0^{+}"));
  /* 0-+ */
  waves.push_back(make_wave(0,0,0,0, "0^{-+}0^{+}"));
  waves.push_back(make_wave(0,0,0,0, "0^{-+}0^{+}"));
  waves.push_back(make_wave(0,0,1,1, "0^{-+}0^{+}"));
  waves.push_back(make_wave(0,0,0,0, "0^{-+}0^{+}"));
  waves.push_back(make_wave(0,0,2,2, "0^{-+}0^{+}"));
  /* 2-+ */
  waves.push_back(make_wave(2,0,2,0, "2^{-+}0^{+}"));
  waves.push_back(make_wave(2,0,2,2, "2^{-+}0^{+}"));
  waves.push_back(make_wave(2,0,1,1, "2^{-+}0^{+}"));
  waves.push_back(make_wave(2,0,1,3, "2^{-+}0^{+}"));
  /* 1-- */
  waves.push_back(make_wave(1,0,1,0, "1^{--}0^{+}"));
  waves.push_back(make_wave(1,0,1,2, "1^{--}0^{+}"));
  waves.push_back(make_wave(1,0,2,1, "1^{--}0^{+}"));
  waves.push_back(make_wave(1,0,2,3, "1^{--}0^{+}"));
  
  double mIso[3] = {0.8,RHO_MASS,F2_MASS};
  double R = 5.0;
  double E_BEAM_LAB = 190;
  double mAsq = POW2(PI_MASS);
  double mBsq = POW2(PROT_MASS);
  double mDsq = POW2(PROT_MASS);
  double m1sq = POW2(PI_MASS);
  double mtRsq = POW2(PI_MASS);
  double tP = -0.1;

  // isobar
  // MIsobar f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, RI);
  // MIsobar rho_iso(RHO_MASS, RHO_WIDTH,  PI_MASS, PI_MASS, S, 5);
  // MIsobarChannel *_iso = new MIsobarChannel(rho_iso, PI_MASS, L);
  // double sth = _iso->sth();

  // Deck table
  uint Ninterp = 100;
  std::vector<std::pair<double, double> > ltable(Ninterp);
  std::pair<double, double> ranges(0.6,2.5);
  TCanvas c1("c1","projections",0,0,1000,1000);
  c1.DivideSquare(waves.size());
  for (uint w = 0; w < waves.size(); w++) {
    for (uint t=0; t < Ninterp; t++) {
      double e = ranges.first + (ranges.second-ranges.first)/(Ninterp-1)*t;  /* triky */
      double s = e*e;
      double m23 = mIso[waves[w].S];  // 2*PI_MASS+(mIso-2*PI_MASS)*(1.-exp(-1./mIso*(e-3*PI_MASS)));
      double val = MAscoli::getProjectedReducedDeck(waves[w].J, waves[w].M, waves[w].L,
                                                    POW2(m23), waves[w].S, R,
                                                    s, tP,
                                                    mtRsq,
                                                    2*PROT_MASS*E_BEAM_LAB,
                                                    mAsq, mBsq, mDsq,
                                                    m1sq);
      ltable[t] = std::make_pair(e, val);
      // std::cout << ltable[t].first << ", " << ltable[t].second << "\n";
    }
    TGraph *gr = new TGraph(Ninterp);
    gr->SetTitle(TString::Format("Deck projection, %s %s#pi %c-wave; M_{3#pi}, GeV",
                                 waves[w].title.c_str(),
                                 (std::vector<std::string>{"(#pi#pi)_S", "#rho", "f_{2}"})[waves[w].S].c_str(),
                                 (std::vector<char>{'S', 'P', 'D', 'F', 'G', 'H'})[waves[w].L]));
    for (uint t=0; t < Ninterp; t++) {
      gr->GetX()[t] = ltable[t].first;
      gr->GetY()[t] = ltable[t].second;
    }
    c1.cd(w+1); gr->Draw("al");
  }
  c1.SaveAs("/tmp/waves.pdf");

  return 1;
}
