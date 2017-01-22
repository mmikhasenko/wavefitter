// Copyright Misha Mikhasenko [12.2016], test
// test_ProjectedDeck.2:
//   In this example I demonstrate that Deck had physical boundary.
//   The amplitude can not be calculate for s > 100 GeV^2


#include "MAscoli.h"
#include "TCanvas.h"
#include "TGraph.h"

#include <utility>

#include "constants.h"
#include "mintegrate.h"
#include "mstructures.h"

#include "MIsobar.h"
#include "MIsobarChannel.h"


int main() {
  uint Jsector = 2;
  uint M = 0;
  uint L = 0;
  uint S = 2;
  double mIso = F2_MASS;
  double R = 5.0;
  double RI = 5.0;
  double E_BEAM_LAB = 190;
  double mAsq = POW2(PI_MASS);
  double mBsq = POW2(PROT_MASS);
  double mDsq = POW2(PROT_MASS);
  double m1sq = POW2(PI_MASS);
  double mtRsq = POW2(PI_MASS);
  double tP = -0.01;

  // isobar
  MIsobar f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, S, RI);
  MIsobarChannel *mCh = new MIsobarChannel(f2_iso, PI_MASS, L);
  MIsobarChannel *_iso = mCh;
  double sth = _iso->sth();

  // Deck table
  uint Ninterp = 100;
  std::vector<std::pair<double, double> > ltable(Ninterp);
  std::pair<double, double> ranges(5, 6.1);
  for (uint t=0; t < Ninterp; t++) {
    double e = ranges.first + (ranges.second-ranges.first)/(Ninterp-1)*t;  /* triky */
    double s = e*e;
    double m23 = 2*PI_MASS+(mIso-2*PI_MASS)*(1.-exp(-1./mIso*(e-3*PI_MASS)));
    double val = MAscoli::getProjectedReducedDeck(Jsector, M, L,
                                                  POW2(m23), S, R,
                                                  s, tP,
                                                  mtRsq,
                                                  2*PROT_MASS*E_BEAM_LAB,
                                                  mAsq, mBsq, mDsq,
                                                  m1sq);
    ltable[t] = std::make_pair(e, val);
    std::cout << ltable[t].first << ", " << ltable[t].second << "\n";
  }

  TGraph gr0(Ninterp); gr0.SetTitle("Deck projection. Border of the physical region; M_{3\pi}, GeV");
  for (uint t=0; t < Ninterp; t++) {
    gr0.GetX()[t] = ltable[t].first;
    gr0.GetY()[t] = ltable[t].second;
  }
  TCanvas c1("c1");
  gr0.Draw("al");
  c1.SaveAs("/tmp/f2PiS.just_deck.pdf");

  return 1;
}
