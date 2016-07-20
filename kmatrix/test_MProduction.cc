// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MmatrixK.h"
#include "MIsobarChannel.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {
  // stable
  MIsobar rho(0.77, 0.15, 0.14, 0.14, 1, 5.);
  MIsobar  f2(1.23, 0.2,  0.14, 0.14, 2, 5.);
  MIsobarChannel rho_q(rho, 0.14);
  MIsobarChannel  f2_q(f2,  0.14);
  rho_q.makeLookupTable(rho_q.sth(), 10., 100);
  f2_q.makeLookupTable(rho_q.sth(), 10., 100);

  rho_q.makeDisperseLookupTable(0, 10., 100);
  f2_q .makeDisperseLookupTable(0, 10., 100);

  std::vector<MChannel*> iset = {&rho_q, &f2_q};
  MmatrixK km(iset, 2);
  MProductionPhysics pr(iset);
  pr.addScattering([&](double s)->b::matrix<cd>{return km.getValue(s);});
  pr.addShortRange();
  // masses
  MParKeeper::gI()->set("m0", 1.502);
  MParKeeper::gI()->set("m1", 1.701);
  // couplings
  MParKeeper::gI()->set("g0", 4.);
  MParKeeper::gI()->set("g1", 1.);
  MParKeeper::gI()->set("h0", 1.);
  MParKeeper::gI()->set("h1", 4.);
  MParKeeper::gI()->set("cr0", 1.);
  MParKeeper::gI()->set("ci0", 0.);
  MParKeeper::gI()->set("cr1", 1.);
  MParKeeper::gI()->set("ci1", 0.);
  MParKeeper::gI()->printAll();

  TCanvas c1("c1");
  combine(
          SET2(draw([&](double s)->double{ auto v = pr.getValue(s); return norm(v(0));},
                    rho_q.sth(), 5.1, 400),
               SetLineColor(kRed), SetLineStyle(2)),
          SET2(draw([&](double s)->double{ auto v = pr.getValue(s); return norm(v(1));},
                    f2_q.sth(), 5.1, 400),
               SetLineColor(kBlack), SetLineStyle(1) ) )
    ->Draw("al");

  c1.SaveAs("/tmp/a.test_MProduction.pdf");

  return 0;
}
