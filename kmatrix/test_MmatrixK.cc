// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MmatrixK.h"
#include "MIsobarChannel.h"
#include "MParKeeper.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {

  // stable
  MIsobar rho(0.77, 0.15, 0.14, 0.14, 1, 5.);
  MIsobar  f2(1.23, 0.2,  0.14, 0.14, 2, 5.);
  MIsobarChannel rho_q(rho, 0.14);
  MIsobarChannel  f2_q(f2,  0.14);

  MmatrixK km({&rho_q, &f2_q}, 2);
  // masses
  MParKeeper::gI()->set("m1", 1.502);
  MParKeeper::gI()->set("m2", 1.701);
  // couplings
  MParKeeper::gI()->set("g1", 4.);
  MParKeeper::gI()->set("g2", 1.);
  MParKeeper::gI()->set("h1", 1.);
  MParKeeper::gI()->set("h2", 4.);
  MParKeeper::gI()->printAll();

  TCanvas c1("c1");
  combine(
          draw([&](double s)->double{ b::matrix<cd> m = km.getValue(s); return norm(m(0, 0)+m(0, 1));},
               rho.sth(), 5.1, 200),
          draw([&](double s)->double{ b::matrix<cd> m = km.getValue(s); return norm(m(1, 0)+m(1, 1));},
               f2.sth(), 5.1, 200) )->Draw("apl");

  c1.SaveAs("/tmp/a.pdf");

  return 0;
}
