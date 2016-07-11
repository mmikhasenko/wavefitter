// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MmatrixK.h"
#include "MIsobarChannel.h"
#include "MProductionPhysics.h"
#include "MTwoBodyChannel.h"
#include "MParKeeper.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {
  // stable
  MIsobar rho(0.77, 0.15, PI_MASS, PI_MASS, 1, 5.);
  // MIsobarChannel rho_q(rho, PI_MASS);
  MTwoBodyChannel rho_pi(rho.GetM(), PI_MASS);

  rho_pi.makeDisperseLookupTable(0, 10., 100);

  std::vector<MChannel*> iset = {&rho_pi};
  MmatrixK km(iset, 1);
  // set parameters
  MParKeeper::gI()->set("m0", 1.502);
  // couplings
  MParKeeper::gI()->set("g0", 4.);

  // production
  MProductionPhysics pr(iset);
  pr.addScattering([&](double s)->b::matrix<cd>{return km.getValue(s);});
  // pr.addShortRange();
  // MParKeeper::gI()->set("c0", 1.);

  double sth0 = rho_pi.sth();
  std::vector<std::function<cd(double)> > getB(1);
  getB[0] = [sth0](double s)->cd{return pow(s-sth0, 2)*exp(-2.*s);};
  pr.addLongRange(getB);
  MParKeeper::gI()->set("Br", 1.);
  MParKeeper::gI()->set("Bi", 0.);
  pr.unitarize(30., 300);

  TCanvas c1("c1");
  combine(
          SET1(draw([&](double s)->double{ auto v = pr.getValue(s); return real(v(0));},
                    sth0, 5.1, 400),
               SetLineColor(kBlack)),
          SET1(draw([&](double s)->double{ auto v = pr.getValue(s); return imag(v(0));},
                    sth0, 5.1, 400),
               SetLineColor(kRed)))
    ->Draw("al");

  c1.SaveAs("/tmp/a.test_Unitarisation.pdf");

  return 0;
}
