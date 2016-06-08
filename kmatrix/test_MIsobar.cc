// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "mstructures.hh"

int main(int argc, char *argv[]) {
  const int Np = 116;
  const double Hlim = 2.5*2.5;

  // stable
  MIsobar rho_s(0.77, 0.15, 0.14, 0.14, 0.14, 1, 5, false);
  MIsobar rho_q(0.77, 0.15, 0.14, 0.14, 0.14, 1, 5, true);
  rho_q.makeLookupTable();

  TCanvas c1("c1");
  draw([&](double e)->double{return rho_s.rho(e*e);}, 0., 2.5)->Draw("alp");
  draw([&](double e)->double{return rho_q.rho(e*e);}, 0., 2.5)->Draw("same");
  c1.SaveAs("/tmp/a.pdf");

  return 0;
}
