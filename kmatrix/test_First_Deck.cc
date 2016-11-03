// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <vector>
#include <fstream>
// #include <>

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "MAscoli.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {

  // Deck parameters
  uint J = 2;
  double mAsq = POW2(PI_MASS);
  double mBsq = POW2(PROT_MASS);
  double mDsq = POW2(PROT_MASS);
  double mtRsq = POW2(PI_MASS);
  double s = 2*190*PROT_MASS;
  double t = -0.01;
  double R = 5.;

  /************************************************************************/
  /* rho pi P - wave */
  double iso_mass = RHO_MASS;
  uint S = 1;
  uint L = 1;
  // calculate value of the function
  double e = 1.1;
  cd vd = MAscoli::getProjectedDeck(mAsq, mBsq, e*e, mDsq, mtRsq, POW2(iso_mass), s, t,
                                   J, 0, L, S, R);
  std::cout << "value is " << vd << "\n";

  return 0;
}
