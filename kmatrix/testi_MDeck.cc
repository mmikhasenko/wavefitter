// Copyright [11.2016] Misha Mikhasenko


#include <MDeck.h>
#include <MAscoli.h>
#include "constants.h"
#include "dFunction.hpp"


int main() {

  double z = 0.9759, phi = 1.107;
  double s1 = 5.0537;
  double R = 5;
  double costheta_pr = -0.434, phi_pr = -3.097;
  double s = 12.2836, t = -0.1;
  double stot = 2*PROT_MASS*190;
  double mAsq = POW2(PI_MASS), mBsq = POW2(PROT_MASS), mDsq = POW2(PROT_MASS);
  double m1sq = POW2(PI_MASS), m2sq = POW2(PI_MASS), m3sq = POW2(PI_MASS);

  std::cout << MDeck::getAmplitude(z, phi,
                                   s1, R,
                                   costheta_pr, phi_pr,
                                   s, t,
                                   POW2(PI_MASS),
                                   stot,
                                   mAsq, mBsq, mDsq,
                                   m1sq, m2sq, m3sq) << "\n";

  std::cout << "MAscoli::psi = " << MAscoli::psi(z,
                                                 s1,
                                                 s, t,
                                                 mAsq,
                                                 m1sq) << "\n";
  std::cout << "d_{1,0}^2(1.1) = " << rpwa::dFunction<double>(4, 2, -2, 1.1) << "\n";
  return 0.0;
}
