// Copyright [11.2016] Misha Mikhasenko


#include <MDeck.h>
#include <MAscoli.h>
#include "constants.h"

#include "TWigner.h"
#include "Math/SpecFuncMathMore.h"


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
  std::cout << "d_{1,0}^2(1.1) = " << Math::WignerD(7, 1, -3, 1.1) << "\n";


  std::cout << "proj = " << MAscoli::getDeck(POW2(PI_MASS), POW2(PROT_MASS), 3.0, POW2(PROT_MASS), POW2(PI_MASS),
                                                       2.0, 2.*PROT_MASS*190, -0.1, 0.5,
                                                       7, -4, 5.0) << "\n";
  std::cout << "clebsch = " << sqrt(2*2+1)*ROOT::Math::wigner_3j(2*3, 2*1, 2*2, 2*0, 2*1, -2*1) << "\n";
  return 0.0;

  // generate output tree
//  TRandom3 ran(12345);
//  TFile fout("/tmp/1.out", "recreate");
//  TTree tout("tout", "events");
//  const uint Nmc = 1e7;
//  for (uint e = 0; e < Nmc; e++) {
//    z = 2*ran.Rndm()-1.;
//    s = 
//  }

}
