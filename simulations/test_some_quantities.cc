// [2016] Misha Mikhasenko
// Discription: for tests

#include "MIsobar.h"
#include "MAscoli.h"
#include "MDeck.h"

#define E_BEAM 16

int main(int argc, char const *argv[]) {
  // masses
  double mpisq = POW2(PI_MASS), mpsq = POW2(PROT_MASS);
  double m1sq = mpisq; //, m2sq = mpisq, m3sq = mpisq;
  double mAsq = mpisq;
  double mBsq = mpsq, mDsq = mpsq;
  // static invariants
  double s0 = 2 * PROT_MASS * E_BEAM + POW2(PROT_MASS) + POW2(PI_MASS);
  double tpr = -0.05;
  // dinamical
  double costheta = 0.2;
  double phi = 0.1;
  double s1 = POW2(0.820);
  double wsq = POW2(1.1);
  // tmin
  double tmin = MAscoli::tmin(mAsq, mBsq, wsq, mDsq, s0);
  std::cout << "tmin " << tmin << '\n';
  double t = tpr + tmin;
  // psi
  double psi = MAscoli::psi(costheta, s1, wsq, t, mAsq, m1sq);
  std::cout << "psi " << psi << '\n';
  // sPiN
  double sPiN = MAscoli::sPionProton(costheta, phi, s1, wsq, t, s0, mAsq, mBsq,
                                     mDsq, m1sq);
  std::cout << "sPiN " << sPiN << '\n';
  return 0;
}
