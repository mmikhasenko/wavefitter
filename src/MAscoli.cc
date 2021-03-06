// Copyright [08.08.2016] Misha Mikhasenko

#include "MAscoli.h"

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "Math/SpecFuncMathMore.h"

#include "constants.h"
// #include "mstructures.h"
#include "TWigner.h"

#include "mintegrate.h"

//***********************************************************
// For the reduced amplitude
// a factor sqrt((2S+1)/(4#pi)) D_{#lambda 0}^{S*}(#theta') is removed
//***********************************************************
double MAscoli::getReducedDeck(double costheta, double phi, double mS1sq,
                               uint S1, int lamS1, double RS1, double wsq,
                               double t, double mtRsq, double stot, double mAsq,
                               double mBsq, double mDsq, double m1sq) {
  if (sqrt(wsq) <= sqrt(mS1sq) + sqrt(m1sq))
    return 0;

  // left break-up momentum for damping
  // double qdsq = LAMBDA(mS1sq, mAsq, mtRsq) / (4 * mS1sq);  // tR OR mtRsq?
  // double qdsq0 = LAMBDA(mS1sq, wPOW2(PI_MASS), POW2(PI_MASS)) / (4 * mS1sq);  // tR OR mtRsq?
  // double ratio =  (1./(RS1 * RS1 * qdsq) + 1.) / (1./(RS1 * RS1 * qdsq0) + 1.);

  // virtuality
  double tR = tPionIsobar(costheta, mS1sq, wsq, t, mAsq, m1sq); // GJ

  double spip =
      sPionProton(costheta, phi, mS1sq, wsq, t, stot, mAsq, mBsq, mDsq, m1sq);
  double amp_pion_proton = spip + (t - mBsq - mDsq - tR - m1sq) / 2;

  double val =
      amp_pion_proton *
      upperPart(costheta, mS1sq, S1, lamS1, RS1, wsq, t, mtRsq, mAsq, m1sq) *
      4 * M_PI * sqrt((2. * S1 + 1.) / (4 * M_PI));  // *  // left after reduction
      // pow(ratio, S1 / 2.);                             // left damping
  return val;
}

double MAscoli::tmin(double m1sq, double m2sq, double m3sq, double m4sq, double s) {
  return m1sq + m3sq -
    2 * (s+m1sq - m2sq) * (s+m3sq-m4sq) / (4*s) +
    2 * sqrt(LAMBDA(s,m1sq,m2sq)*LAMBDA(s,m3sq,m4sq)) / (4*s);
}

//***********************************************************
// projection is done in a standart way
// ((2L+1)/(2J+1))^{1/2} #sum <L0Sl|Jl> ((2J+1)/(4#pi))^{1/2}
// D_{Ml}^{J*}(#theta)
//***********************************************************
double MAscoli::getProjectedReducedDeck(uint J, int M, uint L, double mS1sq,
                                        uint S1, double RS1, double wsq,
                                        double t, double mtRsq, double stot,
                                        double mAsq, double mBsq, double mDsq,
                                        double m1sq) {
  if (M > 1 || M < -1)
    return 0;
  // integrate
  int minSJ = (S1 < J) ? S1 : J;
  double int_val =
      integrate([=](double z)->double {
                  double val = 0;
                  for (int lamS = -minSJ; lamS <= minSJ; lamS++) {
                    double reduced_deck = getReducedDeckProjectionM(
                        M, z, mS1sq, S1, lamS, RS1, wsq, t, mtRsq,
                        stot, mAsq, mBsq, mDsq, m1sq);
                    // calculate all together
                    val +=
                        // clebsch coefficient
                        sqrt(2 * J + 1) *
                        (((L - S1 + (-lamS)) % 2 == 1) ? (-1) : (1)) *
                        ROOT::Math::wigner_3j(2 * L, 2 * S1, 2 * J, 2 * 0,
                                              2 * lamS, -2 * lamS) *
                        // reduced Deck where integral over M has been performed
                        reduced_deck *
                        // projecting D-function
                        sqrt((2. * J + 1.) / (4. * M_PI)) *
                        Math::WignerD(2 * J, 2 * M, 2 * lamS, acos(z));
                  }
                  return val;
                },
                -1., 1.);
  int_val *= sqrt((2. * L + 1.) / (2. * J + 1.)); /*because of normalisation*/

  return int_val;
}

double MAscoli::getProjectedReducedDeck(uint J, int M, bool pos_refl, uint L,
                                        double mS1sq, uint S1, double RS1,
                                        double wsq, double t, double mtRsq,
                                        double stot, double mAsq, double mBsq,
                                        double mDsq, double m1sq) {
  if (M < 0)
    return 0.;
  int eps = pos_refl ? 1 : -1;
  int m1JM = (J + M) % 2 ? -1 : 1;
  double factor = eps * m1JM;
  double val1 = getProjectedReducedDeck(J, M, L, mS1sq, S1, RS1, wsq, t, mtRsq,
                                       stot, mAsq, mBsq, mDsq, m1sq);
  double val2 = getProjectedReducedDeck(J, -M, L, mS1sq, S1, RS1, wsq, t, mtRsq,
                                       stot, mAsq, mBsq, mDsq, m1sq);
  return  (val1 - factor*val2) / (M == 0 ? 2. : sqrt(2.));
}

//***********************************************************
// projection is done in a standart way
// ((2L+1)/(2J+1))^{1/2} #sum <L0Sl|Jl> ((2J+1)/(4#pi))^{1/2}
// D_{Ml}^{J*}(#theta)
//***********************************************************
double MAscoli::getProjectionJMSlam(uint J, int M, int lam, double mS1sq,
                                    uint S1, double RS1, double wsq, double t,
                                    double mtRsq, double stot, double mAsq,
                                    double mBsq, double mDsq, double m1sq) {
  if (M > 1 || M < -1)
    return 0;
  // integrate
  if (abs(lam) > J || abs(lam) > S1) return 0.0;
  double int_val =
      integrate([=](double z)->double {
                  double val = 0;
                  double reduced_deck = getReducedDeckProjectionM(
                        M, z, mS1sq, S1, lam, RS1, wsq, t, mtRsq,
                        stot, mAsq, mBsq, mDsq, m1sq);
                    // calculate all together
                    val +=
                        // reduced Deck where integral over M has been performed
                        reduced_deck *
                        // projecting D-function
                        sqrt((2. * J + 1.) / (4. * M_PI)) *
                        Math::WignerD(2 * J, 2 * M, 2 * lam, acos(z));
                  return val;
                },
                -1., 1.);
  return int_val;
}

//***********************************************************
// perform \int e^{-i M Phi} B(Phi) \diff Phi
// assuming B(phi) = A + B cos Phi
//***********************************************************
double MAscoli::getReducedDeckProjectionM(int M, double costheta,
                                          double mS1sq, uint S1, int lamS1,
                                          double RS1, double wsq, double t,
                                          double mtRsq, double stot,
                                          double mAsq, double mBsq, double mDsq,
                                          double m1sq) {
  // a bit tricky with the reduction because of phi-dependence
  // I know that is in the form (A + B cos Phi)
  // the projection of phi is given by  spip = (A + B cos Phi).
  // \int_{0}^{2pi} \diff Phi exp(-i M Phi) (A + B cos Phi),
  // it gives
  //   -- M = 0 :  2pi A
  //   -- M = \pm 1: pi B = pi ( spip(0) - spip(pi/2) )
  double reduced_deck = 0;
  if (M == 0) {
    reduced_deck =
        2 * M_PI * getReducedDeck(costheta, M_PI / 2., mS1sq, S1, lamS1, RS1,
                                  wsq, t, mtRsq, stot, mAsq, mBsq, mDsq, m1sq);
  } else {
    reduced_deck =
        M_PI * (getReducedDeck(costheta, 0., mS1sq, S1, lamS1, RS1, wsq, t,
                               mtRsq, stot, mAsq, mBsq, mDsq, m1sq) -
                getReducedDeck(costheta, M_PI / 2., mS1sq, S1, lamS1, RS1, wsq,
                               t, mtRsq, stot, mAsq, mBsq, mDsq, m1sq));
  }
  return reduced_deck;
}

//***********************************************************
// still the factor (2S+1) t_S is not in
// one has to sum over S with the factor and unitarity amplitude
//***********************************************************
cd MAscoli::fullDeckTerm(double costheta, double phi, double mS1sq, uint S1,
                         int lamS1, double RS1, double costheta_pr,
                         double phi_pr, double wsq, double t, double mtRsq,
                         double stot, double mAsq, double mBsq, double mDsq,
                         double m1sq) {
  if (sqrt(wsq) <= sqrt(mS1sq) + sqrt(m1sq))
    return 0;
  // calculate pi-p amplitude
  double tR = tPionIsobar(costheta, mS1sq, wsq, t, mAsq, m1sq); // GJ
  // pi-p
  double spip =
      sPionProton(costheta, phi, mS1sq, wsq, t, stot, mAsq, mBsq, mDsq, m1sq);
  double amp_pion_proton = spip + (t - mBsq - mDsq - tR - m1sq) / 2;
  // calculate from three component
  cd one_term =
      amp_pion_proton *
      upperPart(costheta, mS1sq, S1, lamS1, RS1, wsq, t, mtRsq, mAsq, m1sq) *
      Math::WignerD(2 * S1, 2 * lamS1, 0, -phi_pr, acos(costheta_pr), 0.0);  // -phi to conjugate
  return one_term;
}

//***********************************************************
// the quantity depends on phi in a simple way, ~const + ~cos(phi)
// to calculate M = 0 you can just put phi = #pi/2 and get integration factor
// 2*M_PI
// for M = 1 one has to drop the first term, put phi = 0 and multiply to i#pi
//   likely there is a trick: put costheta to 0 and multiply back to sign latter
//***********************************************************
double MAscoli::sPionProton(double costheta, double phi, double mS1sq,
                            double wsq, double t, double stot, double mAsq,
                            double mBsq, double mDsq, double m1sq) {
  if (sqrt(wsq) > sqrt(stot) - sqrt(mDsq)) {
    std::cerr << "Error<MAscoli::sPionProton>: requested values are outside of "
                 "the physical region."
              << "pd is imaginary.\n";
  }
  double pa = sqrt(LAMBDA(wsq, mAsq, t) / (4 * wsq));     // GJ
  double pd = sqrt(LAMBDA(wsq, mDsq, stot) / (4 * wsq));  // GJ
  double u = mAsq + wsq + mBsq + mDsq - stot - t;         // GJ
  double pb = sqrt(LAMBDA(wsq, mDsq, u) / (4 * wsq));     // GJ
  double ed = sqrt(pd * pd + mDsq);                       // GJ
  double p1 = sqrt(LAMBDA(wsq, mS1sq, m1sq) / (4 * wsq)); // break up at GJ
  double e1 = sqrt(p1 * p1 + m1sq);                       // E_{pion1} at GJ
  // pion-proton vertex
  double cos_epsilon =
      (pb * pb - pd * pd - pa * pa) /
      (2. * pd * pa); // some epsilon angle for pion-proton vertex
  if (fabs(cos_epsilon) > 1) {
    std::cerr << "Error<MAscoli::sPionProton>: |cos_epsilon(" << cos_epsilon
              << ")|>1\n";
    return 0.0;
  }
  // Warning: the expression has not been checked!
  double epsilon = acos(cos_epsilon); // angle in triagle
  double spip_int_phi =
      m1sq + mDsq + 2. * ed * e1 -
      2. * pd * p1 * (cos_epsilon * costheta +
                      sin(epsilon) * sqrt(1 - costheta * costheta) * cos(phi));

  return spip_int_phi;
}

double MAscoli::tPionIsobar(double costheta, double mS1sq, double wsq, double t,
                            double mAsq, double m1sq) {
  double pa = sqrt(LAMBDA(wsq, mAsq, t) / (4 * wsq));     // GJ
  double p1 = sqrt(LAMBDA(wsq, mS1sq, m1sq) / (4 * wsq)); // break up at GJ
  double eI = sqrt(p1 * p1 + mS1sq);                      // E_{pion1} at GJ
  double ea = sqrt(pa * pa + mAsq);
  return mAsq + mS1sq - 2. * ea * eI + 2. * pa * p1 * costheta;
}

double MAscoli::upperPart(double costheta, double mS1sq, uint S1, int lamS1,
                          double RS1, double wsq, double t, double mtRsq,
                          double mAsq, double m1sq) {
  if (sqrt(wsq) <= sqrt(mS1sq) + sqrt(m1sq))
    return 0;

  // virtuality
  double tR = MAscoli::tPionIsobar(costheta, mS1sq, wsq, t, mAsq, m1sq);
  // calculation for the main quantity
  double psi = MAscoli::psi(costheta, mS1sq, wsq, t, mAsq, m1sq);

  double val = Math::WignerD(2 * S1, 0, 2 * lamS1, psi) * // strange angle phi
               1. / (mtRsq - tR); // *  // pion propagator
  return val;
}

double MAscoli::psi(double costheta, double mS1sq, double wsq, double t,
                    double mAsq, double m1sq) {
  if (sqrt(wsq) <= sqrt(mS1sq) + sqrt(m1sq))
    return 0;

  double pa = sqrt(LAMBDA(wsq, mAsq, t) / (4 * wsq));     // GJ
  double p1 = sqrt(LAMBDA(wsq, mS1sq, m1sq) / (4 * wsq)); // break up at GJ
  double e1 = sqrt(p1 * p1 + m1sq);                       // E_{pion1} at GJ
  double eI = sqrt(wsq) - e1;                             // E_{isobar} at GJ
  double ea = sqrt(pa * pa + mAsq);
  double psi = atan2(sqrt(mS1sq) * pa * sqrt(1 - costheta * costheta),
                     eI * pa * costheta - p1 * ea);
  return psi;
}
