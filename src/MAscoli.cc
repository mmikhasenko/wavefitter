// Copyright [08.08.2016] Misha Mikhasenko

#include "MAscoli.h"

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "Math/SpecFuncMathMore.h"

#include "constants.h"
// #include "mstructures.h"
#include "dFunction.hpp"

#include "mintegrate.h"

double MAscoli::getDeck(double mAsq, double mBsq, double wsq, double mDsq, double mtRsq,
                        double mIsq, double s, double t, double z,
                        uint S, int lamS, double R) {
  if (sqrt(wsq) <= sqrt(mIsq)+PI_MASS) return 0;
  
  double pa = sqrt(LAMBDA(wsq, mAsq, t)/(4*wsq));
  double pd = sqrt(LAMBDA(wsq, mDsq, s)/(4*wsq));
  double u = mAsq + wsq + mBsq + mDsq - s - t;
  double pb = sqrt(LAMBDA(wsq, mDsq, u)/(4*wsq));
  double ea = sqrt(pa*pa + mAsq);
  double ed = sqrt(pd*pd + mDsq);
  // double eb = sqrt(pb*pb + mBsq);  // not used
  double p1 = sqrt(LAMBDA(wsq, mIsq, POW2(PI_MASS))/(4*wsq));
  double e1 = sqrt(p1*p1 + POW2(PI_MASS));
  double eI = sqrt(wsq) - e1;
  double cos_epsilon = (pb*pb - pd*pd - pa*pa)/(2.*pd*pa);
  double tR = mAsq + mIsq - 2.*ea*eI + 2.*pa*p1*z;
  double spip_int_phi = 2.*M_PI*(POW2(PI_MASS) + mDsq + 2.*ed*e1 - 2.*pd*p1*cos_epsilon*z);

  double qdsq = LAMBDA(mIsq, mAsq, tR)/(4*mIsq);
  // double R = 5;
  double damp = R*R*qdsq/(1.+R*R*qdsq);
  double val = ((S%2 == 0) ? 1. : -1.)*  // (-1)^S
    rpwa::dFunction<double>(2*S, 0, 2*lamS, atan2(sqrt(mIsq)*pa*sqrt(1-z*z), eI*pa*z-p1*ea))*
    1./(mtRsq - tR) *  // pion propagator
    sqrt(2*S+1) *  // kind of normalisation
    spip_int_phi *
    pow(damp, S/2.);  // nuclear vertex
  return val;
}

double MAscoli::getProjectedDeck(double mAsq, double mBsq, double wsq, double mDsq, double mtRsq,
                                 double mIsq, double s, double t,
                                 uint J, uint M, uint L, uint S, double R) {
  if (M!=0) { std::cerr << "Error<>: M!=0 not implemented!" << std::endl; return 1.; }
  // integrate
  int minSJ = (S < J) ? S : J;
  double int_val =
    integrate([&, mAsq, mBsq, wsq, mDsq, mtRsq, mIsq, s, J, L, S, minSJ, R](double z)->double{
      double val = 0;
      for (int lamS = -minSJ; lamS <= minSJ; lamS++)
        val +=
          ROOT::Math::wigner_3j(2*L, 2*S, 2*J, 2*0, 2*lamS, -2*lamS) *
          getDeck(mAsq, mBsq, wsq, mDsq, mtRsq, mIsq, s, t, z, S, -lamS, R) *
          rpwa::dFunction<double> (2*J, 2*M, 2*lamS, acos(z));
      return val;
      }, -1., 1.);
  int_val *= sqrt(2*J+1) * /*because of Clebsch and 3j*/
    sqrt((2*L+1)); /*because of normalisation*/

  return int_val;
}

cd MAscoli::fullDeckTerm(double costheta, double phi,
                         double mS1sq, uint S1, int lamS1, double RS1,
                         double costheta_pr, double phi_pr,
                         double wsq, double t,
                         double mtRsq,
                         double stot,
                         double mAsq, double mBsq, double mDsq,
                         double m1sq) {
  if (sqrt(wsq) <= sqrt(mS1sq)+sqrt(m1sq)) return 0;
  // calculate from three component
  cd one_term =
    sPionProton(costheta, phi, mS1sq, wsq, t, stot, mAsq, mBsq, mDsq, m1sq) *
    upperPart(costheta, mS1sq, S1, lamS1, RS1, wsq, t, mtRsq, mAsq, m1sq) *
    isobarDecayAnglesTerm(costheta_pr, phi_pr, S1, lamS1);
  return one_term;
}


// to calculate M = 0, you can just put phi = M_PI/2 and get integration factor 2*M_PI
double MAscoli::sPionProton(double costheta, double phi,
                           double mS1sq,
                           double wsq, double t,
                           double stot,
                           double mAsq, double mBsq, double mDsq,
                           double m1sq) {
  double pa = sqrt(LAMBDA(wsq, mAsq, t)/(4*wsq));  // GJ
  double pd = sqrt(LAMBDA(wsq, mDsq, stot)/(4*wsq));  // GJ
  double u = mAsq + wsq + mBsq + mDsq - stot - t;     // GJ
  double pb = sqrt(LAMBDA(wsq, mDsq, u)/(4*wsq));  // GJ
  double ed = sqrt(pd*pd + mDsq);                  // GJ
  double p1 = sqrt(LAMBDA(wsq, mS1sq, m1sq)/(4*wsq));  // break up at GJ
  double e1 = sqrt(p1*p1 + m1sq);  // E_{pion1} at GJ
  // pion-proton vertex
  double cos_epsilon = (pb*pb - pd*pd - pa*pa)/(2.*pd*pa);  // some epsilon angle for pion-proton vertex
  // Warning: the expression has not been checked!
  double epsilon = acos(cos_epsilon);  // angle in triagle
  double spip_int_phi = m1sq + mDsq + 2.*ed*e1 - 2.*pd*p1*
    (cos_epsilon*costheta + sin(epsilon)*sqrt(1-costheta*costheta)*cos(phi));

  return spip_int_phi;
}

cd MAscoli::isobarDecayAnglesTerm(double costheta_pr, double phi_pr, int S1, int lamS1) {
  return rpwa::DFunction<cd>(2*S1, 2*lamS1, 0, phi_pr, acos(costheta_pr), 0.0, false);
}

double MAscoli::upperPart(double costheta,
                          double mS1sq, uint S1, int lamS1, double RS1,
                          double wsq, double t,
                          double mtRsq,
                          double mAsq,
                          double m1sq) {
  if (sqrt(wsq) <= sqrt(mS1sq)+sqrt(m1sq)) return 0;

  double pa = sqrt(LAMBDA(wsq, mAsq, t)/(4*wsq));  // GJ
  double p1 = sqrt(LAMBDA(wsq, mS1sq, m1sq)/(4*wsq));  // break up at GJ
  double e1 = sqrt(p1*p1 + m1sq);  // E_{pion1} at GJ
  double eI = sqrt(wsq) - e1;  // E_{isobar} at GJ
  double ea = sqrt(pa*pa + mAsq);
  // virtuality
  double tR = mAsq + mS1sq - 2.*ea*eI + 2.*pa*p1*costheta;  // GJ
  // left break-up momentum for damping
  double qdsq = LAMBDA(mS1sq, mAsq, tR)/(4*mS1sq);  // tR OR mtRsq?
  double damp = RS1*RS1*qdsq/(1.+RS1*RS1*qdsq);
  // calculation for the main quantity
  double psi = atan2(sqrt(mS1sq)*pa*sqrt(1-costheta*costheta), eI*pa*costheta-p1*ea);

  double val =
    rpwa::dFunction<double>(2*S1, 0, 2*lamS1, psi)*  // strange angle phi
    1./(mtRsq - tR) *  // pion propagator
    pow(damp, S1/2.);  // left damping
  return val;
}

