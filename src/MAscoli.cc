// Copyright [08.08.2016] Misha Mikhasenko

#include "MAscoli.h"

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "Math/SpecFuncMathMore.h"

#include "constants.h"
#include "deflib.h"
// #include "mstructures.h"
#include "dFunction.hpp"

#include "mintegrate.h"

double MAscoli::getDeck(double mAsq, double mBsq, double wsq, double mDsq, double mtRsq,
                        double mIsq, double s, double t, double z,
                        uint S, int lamS) {

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
  double val = ((S%2==0) ? 1. : -1.)*  // (-1)^S
    rpwa::dFunction<double>(2*S, 0, 2*lamS, atan2(sqrt(mIsq)*pa*sqrt(1-z*z), eI*pa*z-p1*ea))*
    1./(mtRsq - tR) *  // pion propagator
    sqrt(2*S+1) *  // kind of normalisation
    spip_int_phi;  // nuclear vertex
  return val;
}

double MAscoli::getProjectedDeck(double mAsq, double mBsq, double wsq, double mDsq, double mtRsq,
                                 double mIsq, double s, double t,
                                 uint J, uint M, uint L, uint S) {
  if (M!=0) { std::cerr << "Error<>: M!=0 not implemented!" << std::endl; return 1.; }
  // integrate
  int minSJ = (S < J) ? S : J;
  double int_val =
    integrate([&, mAsq, mBsq, wsq, mDsq, mtRsq, mIsq, s, J, L, S, minSJ](double z)->double{
      double val = 0;
      for (int lamS = -minSJ; lamS <= minSJ; lamS++)
        val +=
          ROOT::Math::wigner_3j(2*L, 2*S, 2*J, 2*0, 2*lamS, -2*lamS) *
          getDeck(mAsq, mBsq, wsq, mDsq, mtRsq, mIsq, s, t, z, S, -lamS) *
          rpwa::dFunction<double> (2*J, 2*M, 2*lamS, acos(z));
      return val;
      }, -1., 1.);
  int_val *= sqrt(2*J+1) * /*because of Clebsch and 3j*/
    sqrt((2*L+1)); /*because of normalisation*/

  return int_val;
}
