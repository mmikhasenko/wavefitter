// Copyright [2016] Misha Mikhasenko
//
// Code is partly copied from ROOTPWA/decayAmplitude/massDependence.cc
// written by Boris Grube
//

#include "MIsobarPiPiS.h"
#include "mintegrate.h"
#include "mathUtils.hpp"

namespace b = boost::numeric::ublas;

MIsobarPiPiS::MIsobarPiPiS() :
  MIsobar(1.1, 0.7, PI_MASS, PI_MASS, 0, 5.),
  _a(2, b::matrix<cd>(2,2)), _c(5, b::matrix<cd>(2,2)), _sP(1, 2) {

  // digitize table given in 
  _a[0](0, 0) =  0.1131;  // AMP Table 1, M solution: f_2^2
  _a[0](0, 1) =  0.0150;  // AMP Table 1, M solution: f_1^3
  _a[0](1, 0) =  0.0150;  // AMP Table 1, M solution: f_1^3
  _a[0](1, 1) = -0.3216;  // AMP Table 1, M solution: f_2^3

  const double f[2] = {0.1968, -0.0154};  // AMP Table 1, M solution: f_1^1 and f_2^1
  _a[1](0, 0) = f[0] * f[0];
  _a[1](0, 1) = f[0] * f[1];
  _a[1](1, 0) = f[1] * f[0];
  _a[1](1, 1) = f[1] * f[1];
  
  _c[0](0, 0) =  0.0337;                // AMP Table 1, M solution: c_11^0
  _c[1](0, 0) = -0.3185;                // AMP Table 1, M solution: c_11^1
  _c[2](0, 0) = -0.0942;                // AMP Table 1, M solution: c_11^2
  _c[3](0, 0) = -0.5927;                // AMP Table 1, M solution: c_11^3
  _c[4](0, 0) =  0.1957;                // AMP Table 1, M solution: c_11^4
  _c[0](0, 1) = _c[0](1, 0) = -0.2826;  // AMP Table 1, M solution: c_12^0
  _c[1](0, 1) = _c[1](1, 0) =  0.0918;  // AMP Table 1, M solution: c_12^1
  _c[2](0, 1) = _c[2](1, 0) =  0.1669;  // AMP Table 1, M solution: c_12^2
  _c[3](0, 1) = _c[3](1, 0) = -0.2082;  // AMP Table 1, M solution: c_12^3
  _c[4](0, 1) = _c[4](1, 0) = -0.1386;  // AMP Table 1, M solution: c_12^4
  _c[0](1, 1) =  0.3010;                // AMP Table 1, M solution: c_22^0
  _c[1](1, 1) = -0.5140;                // AMP Table 1, M solution: c_22^1
  _c[2](1, 1) =  0.1176;                // AMP Table 1, M solution: c_22^2
  _c[3](1, 1) =  0.5204;                // AMP Table 1, M solution: c_22^3
  _c[4](1, 1) = -0.3977;                // AMP Table 1, M solution: c_22^4
  
  _sP(0, 0) = -0.0074;  // AMP Table 1, M solution: s_0
  _sP(0, 1) =  0.9828;  // AMP Table 1, M solution: s_1

  // change parameters according to Kachaev's prescription
  _c[4](0, 0) = 0; // was 0.1957;
  _c[4](1, 1) = 0; // was -0.3977;
  
  _a[0](0, 1) = 0; // was 0.0150
  _a[0](1, 0) = 0; // was 0.0150
  
  _a[1](0, 0) = 0;
  _a[1](0, 1) = 0;
  _a[1](1, 0) = 0;
  _a[1](1, 1) = 0;

  // calculate isobar shape integral
  intU = integrate([&](double u)->double{
      return U(1./u)*1./(u*u)/(2*M_PI);
    }, 0, 1./sth());
  std::cout << "Integral over isobar shape is recalculated. Value is " << intU << "\n";

}

double MIsobarPiPiS::U(double s) const {
  const double cR = 1./35;
  return cR * norm(T(s)) * 1./(8*M_PI)*sqrt(LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS)))/s;
}

cd MIsobarPiPiS::U(cd s) const {
  const cd qPiPi   = sqrt( LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS))/(4.*s) );
  const cd qPi0Pi0 = sqrt( LAMBDA(s, POW2(PI0_MASS), POW2(PI0_MASS))/(4.*s) );

  cd rho00 = (qPiPi + qPi0Pi0) / sqrt(s);  // avarage 2q/m

  const cd scale = (s / (4 * POW2((K_MASS+K0_MASS)/2.))) - 1.;

  cd M00 = 0.0;
  for (uint i = 0; i < _sP.size2(); ++i) {
    const cd fa = 1. / (s - _sP(0, i));
    M00 += fa * _a[i](0,0);
  }
  for (uint i = 0; i < _c.size(); ++i) {
    const cd sc = pow(scale, (int)i);
    M00 += sc *_c[i](0,0);
  }

  double cR = 1./35.;
  cd _U00 = cR*POW2(16*M_PI)/(M00*M00 + rho00*rho00) * 1./(8*M_PI)*sqrt(LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS)))/s;
  return _U00;
}

cd MIsobarPiPiS::T(cd s) const {

  const cd qPiPi   = sqrtPi( LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS))/(4.*s) );
  const cd qPi0Pi0 = sqrtPi( LAMBDA(s, POW2(PI0_MASS), POW2(PI0_MASS))/(4.*s) );

  cd rho00 = (qPiPi + qPi0Pi0) / sqrt(s);  // avarage 2q/m

  const cd scale = (s / (4 * POW2((K_MASS+K0_MASS)/2.))) - 1.;

  cd M00 = 0.0;
  for (uint i = 0; i < _sP.size2(); ++i) {
    const cd fa = 1. / (s - _sP(0, i));
    M00 += fa * _a[i](0,0);
  }
  for (uint i = 0; i < _c.size(); ++i) {
    const cd sc = pow(scale, (int)i);
    M00 += sc *_c[i](0,0);
  }

  cd _T00 = 1./(M00 - cd(0, 1.)*rho00);
  return 16.*M_PI*_T00;
}

cd MIsobarPiPiS::T(double s) const {

  double mass = sqrt(s);
  if (fabs(s - _sP(0, 1)) < 1e-6) s  = POW2(sqrt(s)+1e-6);

  const cd qPiPi   = sqrt( LAMBDA(cd(s,0.), POW2(PI_MASS), POW2(PI_MASS))/(4.*s) );
  const cd qPi0Pi0 = sqrt( LAMBDA(cd(s,0.), POW2(PI0_MASS), POW2(PI0_MASS))/(4.*s) );
  const cd qKK     = sqrt( LAMBDA(cd(s,0.), POW2(K_MASS), POW2(K_MASS))/(4.*s) );
  const cd qK0K0   = sqrt( LAMBDA(cd(s,0.), POW2(K0_MASS), POW2(K0_MASS))/(4.*s) );

  b::matrix<cd> rho(2, 2);
  rho(0, 0) = ((2. * qPiPi) / mass + (2. * qPi0Pi0) / mass) / 2.;
  rho(1, 1) = ((2. * qKK)   / mass + (2. * qK0K0)   / mass) / 2.;
  rho(0, 1) = rho(1, 0) = 0;

  const double scale = (s / (4 * POW2((K_MASS+K0_MASS)/2.))) - 1.;

  b::matrix<cd> M(b::zero_matrix<cd>(2, 2));
  for (uint i = 0; i < _sP.size2(); ++i) {
    const cd fa = 1. / (s - _sP(0, i));
    M += fa * _a[i];
  }
  for (uint i = 0; i < _c.size(); ++i) {
    const cd sc = pow(scale, (int)i);
    M += sc *_c[i];
  }

  // modification: off-diagonal terms set to 0
  M(0, 1) = 0;
  M(1, 0) = 0;

  b::matrix<cd> _T(2, 2);
  rpwa::invertMatrix<cd>(M - cd(0, 1.) * rho, _T);

  const cd amp = _T(0, 0);
  return 16.*M_PI*amp;
}
