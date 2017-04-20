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
  _a(2), _c(5), _sP(1, 2) {

  _a[0] =  0.1131;  // AMP Table 1, M solution: f_2^2

  const double f[2] = {0.1968, -0.0154};  // AMP Table 1, M solution: f_1^1 and f_2^1
  _a[1] = f[0] * f[0];
  
  _c[0] =  0.0337;                // AMP Table 1, M solution: c_11^0
  _c[1] = -0.3185;                // AMP Table 1, M solution: c_11^1
  _c[2] = -0.0942;                // AMP Table 1, M solution: c_11^2
  _c[3] = -0.5927;                // AMP Table 1, M solution: c_11^3
  _c[4] =  0.0;  // ?             // Kachaev modification

  _sP[0] = -0.0074;  // AMP Table 1, M solution: s_0
  _sP[1] =  0.9828;  // AMP Table 1, M solution: s_1

  // calculate isobar shape integral
  intU = 1.;
}

double MIsobarPiPiS::U(double s) const {
  return norm(T(s)) * 1./(8*M_PI)*sqrt(LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS)))/s / intU;
}

cd MIsobarPiPiS::U(cd s) const {
  cd rho00 = sqrt(LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS)))/s;
  const cd scale = (s / (4 * POW2((K_MASS+K0_MASS)/2.))) - 1.;

  cd M00 = 0.0;
  for (uint i = 0; i < _sP.size(); i++) M00 += _a[i] / (s - _sP[i]);
  for (uint i = 0; i < _c.size(); ++i)   M00 += pow(scale, static_cast<int>(i)) *_c[i];

  cd _U00 = POW2(16*M_PI) / (M00*M00 + rho00*rho00) * 1./(8*M_PI)*sqrt(LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS)))/s;
  return _U00 / intU;
}

cd MIsobarPiPiS::ToneVertex(double s) const {
  return T(s) / sqrt(intU);
}

cd MIsobarPiPiS::T(cd s) const {
  cd rho00 = sqrt(LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS)))/s;
  const cd scale = (s / (4 * POW2((K_MASS+K0_MASS)/2.))) - 1.;

  cd M00 = 0.0;
  for (uint i = 0; i < _sP.size(); i++) M00 += _a[i] / (s - _sP[i]);
  for (uint i = 0; i < _c.size(); i++)  M00 += pow(scale, static_cast<int>(i)) *_c[i];

  cd _T00 = 1./(M00 - cd(0, 1.)*rho00);
  return 16.*M_PI*_T00;
}

cd MIsobarPiPiS::T(double s) const {
  if (s < POW2(2*PI_MASS)) return 0.0;

  if (fabs(s - _sP[1]) < 1e-6) s  = POW2(sqrt(s)+1e-6);

  double rho00 = sqrt(LAMBDA(s, POW2(PI_MASS), POW2(PI_MASS)))/s;
  const double scale = (s / (4 * POW2((K_MASS+K0_MASS)/2.))) - 1.;

  double M00 = 0.0;
  for (uint i = 0; i < _sP.size(); i++) M00 += _a[i] / (s - _sP[i]);
  for (uint i = 0; i < _c.size(); i++)  M00 += pow(scale, static_cast<int>(i)) *_c[i];

  cd _T00 = 1./(M00 - cd(0, 1.)*rho00);
  return 16.*M_PI*_T00;
}
