// Copyright [2016] Misha Mikhasenko

#ifndef __M3BODYANGULARBASIS_H__
#define __M3BODYANGULARBASIS_H__

#include <complex>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <functional>

namespace Math {
  std::complex<double> ZJMLS(uint J, int M, uint L, uint S,
                             double thetaI, double phiI,
                             double theta, double phi);
  std::complex<double> ZJMLS_refl(uint J, int M, bool pos_refl, uint L, uint S,
                                  double thetaI, double phiI,
                                  double theta, double phi);

  double integrate3bphs(std::function<double(double, double, double, double, double,
                                           double, double, double, double, double)> funct,
                      uint Npoints,
                      double s, double m1sq, double m2sq, double m3sq);
}
#endif
