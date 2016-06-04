// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_DEFLIB_H_
#define SRC_DEFLIB_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <vector>
#include <utility>

#define LAMBDA(a, b, c) ((a)*(a)+(b)*(b)+(c)*(c)-2.*(a)*(b)-2.*(b)*(c)-2.*(c)*(a))
#define sqrtPi(x) sqrt(x*exp(-cd(0, 1.)*M_PI))*exp(cd(0, 1.)*M_PI/2.)
#define sqrtPhi(x, phi) sqrt(x*exp(-cd(0, 1.)*phi))*exp(cd(0, 1.)*phi)
#define POW2(a) ((a)*(a))

typedef std::complex<double> cd;

#endif  // SRC_DEFLIB_H_
