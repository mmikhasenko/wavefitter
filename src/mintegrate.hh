// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MINTEGRATE_H_
#define SRC_MINTEGRATE_H_

#include <functional>
#include <complex>

typedef std::complex<double> cd;
// integrate real and imaginary part
cd cintegrate(std::function<cd(double)> fint, double a, double b);

// basic integration function
double integrate(std::function<double(double)> fint, double a, double b);

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

// basic integration function
template<class Func >
double tintegrate(const Func &comp, double a, double b) {
  // integrate
  std::function<double(double*, double*)> tfint = [&](double *x, double *p) {return comp(*x);};

  TF1 fdph("fdrho", tfint, a, b, 0);
  ROOT::Math::WrappedTF1 wph(fdph);  ROOT::Math::GaussIntegrator igph;
  igph.SetFunction(wph); igph.SetRelTolerance(1e-8);

  return igph.Integral(a, b);
}

#endif  // SRC_MINTEGRATE_H_
