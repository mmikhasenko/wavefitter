// Copyright [2016] Mikhail Mikhasenko

#include "mintegrate.hh"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

// integrate real and imaginary part
cd cintegrate(std::function<cd(double)> fint, double a, double b) {
  double real_part = integrate([&](double t) -> double {return real(fint(t));}, a, b);
  double imag_part = integrate([&](double t) -> double {return imag(fint(t));}, a, b);
  return cd(real_part, imag_part);
}

// basic integration function
double integrate(std::function<double(double)> fint, double a, double b) {
  // integrate
  std::function<double(double*, double*)> tfint = [&](double *x, double *p) {return fint(*x);};

  TF1 fdph("fdrho", tfint, a, b, 0);
  ROOT::Math::WrappedTF1 wph(fdph);  ROOT::Math::GaussIntegrator igph;
  igph.SetFunction(wph); igph.SetRelTolerance(1e-8);

  return igph.Integral(a, b);
}

