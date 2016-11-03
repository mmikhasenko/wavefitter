// Copyright [2016] Misha Mikhasenko

#ifndef __GKPY_H__
#define __GKPY_H__

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>

#define uPOW2(v) ((v)*(v))
#define uKi(s, m) (sqrt((s)/4.-(m)*(m)))
#define uWM(s, s0) ((sqrt(s)-sqrt(s0-s))/(sqrt(s)+sqrt(s0-s)))

namespace waves {
  
  class GKPY {

   public:
    static double phi1(double s);
    static double phi2(double s);
    static double phi3(double s);
    static double phi4(double s);

    static double phi(double s);
    static std::complex<double> T(double s);
  
   private:
    static const double B0;
    static const double B1;
    static const double B2;
    static const double B3;
    static const double z0;
    static const double d0;
    static const double c;
    static const double B;
    static const double C;
    static const double D;

   private:
    static const double M0;
    static const double M2;

   private:
    static const double phi1_M0;
    static const double phi1_dM0;
    static const double phi3_M2;
    static const double phi3_dM2;

   private:
    static const double pi_mass;
    static const double k_mass;
    static const double eta_mass;
  };

};  // namespace waves

#endif  // __GKPY_H__
