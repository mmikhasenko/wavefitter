// Copyright [2016] Misha Mikhasenko

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "waves.h"

namespace waves {

  const double GKPY::pi_mass = 0.13957;
  const double GKPY::k_mass = 0.496;
  const double GKPY::eta_mass = 0.54751;
  // true consts
  const double GKPY::B0 = 7.14;
  const double GKPY::B1 = -25.3;
  const double GKPY::B2 = -33.2;
  const double GKPY::B3 = -26.2;
  const double GKPY::z0 = pi_mass;
  const double GKPY::d0 = 226.5;
  const double GKPY::c = -81;
  const double GKPY::B = 93.3;
  const double GKPY::C = 48.7;
  const double GKPY::D = -88.3;
  // where to add
  const double GKPY::M0 = 0.85;
  const double GKPY::M2 = 1.42;

  // calculatable
  const double GKPY::phi1_M0 = GKPY::phi1(uPOW2(M0));
  const double GKPY::phi1_dM0 = (GKPY::phi1(uPOW2(M0)+1.e-6)-GKPY::phi1(uPOW2(M0)-1.e-6))/(2.e-6);
  const double GKPY::phi3_M2 = GKPY::phi3(uPOW2(M2));
  const double GKPY::phi3_dM2 = (GKPY::phi3(uPOW2(M2)+1.e-6)-GKPY::phi3(uPOW2(M2)-1.e-6))/(2.e-6);

  // functions
  std::complex<double> GKPY::T(double s) {
    double delta = phi(s)/180.*M_PI;
    return 16*M_PI /sqrt(1.-uPOW2(2*pi_mass)/s) * sin(delta)*exp(std::complex<double>(0., delta));
  }

  double GKPY::phi(double s) {
    if (s <= uPOW2(2*pi_mass)) return 0;
    if (s <= uPOW2(M0)) return phi1(s);
    if (s <= uPOW2(2*k_mass)) return phi2(s);
    if (s <= uPOW2(M2)) return phi3(s);
    return phi4(s);
  }

  double GKPY::phi1(double s) {
    double wm = uWM(s, uPOW2(2*k_mass));
    double cotD = sqrt(s)/(2.*uKi(s, pi_mass))*uPOW2(pi_mass)/(s-uPOW2(z0)/2.)*
      (
       uPOW2(z0)/(pi_mass*sqrt(s))+B0+B1*wm+B2*wm*wm+B3*wm*wm*wm
       );
    double at = atan(1./cotD);
    return 180./M_PI*(at < 0 ? at + M_PI : at);
  }

  double GKPY::phi2(double s) {
    // std::cout << phi1_M0 << ", " << phi1_dM0 << "\n";
    double aK2 = (s > uPOW2(2*k_mass)) ? uKi(s, k_mass) : sqrt(uPOW2(k_mass)-s/4.);
    double K2M = sqrt(uPOW2(k_mass)-uPOW2(M0)/4.);
    return d0*uPOW2(1.-aK2/K2M)+phi1_M0*aK2/K2M*(2.-aK2/K2M)+aK2*(K2M-aK2)*
      (8*phi1_dM0+c*(K2M-aK2)/(k_mass*k_mass*k_mass));
  }

  double GKPY::phi3(double s) {
    double K2 = uKi(s, k_mass);
    return d0+B*uPOW2(K2/k_mass) + C*uPOW2(uPOW2(K2/k_mass))+D*
      ((s > uPOW2(2*eta_mass)) ? uPOW2(uKi(s, eta_mass)/eta_mass) : 0);
  }

  double GKPY::phi4(double s) {
    // std::cout << phi3_M2 << ", " << phi3_dM2 << "\n";
    double a = 360.;
    double alpha = phi3_dM2/(a-phi3_M2);
    double p = phi3_dM2/(a*alpha)*exp(alpha*uPOW2(M2));
    return a*(1.-p*exp(-alpha*s));
  }

}  // namespace waves
