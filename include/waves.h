#ifndef __GKPY_H__
#define __GKPY_H__


#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#define POW2(v) ((v)*(v))
#define Ki(s,m) (sqrt((s)/4.-(m)*(m)))

#define WM(s,s0) ((sqrt(s)-sqrt(s0-s))/(sqrt(s)+sqrt(s0-s)))

#define PI_MASS 0.13957
#define K_MASS 0.496
#define ETA_MASS 0.54751

namespace waves {
  
  class GKPY {


  public:
    static double phi1(double s);
    static double phi2(double s);
    static double phi3(double s);
    static double phi4(double s);

    static double phi(double s);
  
  private:
    static const double B0;
    static const double B1;
    static const double B2;
    static const double B3;
    static const double z0;
    static const double d0;
    static const double c ;
    static const double B ;
    static const double C ;
    static const double D ;

  private:
    static const double M0;
    static const double M2;

  private:
    static const double phi1_M0;
    static const double phi1_dM0;
    static const double phi3_M2;
    static const double phi3_dM2;

  };

} // namespace


#endif

namespace waves {

  // true consts
  const double GKPY::B0 = 7.14;
  const double GKPY::B1 = -25.3;
  const double GKPY::B2 = -33.2;
  const double GKPY::B3 = -26.2;
  const double GKPY::z0 = PI_MASS;
  const double GKPY::d0 = 226.5;
  const double GKPY::c = -81;
  const double GKPY::B = 93.3;
  const double GKPY::C = 48.7;
  const double GKPY::D = -88.3;
  // where to add
  const double GKPY::M0 = 0.85;
  const double GKPY::M2 = 1.42;

  // calculatable
  const double GKPY::phi1_M0 = GKPY::phi1(POW2(M0));
  const double GKPY::phi1_dM0 = (GKPY::phi1(POW2(M0)+1.e-6)-GKPY::phi1(POW2(M0)-1.e-6))/(2.e-6);
  const double GKPY::phi3_M2 = GKPY::phi3(POW2(M2));
  const double GKPY::phi3_dM2 = (GKPY::phi3(POW2(M2)+1.e-6)-GKPY::phi3(POW2(M2)-1.e-6))/(2.e-6);

  // functions
  double GKPY::phi(double s) {
    if (s <= POW2(2*PI_MASS)) return 0;
    if (s <= POW2(M0)) return phi1(s);
    if (s <= POW2(2*K_MASS)) return phi2(s);
    if (s <= POW2(M2)) return phi3(s);
    return phi4(s);
  }

  double GKPY::phi1(double s) {
    double wm = WM(s, POW2(2*K_MASS));
    double cotD = sqrt(s)/(2.*Ki(s, PI_MASS))*POW2(PI_MASS)/(s-POW2(z0)/2.)*
      (
       POW2(z0)/(PI_MASS*sqrt(s))+B0+B1*wm+B2*wm*wm+B3*wm*wm*wm
       );
    double at = atan(1./cotD);
    return 180./M_PI*(at < 0 ? at + M_PI : at);
  }

  double GKPY::phi2(double s) {
    std::cout << phi1_M0 << ", " << phi1_dM0 << "\n";
    double aK2 = (s > POW2(2*K_MASS)) ? Ki(s, K_MASS) : sqrt(POW2(K_MASS)-s/4.);
    double K2M = sqrt(POW2(K_MASS)-POW2(M0)/4.);
    return d0*POW2(1.-aK2/K2M)+phi1_M0*aK2/K2M*(2.-aK2/K2M)+aK2*(K2M-aK2)*
      (8*phi1_dM0+c*(K2M-aK2)/(K_MASS*K_MASS*K_MASS));
  }

  double GKPY::phi3(double s) {
    double K2 = Ki(s, K_MASS);
    return d0+B*POW2(K2/K_MASS) + C*POW2(POW2(K2/K_MASS))+D*
      ((s > POW2(2*ETA_MASS)) ? POW2(Ki(s, ETA_MASS)/ETA_MASS) : 0);
  }

  double GKPY::phi4(double s) {
    std::cout << phi3_M2 << ", " << phi3_dM2 << "\n";
    double a = 360.;
    double alpha = phi3_dM2/(a-phi3_M2);
    double p = phi3_dM2/(a*alpha)*exp(alpha*POW2(M2));
    return a*(1.-p*exp(-alpha*s));
  }

} // namespace
