//Misha Mikhasenko
//14.07.2015

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include "MIsobar.h"
#include "deflib.h"
#include "mintegrate.hh"

#define SCALEX_FOR_CUT 2.
#define SCALEY_FOR_CUT 10.0

#define I_LOOKUP_NPOINT 150
#define I_HLIM 16.0

MIsobar::MIsobar(double Mi, double G0i,
                 double m1i, double m2i, double m3i,
                 int Li, double Ri,
                 bool Quasi):
  M(Mi), G0(G0i), m1(m1i), m2(m2i), m3(m3i), L(Li), R(Ri), quasi(Quasi) {;}

double MIsobar::rho(double s, double s12, double m3_2) {
  return (sqrt(s) > sqrt(s12)+sqrt(m3_2)) ? 1./(8*M_PI)*sqrt( LAMBDA(s, s12, m3_2)/(s*s) ) : 0;
}
cd     MIsobar::rho(cd     s, cd     s12, double m3_2) {
  return 1./(8*M_PI)*sqrt( LAMBDA(s, s12, m3_2)/(s*s) );
}
cd     MIsobar::rhoPi(cd     s, cd     s12, double m3_2) {
  return 1./(8*M_PI)*sqrtPi( LAMBDA(s, s12, m3_2)/(s*s) );
}

double MIsobar::U(double s) const {
  double G = G0 / rho(M*M, m1*m1, m2*m2)  *  rho(s, m1*m1, m2*m2);
  if (L > 0) {
    const double p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
    const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
    G *=  pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);
  }
  return 2.*M*G/(pow(M*M-s, 2)+pow(M*G, 2));
}
cd     MIsobar::U(cd     s) const {
  cd G = G0/rho(M*M, m1*m1, m2*m2)  *  rhoPi(s, m1*m1, m2*m2);
  if (L > 0) {
    const cd p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
    const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
    G *=  pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);
  }
  return -2.0*imag( 1./(M*M-s+cd(0, M)*G) );
  // return 2.*imag( 1./(M*M-s-cd(0, M)*G) );
}

//////////////////////////////////////////////////////////////////////////////
///  //   /  ///      //      ///     ///      ///      ///  ///////     /////
///  //      /////  ////  ///////  //////  //  ///  //  ///  ///////   ///////
///  //  /   /////  ////    /////  /  ///    /////      ///  /////////   /////
///  //  //  /////  ////      ///     ///  /   ///  //  ///      ///     /////
//////////////////////////////////////////////////////////////////////////////

void MIsobar::makeLookupTable() {
  ltable.resize(I_LOOKUP_NPOINT);
  for (int i = 0; i < I_LOOKUP_NPOINT; i++) {
    double s = POW2(m1+m2) + (I_HLIM-POW2(m1+m2))/(I_LOOKUP_NPOINT-1)*i;
    ltable[i].first = s;
    ltable[i].second = CalculateQuasiTwoBody(s);
  }
}

double MIsobar::CalculateQuasiTwoBody(double s) const {
  std::function<double(double)> drho = [&](double s12)->double {
    return 1./(2*M_PI)*U(s12)*rho(s, s12, POW2(m3));
  };
  return integrate(drho, POW2(m1+m2), POW2(sqrt(s)+m3));
}

cd MIsobar::CalculateQuasiTwoBodyStright(cd s) const {
  std::function<cd(double)> drho = [&](double t)->cd {
    cd s12 = POW2(m1+m2)+(POW2(sqrt(s)+m3)-POW2(m1+m2))*t;
    return U(s12)*rho(s, s12, m3);
  };
  return (POW2(sqrt(s)+m3)-POW2(m1+m2))*cintegrate(drho, 0, 1);
}

cd MIsobar::CalculateQuasiTwoBodyEdge(cd s) const {
  std::function<cd(double)> drho = [&](double t)->cd {
    double th1 = POW2(m1+m2);
    cd th2 = POW2(sqrt(s)-m3);
    cd thM(th1 + (real(th2)-th1)/SCALEX_FOR_CUT,
          imag(th2)/SCALEY_FOR_CUT);
    if (t < 1.) {
      cd s12 = th1+(thM-th1)*t;
      return 1./(2*M_PI)*U(s12)*rho(s, s12, m3*m3)  *  (thM-th1);
    } else {
      cd s12 = thM+(th2-thM)*(t-1);
      return 1./(2*M_PI)*U(s12)*rho(s, s12, m3*m3)  *  (th2-thM);
    }
  };
  return cintegrate(drho, 0, 2);
}
