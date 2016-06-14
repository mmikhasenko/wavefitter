//Misha Mikhasenko
//14.07.2015

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include "MIsobar.h"

MIsobar::MIsobar(double Mi, double G0i,
                 double m1i, double m2i,
                 int Li, double Ri):
  M(Mi), G0(G0i), m1(m1i), m2(m2i), L(Li), R(Ri) {;}

double MIsobar::U(double s) const {
  if (s < POW2(m1+m2)) return 0.0;
  double G = G0 / RHO(M*M, m1*m1, m2*m2)  *  RHO(s, m1*m1, m2*m2);
  const double p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
  const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
  G *=  pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);
  return 2.0*M*G/(POW2(M*M-s)+POW2(M*G));
}

cd     MIsobar::U(cd     s) const {
  cd G = G0/RHO(M*M, m1*m1, m2*m2)  *  RHO_PI(s, m1*m1, m2*m2);
  const cd p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
  const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
  G *=  pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);
  // return -2.0*imag( 1./(M*M-s+cd(0, M)*G) );
  // return 2.*imag( 1./(M*M-s-cd(0, M)*G) );
  return 2.0*M*G/(POW2(M*M-s)+POW2(M*G));
}
