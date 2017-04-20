//Misha Mikhasenko
//14.07.2015

#include "TF1.h"
#include "MIsobar.h"
#include "mintegrate.h"

MIsobar::MIsobar(double Mi, double G0i,
                 double m1i, double m2i,
                 uint Li, double Ri):
  M(Mi), G0(G0i), m1(m1i), m2(m2i), L(Li), R(Ri) {
  // calculate isobar shape integral
  intU = 1.0;
}

MIsobar * const MIsobar::setIntU() {
  if (intU != 1.0) return this;
  intU = integrate([&](double u)->double{
    return U(1./u)*1./(u*u)/(2*M_PI);
  }, 0, 1./sth());
  std::cout << "Integral over isobar shape is calculated. Value is " << intU << "\n";
  return this;
}

double MIsobar::IntU() const {
  if (intU == 1.0)
    std::cerr << "Error<> : integral over isobar shape has not been calculated yet, use MIsobar::setIntU()\n";
  return intU;
}


double MIsobar::U(double s) const {
  if (s < POW2(m1+m2)) return 0.0;
  double G = G0 / RHO(M*M, m1*m1, m2*m2)  *  RHO(s, m1*m1, m2*m2);
  const double p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
  const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
  G *=  pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);
  return 2.0*M*G/(POW2(M*M-s)+POW2(M*G)) / intU;
}

cd MIsobar::U(cd s) const {
  const cd p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
  const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
  cd BltWskp = pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);

  cd G = G0/RHO(M*M, m1*m1, m2*m2)  *  RHO(s, m1*m1, m2*m2) * BltWskp;
  return 2.0*M*G/(POW2(M*M-s)+POW2(M*G)) / intU;
}

cd MIsobar::ToneVertex(double s) const {
  if (s < POW2(m1+m2)) return 0.0;
  const double p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
  const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
  double BltWskp = pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);

  double G = G0/RHO(M*M, m1*m1, m2*m2)  *  RHO(s, m1*m1, m2*m2) * BltWskp;
  double gsq = 2*M*G0/RHO(M*M, m1*m1, m2*m2) * BltWskp;
  return sqrt(gsq/intU)/(M*M-s-cd(0., M)*G);  // / intU
}

cd MIsobar::T(cd     s) const {
  const cd p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
  const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
  cd BltWskp = pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);

  cd G = G0/RHO(M*M, m1*m1, m2*m2)  *  RHO(s, m1*m1, m2*m2) * BltWskp;
  cd gsq = 2*M*G0/RHO(M*M, m1*m1, m2*m2) * BltWskp;
  return gsq/(M*M-s-cd(0., M)*G);  // / intU
}

cd MIsobar::T(double s) const {
  const double p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
  const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
  double BltWskp = pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);

  double G = G0/RHO(M*M, m1*m1, m2*m2)  *  RHO(s, m1*m1, m2*m2) * BltWskp;
  double gsq = 2*M*G0/RHO(M*M, m1*m1, m2*m2) * BltWskp;
  return gsq/(M*M-s-cd(0., M)*G);  // / intU;
}
