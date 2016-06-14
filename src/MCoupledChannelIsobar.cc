// Copyright [03.09.2015] Mikhail Mikhasenko

#include "MCoupledChannelIsobar.h"

MCoupledChannelIsobar::MCoupledChannelIsobar(double Mi, double gi, double gpi,
                                             double m1i, double m2i, double m3i,
                                             double m1pi, double m2pi):
  MIsobar(Mi, 0.01, m1i, m2i, m3i), g(gi), gp(gpi), m1p(m1pi), m2p(m2pi) {;}

// channel one. Main channel

double MCoupledChannelIsobar::U(double s12) const {
  double Gtot = (g*g*RHO(s12, m1*m1, m2*m2)+
                 gp*gp*RHO(s12, m1p*m1p, m2p*m2p))/(2.*M);
  return g*g*RHO(s12, m1*m1, m2*m2)/(POW2(M*M-s12)+POW2(M*Gtot));
}

cd MCoupledChannelIsobar::U(cd s12) const {
  cd G  = g*g*RHO(s12, m1*m1, m2*m2)/(2.*M);
  cd Gp = gp*gp*RHO(s12, m1p*m1p, m2p*m2p)/(2.*M);
  cd Gtot = G+Gp;
  // return -2.*imag( G/Gtot * 1./(M*M-s12+cd(0.,M)*Gtot) );
  return 2.*M*G/(POW2(M*M-s12)+POW2(M*Gtot));
}

// channel two. Second channel

double MCoupledChannelIsobar::Up(double s12) const {
  double Gtot = (g*g *RHO(s12, m1*m1, m2*m2)+
                 gp*gp*RHO(s12, m1p*m1p, m2p*m2p))/(2.*M);
  return gp*gp*RHO(s12, m1p*m1p, m2p*m2p)/(POW2(M*M-s12)+POW2(M*Gtot));
}

cd MCoupledChannelIsobar::Up(cd s12) const {
  cd G  = g*g*RHO(s12, m1*m1, m2 *m2)/(2.*M);
  cd Gp = gp*gp*RHO(s12, m1p*m1p, m2p*m2p)/(2.*M);
  cd Gtot = G+Gp;
  return 2.*M*G/(POW2(M*M-s12)+POW2(M*Gtot));
}

/////////////////////////////

double MCoupledChannelIsobar::GetG0() const {
  return (g*g*RHO(M*M, m1*m1, m2*m2) +
          gp*gp*RHO(M*M, m1p*m1p, m2p*m2p))/(2.*M);
}
