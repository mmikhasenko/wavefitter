//Misha Mikhasenko
//03.09.2015

#include "MCoupledChannelIsobar.h"

MCoupledChannelIsobar::MCoupledChannelIsobar(double Mi, double gi, double gpi, 
			double m1i, double m2i, double m3i,
			double m1pi, double m2pi):
  MIsobar(Mi,0.01,m1i,m2i,m3i), g(gi), gp(gpi), m1p(m1pi), m2p(m2pi) {;}

//channel one. Main channel

double MCoupledChannelIsobar::U(double s12) const {
  double Gtot = (g *g *rho(s12,m1 *m1 ,m2 *m2 )+
		 gp*gp*rho(s12,m1p*m1p,m2p*m2p))/(2.*M);
  return g*g*rho(s12,m1*m1,m2*m2)/(pow(M*M-s12,2)+pow(M*Gtot,2));
}

cd MCoupledChannelIsobar::U(cd s12) const {
  cd G  = g *g *rho(s12,m1 *m1 ,m2 *m2 )/(2.*M);
  cd Gp = gp*gp*rho(s12,m1p*m1p,m2p*m2p)/(2.*M);
  cd Gtot = G+Gp;
  return -2.*imag( G/Gtot * 1./(M*M-s12+cd(0.,M)*Gtot) );
}

//channel two. Second channel

double MCoupledChannelIsobar::Up(double s12) const {
  double Gtot = (g *g *rho(s12,m1 *m1 ,m2 *m2 )+
		 gp*gp*rho(s12,m1p*m1p,m2p*m2p))/(2.*M);
  return gp*gp*rho(s12,m1p*m1p,m2p*m2p)/(pow(M*M-s12,2)+pow(M*Gtot,2));
}

cd MCoupledChannelIsobar::Up(cd s12) const {
  cd G  = g *g *rho(s12,m1 *m1 ,m2 *m2 )/(2.*M);
  cd Gp = gp*gp*rho(s12,m1p*m1p,m2p*m2p)/(2.*M);
  cd Gtot = G+Gp;
  return -2.*imag( G/Gtot * 1./(M*M-s12+cd(0.,M)*Gtot) );
}

/////////////////////////////

double MCoupledChannelIsobar::GetG0()  const { 
  return (g *g *rho(M*M,m1 *m1 ,m2 *m2 )+
	  gp*gp*rho(M*M,m1p*m1p,m2p*m2p))/(2.*M);
}
