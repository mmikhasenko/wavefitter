// Copyright [14.07.2015] Misha Mikhasenko

#ifndef __ISOBARS_H__
#define __ISOBARS_H__

#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>

#include "./constants.h"

typedef std::complex<double> cd;

class MIsobar{
 public:
  MIsobar(double M, double G0,
          double m1, double m2, double m3,
          int L = 0, double R = 5/*GeV*/);

 protected:
  double M;  // Breit-Wigner mass
  double G0;  // Breit-Wigner width
  double m1, m2, m3;
  int L;
  double R;

 public:
  double GetM() const { return M; }
  double Mass() const { return M; }
  virtual double GetG0()  const { return G0; }
  virtual double Width() const { return G0; }
  void SetM0(double Mi) { M = Mi; }
  void SetG0(double G0i) { G0 = G0i; }

  double GetL() const { return L; }

 public:
  virtual double U(double s12) const;
  virtual cd     U(cd s)       const;
  double   rho(double s, double s12, double m3_2) const;
  cd       rho(cd     s, cd     s12, double m3_2) const;
  cd     rhoPi(cd     s, cd     s12, double m3_2) const;
  double  rho3(double s) const;
  cd      rho3(cd     s) const;

  // static double Btatt_Weisskopf(int Li, double Ri, double pi);
  // double BlttWsskpf(double pi);

 private:
  double  dph(double *x, double *par) const;
  double rdph(double *x, double *par) const;
  double idph(double *x, double *par) const;
};

//
//
// double MIsobar::BlttWsskpf(double pi, int Li, double Ri) {
//   if(L==0) return 1;
//   if(L==1) return 1./(1+pow(pi*Ri,2));
//   if(L==2) return 1./(3+pow(pi*Ri,2)+pow(pi*Ri,4));
//   return 1./(1+pow(pi*Ri,2*L));
// }
//
// double MIsobar::BlttWsskpf(double pi) {
//   return BlttWsskpf(pi,L,R);
// }
//
// double MIsobar::U(double s) {
//   double p  = lambda(s    ,m1*m1,m2*m2)/(2*sqrt(s));
//   double p0 = lambda(mI*mI,m1*m1,m2*m2)/(2*mI);
//   double bw = BlttWsskpf(p), bw0  = BlttWsskpf(p0);
//   double rho  = 2*p/sqrt(s), rho0 = 2*p0/mI;
//   double G = G0*rho*bw*pow(p,2)/(rho0*bw*pow(p0,2));
//   return 2*M*G/()
// }
//
// double MIsobar::Rho3() {
//  
// }





//
//double UDynamicWidth(double s, double mI, double G0, double m1, double m2, double L, double R) {
//  double p  = lambda(s    ,m1*m1,m2*m2)/(2*sqrt(s));
//  double p0 = lambda(mI*mI,m1*m1,m2*m2)/(2*mI);    
//  double bw  = 1./(1+pow(p *R,2*L)), bw0 = 1./(1+pow(p0*R,2*L));
//  double rho  = 2*p/sqrt(s), rho0 = 2*p0/mI;
//  double G = G0*rho*bw*pow(p,2*L)/(rho0*bw*pow(p0,2*L));
//  double U = 2*mI*G/(pow(mI*mI-s,2)+pow(mI*G,2));
//  return U;
//}
//
//double U_rho(double s) {
//  double U = UDynamicWidth(s,RHO_MASS,RHO_WIDTH,PI_MASS,PI_MASS,1,RHO_R);
//  return U;
//}
//
//double U_f2(double s) {
//  double U = UDynamicWidth(s,F2_MASS,F2_WIDTH,PI_MASS,PI_MASS,2,F2_R);
//  return U;
//}
//
//double d_ph(double *x, double *par) {
//  double s12 = x[0];
//  double s = par[0];
//  double mI=par[1], G0=par[2], R=par[3],
//    m1=par[4], m2=par[5], m3=par[6];
//
//  double p  = lambda(s,s12,m3*m3)/(2*sqrt(s));
//  double rho = 1./(8*M_PI)*2*p/sqrt(s);
//  double U = U_rho(s12,mI,G0,m1,m2,R);
//  return U*rho;
//}
//
//
//double phase_space_rho(double s, 
//		       double mI, double G0, double R,
//		       double m1, double m2, double m3) {
//  TF1 fre("fi", d_ph,pow(m1+m2,2),pow(sqrt(s)-m3,2),7);
//  double pars[] = {real(s),imag(s),k,s0,sl,sr};
//  fre.SetParameters(pars);
//
//  ROOT::Math::WrappedTF1 wre(fre);  ROOT::Math::GaussIntegrator igre;
//  
//  igre.SetFunction(wre); igre.SetRelTolerance(0.01);
//
//  return igre.Integral(0, 1./sr);
//}



#endif
