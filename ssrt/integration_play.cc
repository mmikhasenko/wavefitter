#include <iostream>
#include <complex>
#include <iomanip>

#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
 
#include <constants.h>



double vint(double *x, double *par);

using namespace std; 
typedef std::complex<double> cd;

cd vint_sub(cd sp, cd s, double k, double s0, double sl, double sr);
double re_vint_sub(double *x, double *par);
double im_vint_sub(double *x, double *par);
cd integral(cd s, double k, double s0, double sl, double sr);
cd vf(cd s, double k, double s0, double sl, double sr);

int main(int ac, char **av) {

  if(ac!=2) {cout << "Usage: ./table k"<<endl; return 1;}
  int k = atoi(av[1]);
  //const char *fout = av[2];

  double s_th = pow(RHO_MASS+PI_MASS,2);
  double s_max = pow(2.5,2);
  const int Natt=1000; double step = (s_max-s_th)/Natt;
  double x[Natt],yr[Natt],yi[Natt],ya[Natt];

  /*------------------------------------------------------------------------------*/  
  
  for(int i=0;i<Natt;i++) {
    double s = s_th+step/2+i*step;
    cd val=vf(s,k+1,S1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2))+vf(s,k,S1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));
    x[i]=s;
    yr[i]=real(val); yi[i]=imag(val);
    ya[i]=norm(val);
  }

  TGraph gr(Natt,x,yr);
  TGraph gi(Natt,x,yi);
  TGraph ga(Natt,x,ya);
  TGraph gd(Natt,yr,yi);

  TCanvas can("c1"); can.Divide(2,2);
  can.cd(1); gr.Draw("apl");
  can.cd(2); gi.Draw("apl");
  can.cd(3); ga.Draw("apl"); 
  can.cd(4); gd.Draw("apl");
  can.SaveAs("/tmp/cint.png");

  for(int i=0;i<Natt;i++) 
    printf("%4.8f\t%3.8f\t%4.8f\n",x[i],yr[i],yi[i]);
    //cout<<setw(10)<<x[i]<<setw(10)<<yr[i]<<setw(10)<<yi[i]<<endl;

  return 0;
}

cd vint_sub(cd sp, cd s, double k, double s0, double sl, double sr) {
  cd rho_sp2 = (sp-sl)*(sp-sr)/(sp*sp);
  cd rho_s2  = (s -sl)*(s -sr)/(s*s);
  return (pow(w(sp,s0),k)*sqrt(rho_sp2)-pow(w(s,s0),k)*sqrt(rho_s2))/(sp*(sp-s));
}

double re_vint_sub(double *x, double *par) {
  cd spr(1/x[0],0.), sr(par[0],par[1]);
  return real(vint_sub(spr,sr,par[2],par[3],par[4],par[5]))/(x[0]*x[0]);
}
double im_vint_sub(double *x, double *par) {
  cd spr(1/x[0],0.), sr(par[0],par[1]);
  return imag(vint_sub(spr,sr,par[2],par[3],par[4],par[5]))/(x[0]*x[0]);
}

cd integral(cd s, double k, double s0, double sl, double sr) {

  TF1 fre("fre", re_vint_sub, 0, 1./sr,6);
  TF1 fim("fim", im_vint_sub, 0, 1./sr,6);
  double pars[] = {real(s),imag(s),k,s0,sl,sr};
  fre.SetParameters(pars);
  fim.SetParameters(pars);

  ROOT::Math::WrappedTF1 wre(fre);  ROOT::Math::GaussIntegrator igre;
  ROOT::Math::WrappedTF1 wim(fim);  ROOT::Math::GaussIntegrator igim;
  
  igre.SetFunction(wre); igre.SetRelTolerance(0.01);
  igim.SetFunction(wim); igim.SetRelTolerance(0.01);

  cd res(igre.Integral(0, 1./sr),igim.Integral(0, 1./sr));
  return res;
}

cd vf(cd s, double k, double s0, double sl, double sr) {
  cd rho_s2  = (s -sl)*(s -sr)/(s*s);
  cd mult = pow(w(s,s0),k)*sqrt(rho_s2)/s;
  cd unit(0,1);
  cd alog = 1.-(s+unit*EPSILON)/sr;
  cd first_int = integral(s,k,s0,sl,sr);
  return first_int-mult*log(alog);
}
