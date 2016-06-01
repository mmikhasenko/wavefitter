#include <iostream>
#include <complex>

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

int main() {

  double s_th = pow(RHO_MASS+PI_MASS,2);
  double s_max = 3;
  const int Natt=100000; double step = (s_max-s_th)/Natt;
  double x[Natt],y[Natt];

  /*------------------------------------------------------------------------------*/
  /*------------------------------------------------------------------------------*/  
  
  for(int i=0;i<Natt;i++) {
    double s = s_th+i*step;
    cd val=integral(s,1,0.5,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));
    x[i]=s;
    y[i]=real(val);
  }
  
  /*------------------------------------------------------------------------------*/
  /*------------------------------------------------------------------------------*/
  /*
  TF1 fre("fre", re_vint_sub, 0, 1./s_th,5);
  TF1 fim("fim", im_vint_sub, 0, 1./s_th,5);
  double pars[] = {5,1,0.5,pow(RHO_MASS-PI_MASS,2),s_th};
  fre.SetParameters(pars);
  fim.SetParameters(pars);

  ROOT::Math::WrappedTF1 wre(fre);
  ROOT::Math::WrappedTF1 wim(fim);
  
  ROOT::Math::GaussIntegrator igre;
  ROOT::Math::GaussIntegrator igim;
  
  igre.SetFunction(wre); igre.SetRelTolerance(0.001);
  igim.SetFunction(wim); igim.SetRelTolerance(0.001);
  
  for(int i=0;i<Natt;i++) {
    double s = s_th+i*step;
    pars[0] = s;
    wre.SetParameters(pars);
    wim.SetParameters(pars);
    double re = igre.Integral(0, 1./s_th);
    double im = igim.Integral(0, 1./s_th);
    x[i]=s;
    y[i]=re;

  }
  */
  /*------------------------------------------------------------------------------*/


  TGraph gre(Natt,x,y);

  TCanvas can("c1");
  gre.Draw("apl");
  can.SaveAs("/tmp/c1.png");

  return 0;
}

cd vint_sub(cd sp, cd s, double k, double s0, double sl, double sr) {
  cd rho_sp2 = (sp-sl)*(sp-sr)/(sp*sp);
  cd rho_s2  = (s -sl)*(s -sr)/(s*s);
return (pow(w(sp,s0),k)*sqrt(rho_sp2)-pow(w(s,s0),k)*sqrt(rho_s2))/(sp*(sp-s));
}


double re_vint_sub(double *x, double *par) {
  cd spr(1/x[0],0.), sr(par[0],0.);
  return real(vint_sub(spr,sr,par[1],par[2],par[3],par[4]))/(x[0]*x[0]);
}
double im_vint_sub(double *x, double *par) {
  cd spr(1/x[0],0.), sr(par[0],0.);
  return imag(vint_sub(spr,sr,par[1],par[2],par[3],par[4]))/(x[0]*x[0]);
}

cd integral(cd s, double k, double s0, double sl, double sr) {

  TF1 fre("fre", re_vint_sub, 0, 1./sr,5);
  TF1 fim("fim", im_vint_sub, 0, 1./sr,5);
  double pars[] = {real(s),k,s0,sl,sr};
  fre.SetParameters(pars);
  fim.SetParameters(pars);

  ROOT::Math::WrappedTF1 wre(fre);
  ROOT::Math::WrappedTF1 wim(fim);
  
  ROOT::Math::GaussIntegrator igre;
  ROOT::Math::GaussIntegrator igim;
  
  igre.SetFunction(wre); igre.SetRelTolerance(0.001);
  igim.SetFunction(wim); igim.SetRelTolerance(0.001);

  cd res(igre.Integral(0, 1./sr),igim.Integral(0, 1./sr));
  return res;

}
