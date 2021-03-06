// Copyright [15.07.2015] Mikhail Mikhasenko

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include "NoverD.h"
#include "deflib.h"

#define Sb 3.*3.
// #define w(a,b) (Sb-sqrt(a-b))/(Sb+sqrt(a-b))
#define w(a, b) (1.-sqrt(a-b))/(1.+sqrt(a-b))
#define EPSILON 1e-6

NoverD::NoverD(int k, int p, double s0, cd (*f)(cd), double sth, int NluPh, int NluVf, double hLim):
  _k(k), _p(p), _s0(s0), _f(f), _sth(sth), _hLim(hLim), _NluPh(NluPh), _NluVf(NluVf), lookup_ph(NluPh+1) {
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "---Building model: s0 = " << _s0 << ", k = "
            << _k << ", p = " << _p << ", sth = " << _sth << "---" << std::endl;
  std::cout << "-----------------" << _NluPh << "-x-" << _NluVf << "------------------" << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;
  // fill npar,ppar with zeros
  for (int i = 0; i < _k; i++) npar.push_back(0);
  for (int i = 0; i < _p; i++) ppar.push_back(0);
  // build lookup tables for phasespace
  for (int i=0; i < _NluPh; i++) {
    double u = 1./_sth/_NluPh*(i+1); double s = 1./u;
    lookup_ph[i] = std::make_pair(u, real(_f(cd(s, EPSILON))) );
  }
  double hDouble = 1e6;  cd hValue = _f(cd(hDouble, EPSILON));
  // std::cout << hValue << ",  " << 1./8/M_PI-real(hValue) << std::endl;
  lookup_ph[_NluPh] = std::make_pair(1./hDouble, real(hValue));
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "------------- Ph.sp. has been tabled --------------" << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;
  // build_lookup_tables for integrals
  for (int k = 0; k < _k; k++) {
    std::vector<std::pair<double, cd> > tvi;
    for (int i = 0; i < _NluVf; i++) {
      double s = _sth + (_hLim - _sth)/(_NluVf-1)*i;
      tvi.push_back(std::make_pair(s, vf(s, k)));
      // std::cout << "s = " << s << ", v = " << tvi[tvi.size()-1].second << std::endl;
    }
    lookup_vi.push_back(tvi);
  }
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "--------- Integrals have been tabled --------------" << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;
}

cd NoverD::expand(cd s, double s1, const std::vector<double> &par) {
  cd res = 0;
  for (uint i=0; i < par.size(); i++) res += pow(w(s, s1), i)*par[i];
  return res;
}
cd NoverD::cexpand(cd s, double s1, const std::vector<cd> &par) {
  cd res = 0;
  for (uint i = 0; i < par.size(); i++) res += pow(w(s, s1), i)*par[i];
  return res;
}

cd NoverD::D(double s) const {
  double padd = 0; for (uint i = 0; i < ppar.size(); i++) padd += ppar[i]*pow(s, i+1);
  cd dsum = 0;
  for (uint i = 0; i < npar.size(); i++) dsum += npar[i]*getvi(s, i);
  cd DfullI = 1. - s / (2*M_PI) * dsum   + padd;
  return  DfullI;
}
 
cd NoverD::DI(cd s) const {
  cd padd =0; for (uint i=0; i < ppar.size(); i++) padd += ppar[i]*pow(s, i+1);
  cd dsum = 0;
  for (uint i = 0; i < npar.size(); i++) dsum += npar[i]*vf(s, i);
  cd DfullI = 1. - s / (2*M_PI) * dsum   + padd;
  return  DfullI;
}

cd NoverD::Disc(cd s) const {
  cd rho = _f(s);
  cd DiscJust = cd(0., -1.) * rho * N(s);
  return DiscJust;
}

cd NoverD::DII(cd s) const {
  cd rho = _f(s);
  cd DII = DI(s) + cd(0, 1.) * rho * N(s);
  return DII;
}

// ********************************lookup tables*********************************** //

double NoverD::getrho(double s) const {
  double u = 1./s;
  double step = 1./_sth/_NluPh;
  int Nstep = u/step; if(Nstep == 0) {
    if(u<lookup_ph[_NluPh].first) return lookup_ph[_NluPh].second; 
    else return lookup_ph[0].second+
	   (lookup_ph[_NluPh].second-lookup_ph[0].second) / 
	   (lookup_ph[_NluPh].first -lookup_ph[0].first ) * (u-lookup_ph[0].first);
  }/* real(_f(cd(s,EPSILON))); std::cout<<"Warning: Number of points in phase space is too small. \nThe value will be calculated explicitely. "<< std::endl; */
  double value = lookup_ph[Nstep-1].second+
    (lookup_ph[Nstep].second-lookup_ph[Nstep-1].second) / 
    (lookup_ph[Nstep].first -lookup_ph[Nstep-1].first ) * (u-lookup_ph[Nstep-1].first);
  return value;
}

cd NoverD::getvi(double s, int k) const {
  //return vf(s,k);
  double step = (_hLim - _sth)/(_NluVf-1);
  int Nstep = (s-_sth)/step; if(Nstep>=(_NluVf-1)) {std::cout<<"Warning: Number of points in integrals is too small. \nThe value will be calculated explicitely. "<< std::endl; return vf(s,k);}
  cd value = lookup_vi[k][Nstep].second+
    (lookup_vi[k][Nstep+1].second-lookup_vi[k][Nstep].second) / 
    (lookup_vi[k][Nstep+1].first -lookup_vi[k][Nstep].first ) * (s-lookup_vi[k][Nstep].first);
  return value;
}

cd NoverD::getmax_value_of_integral(int k) {
  if(k>=_k) {std::cerr<<"Error<>NoverD::getmax_value_of_integral: k>=_k"<<std::endl; return 0;}
  cd max=0;
  for(int i=0;i<_NluVf;i++) if(abs(lookup_vi[k][i].second)>abs(max)) max = lookup_vi[k][i].second;
  return max;
}

//********************************************************************************//
//********************************************************************************//
//************************* Implementation of integrals **************************//
//********************************************************************************//
//********************************************************************************//

//********************************Real integrals**********************************//

cd NoverD::vf(double s, int k) const {
  double rho_s  = getrho(s);
  double mult = pow(w(s,_s0),k)*rho_s/s;
  cd alog = 1.-cd(s,EPSILON)/_sth;
  cd first_int = subIntegral(s,k);
  return first_int-mult*log(alog);
}

double NoverD::subIntegral(double s, int k) const {
  TF1 f("fsub",this,&NoverD::dSubI,0,1./_sth,2,"NoverD", "dSubI");
  double pars[] = {s,double(k)};
  f.SetParameters(pars);

  ROOT::Math::WrappedTF1 wra(f);  ROOT::Math::GaussIntegrator ig;
  ig.SetFunction(wra); ig.SetRelTolerance(1e-8);

  return ig.Integral(0, 1./_sth);
}

double NoverD::dSubI(double *x, double *par) const {
  double sp = 1/x[0], s = par[0];  int k(par[1]);
  double rho_sp = getrho(sp);
  double rho_s  = getrho(s );
  double subIntegrand = (pow(w(sp,_s0),k)*rho_sp-pow(w(s,_s0),k)*rho_s)/(sp*(sp-s));
  return subIntegrand/(x[0]*x[0]);
}

//*****************************Complex integrals**********************************//

cd NoverD::vf(cd s, int k) const {
  TF1 fre("fre",this,&NoverD::dvre,0,1./_sth,3,"NoverD", "dvre");
  TF1 fim("fim",this,&NoverD::dvim,0,1./_sth,3,"NoverD", "dvim");
  double pars[] = {real(s),imag(s),double(k)};
  fre.SetParameters(pars);
  fim.SetParameters(pars);

  ROOT::Math::WrappedTF1 wre(fre);  ROOT::Math::GaussIntegrator igre;
  ROOT::Math::WrappedTF1 wim(fim);  ROOT::Math::GaussIntegrator igim;  
  igre.SetFunction(wre); igre.SetRelTolerance(1e-8);
  igim.SetFunction(wim); igim.SetRelTolerance(1e-8);

  cd res(igre.Integral(0, 1./_sth),igim.Integral(0, 1./_sth));
  return res;
}

double NoverD::dvre(double *x, double *par) const {
  double sp = 1./x[0]; cd s(par[0],par[1]);
  int k(par[2]);
  double rho_sp = getrho(sp);
  cd integrand = pow(w(sp,_s0),k) * rho_sp / (sp*(sp-s));
  return real(integrand)/(x[0]*x[0]);
}

double NoverD::dvim(double *x, double *par) const {
  double sp = 1./x[0]; cd s(par[0],par[1]);
  int k(par[2]);
  double rho_sp = getrho(sp);  
  cd integrand = pow(w(sp,_s0),k) * rho_sp / (sp*(sp-s));
  return imag(integrand)/(x[0]*x[0]);
}
