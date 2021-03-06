//Misha Mikhasenko
//15.07.2015

#ifndef __NOVERD_H__
#define __NOVERD_H__

#include <iostream>
#include <vector>
#include <complex>

#include <deflib.h>

typedef std::complex<double> cd;

class NoverD {

 public:
  NoverD(int k, int p, 
	 double s0, cd (*f)(cd), double sth, 
	 int NluPh=1000, int NluVf=1000, double hLim=9.0/*(3GeV)**2*/);
 private:
  int _k, _p;
  double _s0;
  cd (*_f)(cd s);
  double _sth, _hLim;
  // lookupn tables
  int _NluPh, _NluVf;
  std::vector<std::pair<double, double> > lookup_ph;
  std::vector<std::vector<std::pair<double, cd> > > lookup_vi;
  // std::pair<double,cd> **lookup_ph;

 public:
  std::vector<double> npar; 
  std::vector<double> ppar; 

 public:
  cd N(cd s) const {return expand(s,_s0,npar);}
  cd D(double s) const;
  cd DI(cd s) const;
  cd Disc(cd s) const;
  cd DII(cd s) const;
  cd A(double s) const {return N(s)/D(s);}
  cd AI(cd s) const {return N(s)/DI(s);}
  cd AII(cd s) const {return N(s)/DII(s);}
  cd ExtendA(double s, double as1, std::vector<cd> apar) {return cexpand(s,as1,apar)*A(s);}

  double GetThreshold() { return _sth; }

 public:
  cd getvi(double s, int k) const;
  double getrho(double s) const;
  cd getmax_value_of_integral(int k);

 public:
  static cd  expand(cd s, double s1, const std::vector<double> &par);
  static cd cexpand(cd s, double s1, const std::vector<cd> &par);
 
private:

  cd vf(double s, int k) const;
  double subIntegral(double s, int k) const;
  double dSubI(double *x, double *par) const;

  cd vf(cd s, int k) const;
  double dvre(double *x, double *par) const;
  double dvim(double *x, double *par) const;

 
};

#endif
