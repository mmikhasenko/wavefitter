#ifndef __DEFLIB_H__
#define __DEFLIB_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <vector>

#include "TMultiGraph.h"
#include "TH1D.h"

#define LAMBDA(a,b,c) ((a)*(a)+(b)*(b)+(c)*(c)-2.*(a)*(b)-2.*(b)*(c)-2.*(c)*(a)) 
#define sqrtPi(x) sqrt(x*exp(-cd(0,1.)*M_PI))*exp(cd(0,1.)*M_PI/2.)
#define sqrtPhi(x,phi) sqrt(x*exp(-cd(0,1.)*phi))*exp(cd(0,1.)*phi)

typedef std::complex<double> cd;

TMultiGraph *split_data(const TH1D* h, double low, double up);

template <typename Type>
Type getvalue(double M, std::vector< std::pair<double,Type> > &table) {
  const int N = table.size();
  const double lft = table[0].first;
  const double rht = table[N-1].first;
  const double Mstep = table[1].first - lft;
  const int Nsteps = (M - lft)/Mstep;
  if(Nsteps<0 || Nsteps>=N-1) {std::cerr<<"Error!! in getvalue! M = "<<M<<", "<<table[0].second<<std::endl; return 0;}
  const Type value = table[Nsteps].second + 
    ( table[Nsteps+1].second - table[Nsteps].second ) / 
    ( table[Nsteps+1].first  - table[Nsteps].first  ) * (M - table[Nsteps].first); 
  return value;
}

double normPhase(double a);

#endif
