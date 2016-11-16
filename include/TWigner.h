// Copyright [1988] CERNLIB

#ifndef __WignerD_H__
#define __WignerD_H__

#include <complex>

double WignerD(int aj, int am, int an, double beta);
std::complex<double> WignerD(int aj, int am, int an, double alpha, double beta, double gamma);

#endif
