// Copyright [1988] CERNLIB

#ifndef __WignerD_H__
#define __WignerD_H__

#include <complex>

double WignerD(double aj, double am, double an, double beta);
std::complex<double> WignerD(double aj, double am, double an, double alpha, double beta, double gamma);

#endif
