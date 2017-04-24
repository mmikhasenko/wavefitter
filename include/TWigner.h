// Copyright [1988] CERNLIB

#ifndef __WignerD_H__
#define __WignerD_H__

#include <complex>

namespace Math {
  double WignerD(int aj, int am, int an, double beta);
  std::complex<double> WignerD(int aj, int am, int an, double alpha, double beta, double gamma);
  std::complex<double> WignerD_refl(int aj, int am, int an, bool pos_refl,
                                        double alpha, double beta, double gamma);
}
#endif
