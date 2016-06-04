// Copyright [2016] Mikhail Mikhasenko

#ifndef __SRC_MASSOCIATOION_HOLDER_H__
#define __SRC_MASSOCIATOION_HOLDER_H__

#include <deflib.h>
// #include <pair>
#include <iostream>
#include <vector>

#include "TGraphErrors.h"
#include "TH1D.h"

typedef data_points DP;

class MAssociationHolder {
 public:
  MAssociationHolder();

 public:
  void AddMap(const &DP intensity,
              double (*amp)(const double*, const double*, int), double Npar);

 private:
  // intensity
  std::vector< struct {
    const DP & data;
    double (*)(const double *x, const double* pars, int flag) func;
  } >;

  void Print();
};

#endif // __MASSOCIATOION_HOLDER_H__
