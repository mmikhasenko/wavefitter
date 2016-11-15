// Copyright [14.07.2015] Misha Mikhasenko

#ifndef SRC_MDECK_H_
#define SRC_MDECK_H_

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <utility>

#include "TMath.h"
#include "Math/SpecFuncMathMore.h"

#include "deflib.h"
#include "mstructures.h"

#include "dFunction.hpp"

class MIsobar;

class MDeck {
 public:
  static cd getAmplitude(double costheta, double phi,
                         double mS1sq, double RS1,
                         double costheta_pr, double phi_pr,
                         double wsq, double t,
                         double mtRsq,
                         double stot,
                         double mAsq, double mBsq, double mDsq,
                         double m1sq, double m2sq, double m3sq);

  static cd getAmplitude(double p1x, double p1y, double p1z,  // pi-
                         double p2x, double p2y, double p2z,  // pi+
                         double t, double mtRsq, double stot,
                         double mAsq, double mBsq, double mDsq,
                         double m1sq, double m2sq, double m3sq);
};

#endif  // SRC_MDECK_H_
