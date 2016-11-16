// Copyright [14.07.2015] Misha Mikhasenko

#ifndef SRC_MDECK_H_
#define SRC_MDECK_H_

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <utility>

#include "deflib.h"

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

  static uint fromLabToGJ(double pAx, double pAy, double pAz, double mAsq,
                          double pBx, double pBy, double pBz, double mBsq,
                          double p1x, double p1y, double p1z, double m1sq,  // pi-
                          double p2x, double p2y, double p2z, double m2sq,  // pi+
                          double p3x, double p3y, double p3z, double m3sq,  // pi-
                          double mDsq = 0.);
};

#endif  // SRC_MDECK_H_
