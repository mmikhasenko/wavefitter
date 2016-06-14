// Copyright [2016] Misha Mikhasenko

#ifndef SRC_MTWOBODYCHANNEL_H_
#define SRC_MTWOBODYCHANNEL_H_

#include <deflib.h>

#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <utility>

#include "MChannel.h"

class MTwoBodyChannel : public MChannel {
 public:
  MTwoBodyChannel(double im1, double im2,
                 int L = 0, double R = 5) : MChannel(L, R), m1(im1), m2(im2) {;}

 protected:
  double m1;
  double m2;

 public:
  virtual double sth() const {return POW2(m1+m2);}

  virtual double rho(double s) const {
    if (s < sth()) { return 0;
    } else { return RHO(s, m1*m1, m2*m2); }
  }
};

#endif  // SRC_MTWOBODYCHANNEL_H_
