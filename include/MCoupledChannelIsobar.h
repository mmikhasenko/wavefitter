// Copyright [14.07.2015] Mikhail Mikhasenko

#ifndef SRC_MCOUPLEDCHANNELISOBAR_H_
#define SRC_MCOUPLEDCHANNELISOBAR_H_

#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>

#include "constants.h"
#include "deflib.h"
#include "MIsobar.h"

class MCoupledChannelIsobar : public MIsobar {
 public:
  MCoupledChannelIsobar(double M, double g1, double g2,
                        double m1, double m2, double m3,
                        double m1p, double m2p);

 protected:
  double g;
  double gp;

  double m1p;
  double m2p;

 public:
  virtual double GetG0()  const;
  virtual double Width() const { return GetG0(); }
  void SetG1(double g1) { g = g1; }
  void SetG2(double g2) { gp = g2; }

 public:
  virtual double U(double s12) const;
  virtual cd     U(cd s) const;
  double Up(double s12) const;
  cd Up(cd s) const;
};

#endif  // SRC_MCOUPLEDCHANNELISOBAR_H_
