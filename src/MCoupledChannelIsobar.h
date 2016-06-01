//Misha Mikhasenko
//14.07.2015

#ifndef __COUPLEDCHANNELISOBARS_H__
#define __COUPLEDCHANNELISOBARS_H__

#include "MIsobar.h"

#include "constants.h"

#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>


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
  void SetG1(double g1) { g=g1; } 
  void SetG2(double g2) { gp=g2; } 
 
 public:

  virtual double U (double s12) const;
  virtual cd     U (cd s      ) const;
  double Up(double s12) const;
  cd     Up(cd s      ) const;
  

};


#endif
