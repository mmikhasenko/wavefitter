// Copyright [14.07.2015] Misha Mikhasenko

#ifndef SRC_MISOBARCHANNEL_H_
#define SRC_MISOBARCHANNEL_H_

#include <constants.h>
#include <deflib.h>

#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <utility>

#include "MChannel.h"
#include "MIsobar.h"

class MIsobarChannel : public MChannel {
 public:
  MIsobarChannel(MIsobar &iso,
                 double m3,
                 int L = 0, double R = 5);

 private:
  const MIsobar &_iso;

 protected:
  double _m3;

 public:
  const MIsobar &getIsobar() const {return _iso;}
  double getBachelorMass() const {return _m3;}

  double sth() const {return POW2(sqrt(_iso.sth())+_m3);}

  double CalculateQuasiTwoBody(double s) const;
  double InterpolateQuasiTwoBody(double s) const;

  cd CalculateQuasiTwoBodyStright(cd s) const;
  cd CalculateQuasiTwoBodyEdge(cd s) const;

  void makeLookupTable(double from, double to, uint Npoints);

  double rho(double s) const {
    if (s < sth()) return 0;
    if (ltable.size()) { return InterpolateQuasiTwoBody(s);
    } else { return CalculateQuasiTwoBody(s); }
  }
  // cd rho(cd s) const { return CalculateQuasiTwoBodyStright(s); }
  cd rho(cd s) const { return CalculateQuasiTwoBodyEdge(s); }

 private:
  std::vector<std::pair<double, double> > ltable;
};

#endif  // SRC_MISOBARCHANNEL_H_
