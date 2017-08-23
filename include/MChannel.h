// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MCHANNEL_H_
#define SRC_MCHANNEL_H_

#include <deflib.h>

#include <vector>
#include <utility>
#include "BlttWsskpf.h"

class MChannel {
 public:
  MChannel(uint iL, double iR) : L(iL), R(iR) {;}

 protected:
  uint L;
  double R;

 public:
  virtual double sth() const = 0;
  virtual double rho(double s) const = 0;
  virtual cd rho(cd s) const = 0;
  virtual double p(double s) const {return 4*M_PI*sqrt(s)*rho(s);}
  virtual cd p(cd s) const {return 4*M_PI*sqrt(s)*rho(s);}
  virtual double DumpC(double s) const { double q = p(s); return sqrt(FFMod::BlttWskpf[L](R*R*q*q)); }
  virtual     cd DumpC(cd s)     const {     cd q = p(s); return sqrt(FFMoc::BlttWskpf[L](R*R*q*q)); }
  virtual cd rholtilde(double s) const {
    if (dtable.size()) { return InterpolateRhoLtilda(s);
    } else { return DisperceRhoLtilda(s); }
  }
  virtual cd rholtilde(cd s) const { return DisperceRhoLtilda(s); }

  void makeDisperseLookupTable(double from, double to, uint Npoints);
  cd DisperceRhoLtilda(double s) const;
  cd DisperceRhoLtilda(cd s) const;
  cd InterpolateRhoLtilda(double s) const;

 public:
  uint GetL() const { return L; }
  uint GetR() const { return R; }
  void Print();

 private:
  std::vector<std::pair<double, cd> > dtable;
};

#endif  // SRC_MCHANNEL_H_
