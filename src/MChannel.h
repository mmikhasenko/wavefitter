// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MCHANNEL_H_
#define SRC_MCHANNEL_H_

#include <deflib.h>

#include <vector>
#include <utility>

class MChannel {
 public:
  MChannel(uint iL, double iR) : L(iL), R(iR) {;}

 protected:
  uint L;
  uint R;

 public:
  virtual double sth() const = 0;
  virtual double rho(double s) const = 0;
  virtual double p(double s) const {return 4*M_PI*sqrt(s)*rho(s);}
  virtual double DumpC(double s) const { double q = p(s); return pow(q*q/(1./(R*R)+q*q), L); }
  virtual cd rholtilde(double s) const {
    if (dtable.size()) { return InterpolateRhoLtilda(s);
    } else { return DisperceRhoLtilda(s); }
  }

  void makeDisperseLookupTable(double from, double to, uint Npoints);
  cd DisperceRhoLtilda(double s) const;
  cd InterpolateRhoLtilda(double s) const;

 public:
  void Print();

 private:
  std::vector<std::pair<double, cd> > dtable;
};

#endif  // SRC_MCHANNEL_H_
