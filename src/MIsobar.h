// Copyright [14.07.2015] Misha Mikhasenko

#ifndef SRC_MISOBAR_H_
#define SRC_MISOBAR_H_

#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <utility>

#include "constants.h"
#include "deflib.h"

class MIsobar{
 public:
  MIsobar(double M, double G0,
          double m1, double m2, double m3,
          int L = 0, double R = 5/*GeV*/,
          bool quasi = false);

 protected:
  double M;  // Breit-Wigner mass
  double G0;  // Breit-Wigner width
  double m1, m2, m3;
  int L;
  double R;
  bool quasi;

 public:
  double GetM() const { return M; }
  double Mass() const { return M; }
  virtual double GetG0()  const { return G0; }
  virtual double Width() const { return G0; }
  void SetM0(double Mi) { M = Mi; }
  void SetG0(double G0i) { G0 = G0i; }

  double GetL() const { return L; }

  inline double GetTwoBodySth() const { return POW2(M+m3); }
  inline double GetQuasiTwoBodySth() const { return POW2(m1+m2+m3); }
  inline double sth() const { return (quasi) ? GetQuasiTwoBodySth() : GetTwoBodySth(); }

  double GetTwoBodyBreakUpMom(double s) const { return 2*sqrt(LAMBDA(s, POW2(M), POW2(m3))/s); }
  double GetQuasiTwoBodyBreakUpMom(double s) const { return 4*M_PI*sqrt(s)*GetQuasiTwoBodyRho(s); }
  double p(double s) const { return (quasi) ? GetQuasiTwoBodyBreakUpMom(s) : GetTwoBodyBreakUpMom(s); }
  double DumpC(double s) const { double q = p(s); return pow(q*q/(1./(R*R)+q*q), L); }

  double CalculateQuasiTwoBody(double s) const;
  double InterpolateQuasiTwoBody(double s) const {
    if (!ltable.size()) return -1;
    return getvalue(s, ltable);
  }

  cd CalculateQuasiTwoBodyStright(cd s) const;
  cd CalculateQuasiTwoBodyEdge(cd s) const;

  void makeLookupTable();

 public:
  virtual double U(double s12) const;
  virtual cd     U(cd s)       const;

  // phase space
  static double rho(double s, double s12, double m3_2);
  static cd     rho(cd     s, cd     s12, double m3_2);
  static cd   rhoPi(cd     s, cd     s12, double m3_2);

  double   rho(double s) const { if (quasi) { return GetQuasiTwoBodyRho(s);} else {return GetTwoBodyRho(s); } }
  double   GetTwoBodyRho(double s) const {return (s > GetTwoBodySth()) ? rho(s, POW2(M), POW2(m3)) : 0.;}
  double   GetQuasiTwoBodyRho(double s) const {
    if (s < GetQuasiTwoBodySth()) return 0;
    if (ltable.size()) { return InterpolateQuasiTwoBody(s);
    } else { return CalculateQuasiTwoBody(s); }
  }

 private:
  std::vector<std::pair<double, double> > ltable;
  std::vector<std::pair<double, cd> > dtable;
};

#endif  // SRC_MISOBAR_H_
