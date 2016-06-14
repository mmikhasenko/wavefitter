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

class MIsobar {
 public:
  MIsobar(double M, double G0,
          double m1, double m2,
          int L = 0, double R = 5/*GeV*/);

 protected:
  double M;  // Breit-Wigner mass
  double G0;  // Breit-Wigner width
  double m1, m2;
  int L;
  double R;

 public:
  inline double GetM() const { return M; }
  inline double Mass() const { return M; }
  virtual double GetG0()  const { return G0; }
  virtual double Width() const { return G0; }
  inline void SetM0(double Mi) { M = Mi; }
  inline void SetG0(double G0i) { G0 = G0i; }

  inline double GetL() const { return L; }

  inline double sth() const { return POW2(m1+m2); }

 public:
  virtual double U(double s12) const;
  virtual cd     U(cd s)       const;

 private:
  std::vector<std::pair<double, double> > ltable;
};

#endif  // SRC_MISOBAR_H_
