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

class MHelicityDeck {
 public:
  MHelicityDeck(double m1sq, double m2sq, double m3sq, double m4sq, double mtsq,
        uint J, uint L,
        uint Sp, int lamP,
        uint S , int lamS,
        double R) :
  _m1sq(m1sq), _m2sq(m2sq), _m3sq(m3sq), _m4sq(m4sq), _mtsq(mtsq),
    _J(J), _L(L), _Sp(Sp), _lamP(lamP), _S(S), _lamS(lamS), _R(R), _ltable(0) { ; }

 private:
  double _m1sq;
  double _m2sq;
  double _m3sq;
  double _m4sq;
  double _mtsq;

 private:
  uint _J, _L;
  uint _Sp; int _lamP;
  uint _S;  int _lamS;
  double _R;

 private:
  std::vector<std::pair<double, double> > _ltable;

 public:
  static double getDeck(double m1sq, double m2sq, double m3sq, double m4sq, double mtsq,
                        double s, double z,
                        uint Sp, int lamP,
                        uint S , int lamS, double R);

  static double getProjectedDeck(double m1sq, double m2sq, double m3sq, double m4sq, double mtsq,
                                 double s,
                                 uint J, uint L,
                                 uint Sp, int lamP,
                                 uint S, double R);

 public:
  inline cd getValue(double s, double z) const { return getDeck(_m1sq, _m2sq, _m3sq, _m4sq, _mtsq, s, z, _Sp, _lamP, _S, _lamS, _R); }
  inline cd getValue(double s, double z, double s3) const { return getDeck(_m1sq, _m2sq, _m3sq, _m4sq, _mtsq, s, z, _Sp, _lamP, _S, _lamS, _R); }
  inline cd getValue(double s, double z, double s3, double t) const { return getDeck(_m1sq, t, s3, _m4sq, _mtsq, s, z, _Sp, _lamP, _S, _lamS, _R); }

 public:
  inline double getProjection(double s) const  { return getProjectedDeck(_m1sq, _m2sq, _m3sq, _m4sq, _mtsq, s, _J, _L, _Sp, _lamP, _S, _R); }
  inline double getProjection(double s, double s3) const { return getProjectedDeck(_m1sq, _m2sq, s3, _m4sq, _mtsq, s, _J, _L, _Sp, _lamP, _S, _R); }
  inline double getProjection(double s, double s3, double t) const {
    return getProjectedDeck(_m1sq, t, s3, _m4sq, _mtsq, s, _J, _L, _Sp, _lamP, _S, _R); }

 public:
  void setLookupTable(const std::vector<std::pair<double, double> > &ltable);
  void makeLookupTable(const MIsobar &iso, double m3, double from, double to, uint Npoints);
  double getPrecalculated(double s) const;
  const std::vector<std::pair<double, double> > & getLookupTable() const { return _ltable; }

 public:
  inline double sth() const { return POW2(sqrt(_m3sq)+sqrt(_m4sq)); }
  inline uint J() const { return _J; }
  inline uint L() const { return _L; }
  inline uint Sp() const { return _Sp; }
  inline uint lamP() const { return _lamP; }
  inline uint S() const { return _S; }
  inline uint lamS() const { return _lamS; }

 public:
  void Print() const;
};

#endif  // SRC_MDECK_H_
