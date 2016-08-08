// Copyright [08.08.2016] Misha Mikhasenko

#ifndef SRC_MASCOLI_H_
#define SRC_MASCOLI_H_

#include <iostream>
#define _USE_MATH_DEFINES

class MAscoli {
 public:
 MAscoli(double mAsq, double mBsq, double wsq, double mDsq, double mtRsq, double mIsq,
         double s, double t,
         uint S, int lamS, double R) :
  _mAsq(mAsq), _mBsq(mBsq), _wsq(wsq), _mDsq(mDsq), _mtRsq(mtRsq), _mIsq(mIsq),
    _s(s), _t(t), _S(S), _lamS(lamS), _R(R)/*, _table(0)*/ { ; }

 private:
  double _mAsq;
  double _mBsq;
  double _wsq;
  double _mDsq;
  double _mtRsq;
  double _mIsq;

  double _s;
  double _t;

 private:
  uint _J, _M, _L, _S;
  int _lamS;

  double _R;

// private:
//  std::vector<std::pair<double, double> > _ltable;

 public:
  static double getDeck(double mAsq, double mBsq, double wsq, double mDsq, double mtRsq,
                        double mIsq, double s, double t, double z,
                        uint S, int lamS, double R);

  static double getProjectedDeck(double mAsq, double mBsq, double wsq, double mDsq, double mtRsq,
                                 double mIsq, double s, double t,
                                 uint J, uint M, uint L, uint S, double R);

  inline double getValue(double wsq, double z) const {
    return getDeck(_mAsq, _mBsq, wsq, _mDsq, _mtRsq,
                   _mIsq, _s, _t, z,
                   _S, _lamS, _R);
  }

  inline double getProjection(double wsq, uint J, uint L) const  {
    return getProjectedDeck(_mAsq, _mBsq, wsq, _mDsq, _mtRsq,
                            _mIsq, _s, _t,
                            J, 0, L, _S, _R);
  }
  inline double getProjection(double wsq, double s1, uint J, uint L) const  {
    return getProjectedDeck(_mAsq, _mBsq, wsq, _mDsq, _mtRsq,
                            s1, _s, _t,
                            J, 0, L, _S, _R);
  }

};

#endif  // SRC_MASCOLI_H_
