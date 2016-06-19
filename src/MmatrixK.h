// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MMATRIXK_H_
#define SRC_MMATRIXK_H_

#include <deflib.h>

#include <vector>
#include <iostream>

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include <MChannelPhysics.h>
#include <MChannel.h>

namespace b = boost::numeric::ublas;

class MmatrixK : public MChannelPhysics<b::matrix<cd> > {
 public:
  explicit MmatrixK(const std::vector<MChannel*> &channels, uint Npoles = 0);
  MmatrixK(uint Nchannels, uint Npoles);

  void SetNpoles(uint Npoles);

  b::matrix<cd> getSSInverseValue(cd s);  // second sheet value

 private:
  uint _Np;
  // matrix-vector structures
  b::matrix<cd> _T;
  // parameters keepers
  std::vector<uint> _mass;
  std::vector<uint> _coupling;

 private:
  void calculate(double s) {tmpl_calculate<double>(s);}
  void calculate(cd s)     {tmpl_calculate<cd>(s);}

  template<typename sType>
    void tmpl_calculate(sType s);

 public:
  uint getNp () const {return _Np ;}

 public:
  void Print();
};

#endif  // SRC_MMATRIXK_H_
