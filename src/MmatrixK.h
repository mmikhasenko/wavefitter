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

namespace b = boost::numeric::ublas;

class MmatrixK : public MChannelPhysics<b::matrix<cd> > {
 public:
  explicit MmatrixK(const std::vector<MIsobar*> &channels, uint Npoles = 0);
  MmatrixK(uint Nchannels, uint Npoles);

  void SetNpoles(uint Npoles);

 private:
  uint _Np;
  // matrix-vector structures
  b::matrix<cd> _T;
  // parameters keepers
  std::vector<uint> _mass;
  std::vector<uint> _coupling;

 private:
  void calculate(double s);

 public:
  uint getNp () const {return _Np ;}

 public:
  void Print();
};

#endif  // SRC_MMATRIXK_H_
