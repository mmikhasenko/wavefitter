// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MMATRIXK_H_
#define SRC_MMATRIXK_H_

#include <deflib.h>

#include <vector>
#include <iostream>
#include <string>

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
  void addPole(const std::string &masssq_name, const std::string &par_name);
  void addBackground(const std::string &msq_name, const std::string &par_name);

  b::matrix<cd> getSSvalue(cd s);         // second sheet value
  b::matrix<cd> getSSInverseValue(cd s);  // second sheet value
  b::matrix<cd> getK(cd s);
  b::matrix<cd> getSSdenominator(cd s);
  b::matrix<cd> getFSdenominator(cd s);

 private:
  // matrix-vector structures
  b::matrix<cd> _T;
  // poles parameters keepers
  std::vector<uint> _masssq;
  std::vector<uint> _coupling;
  // background parameters keepers
  std::vector<uint> _bmasssq;
  std::vector<uint> _bcs;

 private:
  inline void calculate(double s) {tmpl_calculate<double>(s);}
  inline void calculate(cd s)     {tmpl_calculate<cd>(s);}

  template<typename sType>
    void tmpl_calculate(sType s);

 public:
  uint getNp () const {return _masssq.size();}

 public:
  void Print();
};

#endif  // SRC_MMATRIXK_H_
