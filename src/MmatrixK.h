// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MMATRIXK_H_
#define SRC_MMATRIXK_H_

#include <deflib.h>

#include <vector>
#include <iostream>

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>


#include <MModel.h>

namespace b = boost::numeric::ublas;

class MmatrixK : public MModel {
 public:
  // constructer
  MmatrixK(int Nch,
           int Npoles,
           std::vector<int> amap,
           b::matrix<cd> (*background)(double s, double *pars) =
           [](double, double *) -> b::matrix<cd> { return b::zero_matrix<cd>(3);} );

  void setPhSp(int i, double (&ph)(double s));
  void setProd(int i, cd (*prod)(double, const double *));

 private:
  int _Nch;
  int _Np;
  int _Na;
  std::vector<double (*)(double)> _fph;
  b::matrix<cd> (*_bkgr)(double, double*);
  // matrix-vector structures
  b::matrix<cd> _T;
  b::vector<cd> _alpha;
  b::vector<cd> _A;
  // parameters keepers
  std::vector<double> _mass;
  std::vector<double> _coupling;
  std::vector<double> _apar;
  // production machanism keeper
  std::vector<cd (*)(double, const double*)> _Aprod;
  std::vector<int> _amap;

 private:
  const b::vector<cd>& getA(double s, const double *pars);
  const b::vector<cd>& getA(double s);
  void calculateT(double s);
  void calculateA(double s);

 private:
  double last_s;
  bool need_to_recalculate_T;
  bool need_to_recalculate_A;

 public:
  void setPars(const double *pars);
  cd getA(int i, double s, const double *pars);
  cd getA(int i, double s);
  const b::matrix<cd>& getT(double s);
  // get constatns bach
  int getNch() const {return _Nch;}
  int getNp () const {return _Np ;}
  int getNa () const {return _Na ;}
  int getNpar () const {return _Np*(_Nch+1)+_Na;}

 public:
  void Print();
};

#endif  // SRC_MMATRIXK_H_
