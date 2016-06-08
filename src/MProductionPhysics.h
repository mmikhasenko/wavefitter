// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MPRODUCTIONPHYSICS_H_
#define SRC_MPRODUCTIONPHYSICS_H_

#include <vector>
#include <iostream>
#include <utility>

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "MChannelPhysics.h"
#include "MIsobar.h"

#define UB_LOOKUP_NPOINT 150

namespace b = boost::numeric::ublas;

class MProductionPhysics : public MChannelPhysics<b::vector<cd> >  {
 public:
  explicit MProductionPhysics(const std::vector<MIsobar*> &channels);

 private:
  void calculate(double s);

 private:
  std::vector<int> _fC;
  std::vector<int> _fB;

 private:
  std::function<b::matrix<cd>(double)> _getT;
  std::vector<std::function<cd(double)> > _getB;

 private:
  std::vector<std::array<std::pair<double, cd>, UB_LOOKUP_NPOINT> > _ubLookup;

 public:
  void addLongRange(const std::vector<std::function<cd(double)> > &getB);
  void addShortRange();
  void unitarize();
  void addScattering(std::function<b::matrix<cd>(double)> getT);
};

#endif  // SRC_MPRODUCTIONPHYSICS_H_
