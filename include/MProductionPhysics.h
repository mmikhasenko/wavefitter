// Copyright [2016] Mikhail Mikhasenko

#ifndef __MPRODUCTIONPHYSICS_H__
#define __MPRODUCTIONPHYSICS_H__

#include <vector>
#include <iostream>
#include <utility>
#include <string>

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "MChannelPhysics.h"
#include "MChannel.h"

namespace b = boost::numeric::ublas;

class MProductionPhysics : public MChannelPhysics<b::vector<cd> >  {
 public:
  explicit MProductionPhysics(const std::vector<MChannel*> &channels);

 private:
  virtual void calculate(double s);
  virtual void calculate(cd s) {;}

 private:
  std::pair<uint, uint> _fB;
  std::vector<std::vector<std::pair<uint, uint> > > _fC;
  bool _fCloc;

 private:
  
  std::function<b::matrix<cd>(double)>    _getT;
  std::vector<std::function<cd(double)> > _getB;
  std::function<cd(double)> _smap;

 private:
  std::vector<std::vector<std::pair<double, cd> > > _ubLookup;

 public:
  void addLongRange(const std::vector<std::function<cd(double)> > &getB, std::string par_name = "B");
  void addShortRange(const std::vector<std::string> powers, std::function<cd(double)> getC);
  void addShortRange(std::string par_name = "c");
  void unitarize(double to = POW2(3.0), uint Npoints = 200);
  void addScattering(std::function<b::matrix<cd>(double)> getT);
  void lockShortRangePhase();
};

#endif  // __MPRODUCTIONPHYSICS_H__
