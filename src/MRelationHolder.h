// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MRELATIONHOLDER_H_
#define SRC_MRELATIONHOLDER_H_

#include <deflib.h>
#include <iostream>
#include <vector>
#include <functional>

#include "mstructures.hh"

class MRelationHolder {
 public:
  static MRelationHolder *getInstance();
  static MRelationHolder *gI();
 private:
  MRelationHolder();

 public:
  // relation beetween data and function
  void AddRelation(const DP & intensity, std::function<double(double)> func);

 public:
  double CalculateChi2();

 private:
  static MRelationHolder *_ref;
  // intensity
  std::vector<relation> store;

  void Print();
};

#endif  // SRC_MRELATIONHOLDER_H_
