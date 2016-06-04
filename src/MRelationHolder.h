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
  void AddRelation(const DP & intensity, std::function<double(double)> func);

 public:
  double CalculateChi2();

 private:
  // intensity
  std::vector<relation> store;

  void Print();
};

#endif  // SRC_MRELATIONHOLDER_H_
