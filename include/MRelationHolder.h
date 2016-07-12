// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MRELATIONHOLDER_H_
#define SRC_MRELATIONHOLDER_H_

#include <deflib.h>
#include <iostream>
#include <vector>
#include <functional>

#include "mstructures.h"

class MRelationHolder {
 public:
  static MRelationHolder *getInstance();
  static MRelationHolder *gI();

 private:
  MRelationHolder();

 public:
  uint Nrels() const { return store.size(); }

  // relation beetween data and function
  void AddRelation(const DP & intensity, std::function<double(double)> func);
  const relation & GetRelation(uint iR) {
    if (iR >= store.size()) {
      std::cout << "Error<GetRelation> iR>=_store.size()" << std::endl;
      return store[0];
    }
    return store[iR];
  }
  bool relationStatus(uint iR) {
    if (iR >= store.size()) {
      std::cout << "Error<status(uint iR)>=_store.size()" << std::endl;
      return false;
    }
    return status[iR];
  }

 public:
  double CalculateChi2();
  void passiveAll() {for (uint i=0; i < status.size(); i++) status[i] = false;}
  void activateRelation(uint i) {
    if (i >= store.size()) {std::cerr << "Error<void activateRelation>\n"; return;}
    status[i] = true;
  }

 private:
  static MRelationHolder *_ref;
  // intensity
  std::vector<relation> store;
  std::vector<bool> status;

 public:
  void Print() const;
};

#endif  // SRC_MRELATIONHOLDER_H_
