// Copyright [2016] Mikhail Mikhasenko

#include "MRelationHolder.h"

MRelationHolder *MRelationHolder::_ref = 0;

MRelationHolder::MRelationHolder() {
  std::cout << "--------MRelationHolder.Constructor--------" << "\n";
}

MRelationHolder *MRelationHolder::getInstance() {
  if (!_ref) _ref = new MRelationHolder();
  return _ref;
}
MRelationHolder *MRelationHolder::gI() {
  if (!_ref) _ref = new MRelationHolder();
  return _ref;
}

void MRelationHolder::AddRelation(const DP &intensity,
                                  std::function<double(double)> func) {
  relation r {intensity, func};
  store.push_back(r);
}

double MRelationHolder::CalculateChi2() {
  double chi2 = 0;
  for (const relation & rel : store) {
    const DP & dps = rel.data;
    for (auto && dp : dps.data) {
      if (dp.x > dps.lrange && dp.x < dps.rrange) {
        double yf = rel.func(dp.x);  // 0 is calculation option
        double diff = dp.y - yf;
        double sigma = dp.dy;
        if (dp.dy < 0) std::cerr << "Something is completely wrong!" << std::endl;
        chi2 += POW2(diff/sigma);
      }
    }
  }
  return chi2;
}

void MRelationHolder::Print() {
  std::cout << "Relation content:" << "\n";
  for (auto && rel : store) {
    std::cout << "\t" << rel.data.title << "\n";
  }
  std::cout << "\n";
}
