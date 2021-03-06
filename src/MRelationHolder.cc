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
  status.push_back(true);
}

uint MRelationHolder::ExtractChi2Array(uint iR, double* arr) {
  if (iR >= store.size()) {
    std::cout << "Error<MRelationHolder::ExtractChi2Array> : _store.size()" << std::endl;
    return 0;
  }
  const relation & rel = store[iR];
  const DP & dps = rel.data;
  uint arr_index = 0;
  for (auto && dp : dps.data) {
    if (dp.x > dps.lrange && dp.x < dps.rrange) {
      if (arr == 0) { arr_index++; continue; }
      double yf = rel.func(dp.x);  // 0 is calculation option
      double diff = dp.y - yf;
      double sigma = dp.dy;
      arr[arr_index++] = POW2(diff/sigma);
    }
  }
  return arr_index;
}


double MRelationHolder::CalculateChi2() {
  double chi2 = 0;
  for (uint i=0; i < store.size(); i++) {
    if (!status[i]) continue;
    const relation & rel = store[i];
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

void MRelationHolder::Print() const {
  std::cout << "Relation content:" << "\n";
  for (uint i=0; i < store.size(); i++) {
    std::cout << "\t" << store[i].data.title << ": " << status[i] <<"\n";
  }
  std::cout << "\n";
}
