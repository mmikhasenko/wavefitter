// Copyright [2016] Mikhail Mikhasenko

#include <iostream>

#include "MModel.h"

MModel::MModel() {
}

void MModel::SetParameters(const double* pars) {
  for (int i = 0; i < pool_pointer.size(); i++) *pool_pointer[i] = pars[i];
}

uint MModel::CreatePool(std::vector<std::string> names) {
  pool_pointer.clear();
  for (const std::string & it : names) {
    auto m = map_to_pars.find(it);
    if (m != map_to_pars.end()) pool_pointer.push_back(m->second);
  }
  return pool_pointer.size();
}

void MModel::PrintParameters() {
  for (const double* m : pool_pointer) std::cout << *m << "\n";
}

// if ('*' is found) {
//   for (auto && m : map_to_pars)
//     if (std::regex_match(m.first, std::regex(it))) pool_pointer.push_back(m.second);
// }
