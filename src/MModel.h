// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MMODEL_H_
#define SRC_MMODEL_H_

#include <string>
#include <vector>
#include <map>

class MModel {
 public:
  MModel();

 protected:
  std::map<std::string, double*> map_to_pars;
  std::vector<double*> pool_pointer;

 public:
  void SetParameters(const double *pars);
  void PrintParameters();

  uint CreatePool(std::vector<std::string> names);
};

#endif  // SRC_MMODEL_H_
