// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MPARKEEPER_H_
#define SRC_MPARKEEPER_H_

#include <string>
#include <vector>
#include <map>

class MParKeeper {
 private:
  MParKeeper();

 public:
  static MParKeeper *getInstance();
  static MParKeeper *gI();
  double get(int i) const;
  double get(std::string pname) const;
  void set(int i, double v);
  int add(std::string pname, double v0 = 0);

 private:
  static MParKeeper *_ref;
  std::map<std::string, int> _map;
  std::vector<double> _pars;
};

#endif  // SRC_MPARKEEPER_H_
