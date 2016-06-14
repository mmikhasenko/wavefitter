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
  double get(uint i) const;
  double get(std::string pname) const;
  void set(uint i, double v);
  void set(std::string pname, double v);
  uint add(std::string pname, double v0 = 0);

  void setPool(const double *pars);
  void makePool(std::vector<std::string> names);
  void makePool(std::vector<uint> names);

  inline uint nPars() const {return _pars.size();}

  void printAll();

 private:
  static MParKeeper *_ref;
  std::map<std::string, uint> _map;
  std::vector<double> _pars;

  std::vector<uint> pool;
};

#endif  // SRC_MPARKEEPER_H_
