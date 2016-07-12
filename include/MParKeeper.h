// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MPARKEEPER_H_
#define SRC_MPARKEEPER_H_

#include <string>
#include <vector>
#include <map>
#include <utility>

class MParKeeper {
 private:
  MParKeeper();

 public:
  static MParKeeper *getInstance();
  static MParKeeper *gI();

  uint add(const std::string &pname, double v0 = 0., double lrange = -1., double rrange = 1.);
  inline uint nPars() const {return _pars.size();}
  /* main */
  double get(uint i) const {check(i, _pars.size()); return _pars[i];}
  void set(uint i, double v) {check(i, _pars.size()); _pars[i] = v;}
  void setRandom(uint i);
  void setRange(uint i, double v1, double v2);
  // name <-> index relations
  const std::string & getName(uint i) const {check(i, _pars.size()); return _map[i].name;}
  uint getIndex(const std::string &name) const;
  // methods called by name
  const std::vector<double> & get() const {return _pars;}
  inline double get(const std::string &pname) const {return get(getIndex(pname));}
  inline void   set(const std::string &pname, double v) {set(getIndex(pname), v);}
  inline void setRandom(const std::string &pname) {setRandom(getIndex(pname));}
  inline void setRange(const std::string &pname, double v1, double v2)  {setRange(getIndex(pname), v1, v2);}

  /* pool */
  void makePool(const std::vector<std::string> &names);
  void makePool(const std::vector<uint> &indexes);
  inline uint pnPars() const {return _pool.size();}

  const std::vector<uint> & getPool() const {return _pool;}
  void randomizePool() { for (uint & it : _pool) setRandom(it);}
  void pset(const double *pars);
  /* operations with pool */
  inline double pget(uint i) const {check(i, _pool.size()); return get(_pool[i]);}
  inline void pset(uint i, double v) {check(i, _pool.size()); return set(_pool[i], v);}
  inline void psetRandom(uint i) {check(i, _pool.size()); return setRandom(_pool[i]);}
  inline void psetRange(uint i, double v1, double v2) {check(i, _pool.size()); return setRange(_pool[i], v1, v2);}
  // name <-> index relations
  const std::string &pgetName(uint i) const {check(i, _pool.size()); return getName(_pool[i]);}
  uint pgetIndex(const std::string &name) const;

 private:
  static MParKeeper *_ref;

  std::vector<double> _pars;
  typedef struct {
    std::string name;
    uint index;
    double lrange;
    double rrange;
  } parameter_description;
  std::vector<parameter_description> _map;

  std::vector<uint> _pool;

 public:
  void check(uint i, uint N, const char *message = "DEFAULT") const;
  void printAll();
};

#endif  // SRC_MPARKEEPER_H_
