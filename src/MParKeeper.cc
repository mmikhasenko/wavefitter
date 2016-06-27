// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <cstdlib>
#include <algorithm>

#include "MParKeeper.h"
	
#define ERROR_VALUE -11.11

MParKeeper *MParKeeper::_ref = 0;

MParKeeper::MParKeeper() {
  std::cout << "--------MParKeeper.Constructor--------" << "\n";
}

MParKeeper *MParKeeper::getInstance() {
  if (!_ref) _ref = new MParKeeper();
  return _ref;
}
MParKeeper *MParKeeper::gI() {
  if (!_ref) _ref = new MParKeeper();
  return _ref;
}

uint MParKeeper::add(const std::string &pname, double v0, double lrange, double rrange) {
  _pars.push_back(v0);
  _map.push_back({std::string(pname), uint(_pars.size()-1), lrange, rrange});
  return _map.size()-1;
}

// name <-> index relations
uint MParKeeper::getIndex(const std::string &name) const {
  auto it = std::find_if(_map.begin(), _map.end(),
            [&](const parameter_description &d)->bool{
                           return (d.name == name) ? true : false;
                         });
  if(it == _map.end()) {std::cerr << "Error<MParKeeper::getIndex>: no name \"" << name << "\" found\n" ;return 0;}
  return it->index;
}

void MParKeeper::setRange(uint i, double v1, double v2) {
  check(i, _pars.size(), "setRange(i)"); 
  _map[i].lrange = v1;
  _map[i].rrange = v2;
}

/* operations with main set */
void MParKeeper::setRandom(uint i) {
  check(i, _pars.size(), "setRandom");
  double r01 = 1.*std::rand()/RAND_MAX;
  _pars[i] = _map[i].lrange+(_map[i].rrange-_map[i].lrange)*r01;
}

void MParKeeper::pset(const double *pars) {
  for (uint i=0; i < _pool.size(); i++) _pars[_pool[i]] = pars[i];
}

uint MParKeeper::pgetIndex(const std::string &name) const {
  auto it = std::find_if(_pool.begin(), _pool.end(),
                         [&, name](uint i)->bool{
                           return (_map[i].name == name) ? true : false;
                         });
  if(it == _pool.end()) {std::cerr << "Error<MParKeeper::pgetIndex>: no name \"" << name << "\" found in the pool\n" ;return 0;}
  return it - _pool.begin();
}


/*********************** P O O L ************************/

void MParKeeper::makePool(const std::vector<uint> &indexes) {
  const uint size = indexes.size();
  _pool.resize(size);
  for (uint i=0; i < size; i++) _pool[i] = indexes[i];
}
void MParKeeper::makePool(const std::vector<std::string> &names) {
  std::vector<uint> indexes(names.size());
  std::transform(names.begin(), names.end(), indexes.begin(),
                 [&](const std::string &s)->uint{return getIndex(s);});
  makePool(indexes);
}

void MParKeeper::check(uint i, uint N, const char *message) const {
  if (i >= N) std::cerr << "Error<MParKeeper::check in " << message << ">: i>=N\n";
}
void MParKeeper::printAll() {
  for (auto && it : _map) std::cout << it.index << ") "
                                    << it.name << ": "
                                    << _pars[it.index] << " in (" << it.lrange
                                    << ", " << it.rrange << ")\n";
  std::cout << "----> Pool: ";
  for (uint & it : _pool) std::cout << it << " "; 
  std::cout << "\n";
}
