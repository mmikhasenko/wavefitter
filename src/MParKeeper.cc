// Copyright [2016] Mikhail Mikhasenko

#include "MParKeeper.h"
#include <iostream>

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

double MParKeeper::get(uint i) const {
  if (i > _pars.size()) { std::cerr << "Error<MParKeeper::get>: i>=size\n"; return ERROR_VALUE; }
  return _pars[i];
}

double MParKeeper::get(std::string pname) const {
  auto it = _map.find(pname);
  if (it == _map.end()) { std::cerr << "Error<MParKeeper::get(string)>: not found in the map\n"; return ERROR_VALUE; }
  return _pars[it->second];
}

void MParKeeper::set(uint i, double v) {
  if (i > _pars.size()) { std::cerr << "Error<MParKeeper::set>: i>=size\n"; return; }
  _pars[i] = v;
}

void MParKeeper::set(std::string pname, double v) {
  auto it = _map.find(pname);
  if (it == _map.end()) { std::cerr << "Error<MParKeeper::set(string)>: " << pname << " not found in the map\n"; }
  _pars[it->second] = v;
}

uint MParKeeper::add(std::string pname, double v0) {
  _pars.push_back(v0);
  uint index = _pars.size()-1;
  _map[pname] = index;
  if (_map.size() != _pars.size()) { std::cerr << "Error<MParKeeper::add>: _map, _pars mismatch\n"; return 0; }
  return index;
}

void MParKeeper::setPool(const double *pars) {
  for (uint i=0; i < pool.size(); i++) pool[i] = pars[i];
}

void MParKeeper::makePool(std::vector<std::string> names) {
  const uint size = names.size();
  pool.resize(size);
  for (uint i=0; i < size; i++) {
    auto it = _map.find(names[i]);
    if (it == _map.end()) { std::cerr << "Error<MParKeeper::makePool(string)>: not found in the map\n"; return; }
    pool[i] = it->second;
  }
}
void MParKeeper::makePool(std::vector<uint> indexes) {
  const uint size = indexes.size();
  pool.resize(size);
  for (uint i=0; i < size; i++) pool[i] = indexes[i];
}

void MParKeeper::printAll() {
  for (auto && it : _map) std::cout << it.first << ": " << _pars[it.second] << "\n";
}
