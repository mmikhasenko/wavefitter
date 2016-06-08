// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <sstream>

#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mintegrate.hh"
#include "deflib.h"

#define UB_HLIM 4.0

MProductionPhysics::MProductionPhysics(const std::vector<MIsobar*> &channels) :
  MChannelPhysics<b::vector<cd> > (channels),
  _fB(0), _fC(0) {
  std::cout << "------------------------------------\n";
  std::cout << "MProductionPhysics instance is created!" << std::endl;
}

/* add non-resonance production, it is functions to be called */
/* when it has been called then _getB is not empty */
void MProductionPhysics::addLongRange(const std::vector<std::function<cd(double)> > &getB) {
  if (getB.size() != _Nch) std::cerr << "Error<void MProductionPhysics::addLongRange>\n";

  _getB.resize(_Nch);
  for (uint i = 0; i < _Nch; i++) _getB[i] = getB[i];
  _fB.resize(_Nch);
  // fill long range term
  for (uint i = 0; i < _Nch; i++) {
    std::ostringstream ss; ss << "b" << i;
    _fB[i] = MParKeeper::gI()->add(ss.str(), 0.0);
  }
}

/* add directe resonance production, just set of parameters - production constants */
/* when it has been called then _fC is not empty */
void MProductionPhysics::addShortRange() {
  // fill short range term
  _fC.resize(_Nch);
  for (uint i = 0; i < _Nch; i++) {
    std::ostringstream ss; ss << "c" << i;
    _fC[i] = MParKeeper::gI()->add(ss.str(), 0.0);
  }
}

/* calculate unitarisation term and put it to a lookup table */
/* when it has been called then _ubLookup is not empty */
void MProductionPhysics::unitarize() {
  if (_getB.size() == 0 || _getT == 0) {
    std::cerr << "Error<MProductionPhysics::unitarize> You are trying to make unitarisation before T and B are defined\n";
    return;
  }
  // Fill lookup tables
  _ubLookup.resize(_Nch);
  for (uint i = 0; i < _Nch; i++) {
    double sth = _iso[i]->sth();
    for (uint t = 0; t < UB_LOOKUP_NPOINT; t++) {
      double s = sth + (UB_HLIM - sth)/(UB_LOOKUP_NPOINT-1)*t;
      _ubLookup[i][t].first = s;
      _ubLookup[i][t].second = 1./(2.*M_PI)*cintegrate([&](double uprime)->cd{
          return _iso[i]->rho(1./uprime)*_getB[i](1./uprime)*_iso[i]->DumpC(1./uprime)/(1./uprime-s);},
        0, 1./sth);
    }
  }
}

/* calculate scattering matrix, final state interaction dynamics */
/* when it has been called then _getT = 0 */
void MProductionPhysics::addScattering(std::function<b::matrix<cd>(double)> getT) {
  _getT = getT;
}

void MProductionPhysics::calculate(double s) {
  // value I gonna calculate
  _value =  b::vector<cd>(_Nch, 0);

  b::vector<cd> bvect;

  // pure Deck
  if (_getB.size()) {
    b::vector<cd> B(_Nch);
    for (uint i = 0; i < _Nch; i++)
      B(i) = _getB[i](s);  // getvalue(s, [i].data(), table[i].size());
    // get parameters
    for (uint i = 0; i < _Nch; i++) bvect[i] = MParKeeper::gI()->get(_fB[i]);
    // add to value
    _value += element_prod(B, bvect);
  }
  if (!_fC.size() && !_ubLookup.size()) return;  // if no unitarisation, no direct production -> return

  // scattering
  b::matrix<cd> CThat(_Nch, _Nch);
  const b::matrix<cd> &T = _getT(s);
  b::vector<cd> DumpC(_Nch);
  for (uint i = 0; i < _Nch; i++) DumpC(i) = 1.;  // _iso[i]->DumpC(s);

  for (uint i = 0; i < _Nch; i++)
    for (uint j = 0; j < _Nch; j++) CThat(i, j) = T(i, j)*DumpC(j);
  // direct production
  if (_fC.size()) {
    b::vector<cd> cvect(_Nch);
    for (uint i = 0; i < _Nch; i++) cvect[i] = MParKeeper::gI()->get(_fC[i]);
    // add to value
    _value += prod(CThat, cvect);
  }

  // unitarisation
  if (_ubLookup.size()) {
    b::vector<cd> UB(_Nch);
    for (uint i = 0; i < _Nch; i++) UB(i) = getvalue(s, _ubLookup[i].data(), _ubLookup[i].size());
    b::vector<cd> CThatUB = prod(CThat, UB);
    // add to value
    _value += element_prod(CThatUB, bvect);
  }
}
