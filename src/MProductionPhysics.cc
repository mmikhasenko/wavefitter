// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <sstream>

#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mintegrate.h"
#include "deflib.h"

MProductionPhysics::MProductionPhysics(const std::vector<MChannel*> &channels) :
  MChannelPhysics<b::vector<cd> > (channels), _fB(0), _fC(0), _fCloc(false), _ubLookup(0) {
  std::cout << "------------------------------------\n";
  std::cout << "MProductionPhysics instance is created!" << std::endl;
}

/* add non-resonance production, it is functions to be called */
/* when it has been called then _getB is not empty */
void MProductionPhysics::addLongRange(const std::vector<std::function<cd(double)> > &getB,
                                      std::string par_name) {
  if (getB.size() != _Nch) std::cerr << "Error<void MProductionPhysics::addLongRange> : getB.size() != _Nch\n";

  _getB.resize(_Nch);
  for (uint i = 0; i < _Nch; i++) _getB[i] = getB[i];

  // fill long range term
  std::ostringstream ss; ss << par_name << "r";
  std::ostringstream iss; iss << par_name << "i";
  uint iRealPar = MParKeeper::gI()->add(ss.str(), 1.0, -10., 10.);
  uint iImagPar = MParKeeper::gI()->add(iss.str(), 0.0, -10., 10.);
  _fB.resize(_Nch);
  for (uint i = 0; i < _Nch; i++) {
    _fB[i].first  = iRealPar;
    _fB[i].second = iImagPar;
  }

  std::cout << "-----------> Long range production is added!\n";
}

void MProductionPhysics::addLongRangeSeparated(const std::vector<std::function<cd(double)> > &getB,
                                               std::string par_name) {
  if (getB.size() != _Nch) std::cerr << "Error<void MProductionPhysics::addLongRange> : getB.size() != _Nch\n";

  _getB.resize(_Nch);
  for (uint i = 0; i < _Nch; i++) _getB[i] = getB[i];

  _fB.resize(_Nch);
  for (uint i = 0; i < _Nch; i++) {
    // fill long range term
    std::ostringstream ss; ss << par_name << "r" << i;
    std::ostringstream iss; iss << par_name << "i" << i;
    _fB[i].first  = MParKeeper::gI()->add(ss.str(), 1.0, -10., 10.);
    _fB[i].second = MParKeeper::gI()->add(iss.str(), 0.0, -10., 10.);
  }

  std::cout << "-----------> Long range production is added!\n";
}

/* add directe resonance production, just set of parameters - production constants */
/* when it has been called then _fC is not empty */
void MProductionPhysics::addShortRange(std::string par_name) {
  _smap = [](double s)->cd{return s;};
  // add parameters to MParKeeper
  _fC.resize(1);
  _fC[0].resize(_Nch);
  for (uint i = 0; i < _Nch; i++) {
    std::ostringstream ss; ss << par_name << "r" << i;
    _fC[0][i].first  = MParKeeper::gI()->add(ss.str(), 0.0, -10., 10.);
    std::ostringstream iss; iss << par_name << "i" << i;
    _fC[0][i].second = MParKeeper::gI()->add(iss.str(), 0.0, -10., 10.);
  }
  std::cout << "-----------> Short range production is added!\n";
}

void MProductionPhysics::addShortRange(const std::vector<std::string> powers,
                                       std::function<cd(double)> smap) {
  _smap = smap;
  // copy powers array and add "r" and "i" part
  for (const auto & it : powers) {
    std::vector<std::pair<uint, uint> > vtmp(_Nch);
    for (uint i = 0; i < _Nch; i++) {
      std::ostringstream ss; ss << it << "r" << i;
      vtmp[i].first  = MParKeeper::gI()->add(ss.str(), 0.0, -10., 10.);
      std::ostringstream iss; iss << it << "i" << i;
      vtmp[i].second = MParKeeper::gI()->add(iss.str(), 0.0, -10., 10.);
    }
    _fC.push_back(vtmp);
  }
}

void MProductionPhysics::lockShortRangePhase() {
  _fCloc = true;
}

/* calculate unitarisation term and put it to a lookup table */
/* when it has been called then _ubLookup is not empty */
void MProductionPhysics::unitarize(double to, uint Npoints) {
  if (_getB.size() == 0 || _getT == 0) {
    std::cerr << "Error<MProductionPhysics::unitarize> You are trying to make unitarisation before T and B are defined\n";
    return;
  }
  // Fill lookup tables
  _ubLookup.resize(_Nch);
  for (uint i = 0; i < _Nch; i++) {
    _ubLookup[i].resize(Npoints);
    double sth = _iso[i]->sth();
    for (uint t = 0; t < Npoints; t++) {
      double s = sth + (to - sth)/(Npoints-1)*t;
      _ubLookup[i][t].first = s;
      cd value = 1./(2.*M_PI)*cintegrate([&](double up)->cd{
          double sp = 1./up;
          return _iso[i]->rho(sp)*_getB[i](sp)*_iso[i]->DumpC(sp)/(sp-s-cd(0., 1.e-4)) * 1./(up*up);
        }, 0, 1./sth);
      _ubLookup[i][t].second = value;
    }
  }
  std::cout << "-----------> Unitarisation is done. Lookup tables are built!\n";
}

/* calculate scattering matrix, final state interaction dynamics */
/* when it has been called then _getT = 0 */
void MProductionPhysics::addScattering(std::function<b::matrix<cd>(double)> getT) {
  _getT = getT;
}

void MProductionPhysics::calculate(double s) {
  // value I gonna calculate
  _value =  b::vector<cd>(_Nch, 0);

  // pure Deck
  b::vector<cd> fB(_Nch);
  if (_getB.size()) {
    b::vector<cd> B(_Nch);
    for (uint i = 0; i < _Nch; i++) {
      fB[i] = cd(MParKeeper::gI()->get(_fB[i].first), MParKeeper::gI()->get(_fB[i].second));  // get parameters
      B(i) = _getB[i](s) * fB[i];
    }
    // add to value
    _value += B;
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
    b::vector<cd> cvect(_Nch, 0);
    for (uint w = 0; w < _fC.size(); w++) {
      for (uint i = 0; i < _Nch; i++) {
        // work around !!!!!!
        double imag0 = MParKeeper::gI()->get(_fC[0][i].second);
        double real0 = MParKeeper::gI()->get(_fC[0][i].first);
        if (_fCloc && w != 0 && real0 != 0.0) {
          MParKeeper::gI()->set(_fC[w][i].second,
                                MParKeeper::gI()->get(_fC[w][i].first) *
                                imag0 / real0);
        }
        // work around !!!!!!
        cvect[i] += cd(MParKeeper::gI()->get(_fC[w][i].first),
                       MParKeeper::gI()->get(_fC[w][i].second)) *
          ((w == 0) ? 1 : ((w == 1) ?  _smap(s) : pow(_smap(s), w)));
      }
    }
    _value += prod(CThat, cvect);
  }

  // unitarisation
  if (_ubLookup.size()) {
    b::vector<cd> UB(_Nch);
    for (uint i = 0; i < _Nch; i++) UB(i) = getvalue(s, _ubLookup[i].data(), _ubLookup[i].size());
    b::vector<cd> CThatUB = prod(CThat, UB);
    // add to value
    for (uint j = 0; j < _Nch; j++) _value(j) += CThatUB(j) * fB(j);
  }
}
