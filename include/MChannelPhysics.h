// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MCHANNELPHYSICS_H_
#define SRC_MCHANNELPHYSICS_H_

#include <deflib.h>

#include <vector>
#include <iostream>
#include <map>
#include <utility>

#include "MChannel.h"

#define MCHANNELPHYSICS_LARGE_SIZE 250

template<class typeT>
class MChannelPhysics {
 public:
  explicit MChannelPhysics(uint Nchannels);
  explicit MChannelPhysics(const std::vector<MChannel*> &channels);

  void setPhSp(uint i, double (&ph)(double));

  typeT getValue(double s);
  typeT getValue(cd s);

  void RecalculateNextTime() {_smap.clear();}

 protected:
  uint _Nch;
  std::vector<MChannel*> _iso;

  virtual void calculate(double s) = 0;
  virtual void calculate(cd s) = 0;

  std::map<double, typeT> _smap;

  typeT _value;

 public:
  inline uint getNch() const {return _Nch;}
  void Print();
};

#define ERROR_VALUE -22.22

template<class typeT> MChannelPhysics<typeT>::MChannelPhysics(uint Nchannels) :
_Nch(Nchannels), _iso(Nchannels) {
}

template<class typeT> MChannelPhysics<typeT>::MChannelPhysics(const std::vector<MChannel*> &channels) :
_Nch(channels.size()), _iso(channels.size()) {
  // get information about channels
  std::cout << "-> Template MChannelPhysics constructor!\n";
  for (uint i = 0; i < _Nch; i++) _iso[i] = channels[i];
}

template<class typeT> void MChannelPhysics<typeT>::Print() {
  std::cout << "Number of channels: " << _Nch << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
}

template<class typeT> typeT MChannelPhysics<typeT>::getValue(double s) {
  if (_smap.find(s) == _smap.end()) { calculate(s); _smap[s] = _value;
  } else { _value = _smap[s]; }
  // Awere user if he forgot to recalculate
  // if (_smap.size() > MCHANNELPHYSICS_LARGE_SIZE)
  //   std::cout << "Warning<MChannelPhysics<typeT>::getValue>: size of s-store pool is too large. \n"
  //             << "Please check that you are refreshing amplitude calculator in time or increase MCHANNELPHYSICS_LARGE_SIZE.\n";
  // This check is consumering, should be removed later
  //  if (_smap[s] != _value) {
  //    std::cout << "Warning<MChannelPhysics<typeT>::getValue>: check your code! "
  //              << _smap[s] << "!=" << _value << "\n";
  //  }
  return _value;
}

template<class typeT> typeT MChannelPhysics<typeT>::getValue(cd s) {
  calculate(s);
  return _value;
}

#endif  // SRC_MCHANNELPHYSICS_H_
