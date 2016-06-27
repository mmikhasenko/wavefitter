// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MCHANNELPHYSICS_H_
#define SRC_MCHANNELPHYSICS_H_

#include <deflib.h>

#include <vector>
#include <iostream>

#include "MChannel.h"

template<class typeT>
class MChannelPhysics {
 public:
  explicit MChannelPhysics(uint Nchannels);
  explicit MChannelPhysics(const std::vector<MChannel*> &channels);

  void setPhSp(uint i, double (&ph)(double));

  typeT getValue(double s);
  typeT getValue(cd s);

 protected:
  uint _Nch;
  std::vector<MChannel*> _iso;

  virtual void calculate(double s) = 0;
  virtual void calculate(cd s) = 0;

  double last_s;
  bool need_for_recalculation;

  typeT _value;

 public:
  inline uint getNch() const {return _Nch;}
  void Print();
};

#define ERROR_VALUE -22.22

template<class typeT> MChannelPhysics<typeT>::MChannelPhysics(uint Nchannels) :
_Nch(Nchannels), _iso(Nchannels),
  last_s(ERROR_VALUE),
  need_for_recalculation(true) {
}

template<class typeT> MChannelPhysics<typeT>::MChannelPhysics(const std::vector<MChannel*> &channels) :
_Nch(channels.size()), _iso(channels.size()),
  need_for_recalculation(true),
  last_s(ERROR_VALUE) {
  // get information about channels
  std::cout << "-> Template MChannelPhysics constructor!\n";
  for (uint i = 0; i < _Nch; i++) _iso[i] = channels[i];
}

template<class typeT> void MChannelPhysics<typeT>::Print() {
  std::cout << "Number of channels: " << _Nch << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
}

template<class typeT> typeT MChannelPhysics<typeT>::getValue(double s) {
  if (s != last_s || need_for_recalculation) calculate(s);
  return _value;
}

template<class typeT> typeT MChannelPhysics<typeT>::getValue(cd s) {
  need_for_recalculation = true;
  calculate(s);
  return _value;
}

#endif  // SRC_MCHANNELPHYSICS_H_