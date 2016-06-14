// Copyright [2016] Mikhail Mikhasenko

#include "MmatrixK.h"
#include "MParKeeper.h"
#include "MatrixInverse.h"

#include <deflib.h>

MmatrixK::MmatrixK(uint Nch,
                   uint Npoles) :
  // mother class
  MChannelPhysics<b::matrix<cd> >(Nch),
  _Np(0) {
  // allocate poles parameters
  SetNpoles(Npoles);
  // init matrix
  _value = b::matrix<cd>(_Nch, _Nch);
  std::cout << "--------------------------------\n";
  std::cout << "MmatrixK instance is created!" << std::endl;
}

MmatrixK::MmatrixK(const std::vector<MIsobarChannel*> &channels, uint Npoles) :
  // mother class
  MChannelPhysics<b::matrix<cd> >(channels),
  _Np(0) {
  // allocate poles parameters
  if (Npoles) SetNpoles(Npoles);
  // init matrix
  _value = b::matrix<cd>(_Nch, _Nch);
  std::cout << "MmatrixK instance is created!" << std::endl;
  std::cout << "  channel parameters are set!" << std::endl;
}

void MmatrixK::SetNpoles(uint Npoles) {
  if (_Np != 0) {std::cerr << "Error<MmatrixK::SetNpoles>: try to change nomber of poles!" << std::endl; return;}
  _Np = Npoles;
  _mass.resize(_Np);
  _coupling.resize(_Nch*Npoles);

  // Fill parameters map
  for (uint ipole = 0; ipole < Npoles; ipole++) {
    // generate name
    std::ostringstream mname; mname << "m" << ipole + 1;
    _mass[ipole] = MParKeeper::gI()->add(mname.str(), 1.702);
    for (uint jch = 0; jch < _Nch; jch++) {
      std::ostringstream gname; gname << char(103+ipole + ((ipole >= 7) ? 1 : 0))
                                      << jch+1;
      _coupling[ipole*_Nch+jch] = MParKeeper::gI()->add(gname.str(), 0.0);
    }
  }
  std::cout << "parameters for " << _Np << " poles are allocated!" << std::endl;
}

void MmatrixK::calculate(double s) {
  // clear K
  b::symmetric_matrix<cd, b::upper> K = b::zero_matrix<cd>(_Nch);  // (_Nch) K*=0.;
  for (uint i = 0; i < _Np; i++) {
    b::symmetric_matrix<cd, b::upper> km(_Nch);
    double gpart[_Nch];
    for (uint j = 0; j < _Nch; j++) gpart[j] = MParKeeper::gI()->get(_coupling[i*_Nch+j]);
    double mass = MParKeeper::gI()->get(_mass[i]);
    for (uint j = 0; j < _Nch; j++)
      for (uint t = j; t < _Nch; t++)
        km(j, t) = gpart[j]*gpart[t]/(mass*mass-s);
    K += km;
  }
  // ph.sp.matrix
  b::symmetric_matrix<cd, b::upper> mrho(_Nch);
  // possibly add some background
  for (uint i = 0; i < _Nch; i++)
    for (uint j = 0; j < _Nch; j++)
      mrho(i, j) = (i == j) ? _iso[i]->rho(s) : 0.0;

  b::matrix<cd> irhoK = cd(0, 1)*prod(mrho, K);
  b::matrix<cd> din = b::identity_matrix<cd>(_Nch) - irhoK;

  bool sing = false;
  b::matrix<cd> din_inv = gjinverse(din, sing);
  // std::cout << din_inv << std::endl;
  if (sing) {std::cerr << "ERROR: SINGULAR" << std::endl;}  // exit(); }

  // finally
  _value = prod(K, din_inv);
  // change the flag
  last_s = s; need_for_recalculation = false;
}

void MmatrixK::Print() {
  std::cout << "------------ " << _Nch << "-channel K-matrix --------------" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Number of poles: " << _Np << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
}
