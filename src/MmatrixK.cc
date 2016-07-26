// Copyright [2016] Mikhail Mikhasenko

#include "MmatrixK.h"
#include "MParKeeper.h"
#include "MatrixInverse.h"

#include <deflib.h>

MmatrixK::MmatrixK(uint Nch,
                   uint Npoles) :
  // mother class
  MChannelPhysics<b::matrix<cd> >(Nch) {
  // allocate poles parameters
  SetNpoles(Npoles);
  // init matrix
  _value = b::matrix<cd>(_Nch, _Nch);
  std::cout << "--------------------------------\n";
  std::cout << "MmatrixK(" << _Nch << "x" << _Nch << ") instance is created!" << std::endl;
}

MmatrixK::MmatrixK(const std::vector<MChannel*> &channels, uint Npoles) :
  // mother class
  MChannelPhysics<b::matrix<cd> >(channels),
  _bmass(0), _bcs(0) {
  // allocate poles parameters
  if (Npoles) SetNpoles(Npoles);
  // init matrix
  _value = b::matrix<cd>(_Nch, _Nch);
  std::cout << "MmatrixK(" << _Nch << "x" << _Nch << ") instance is created!" << std::endl;
  std::cout << "  channel parameters are set!" << std::endl;
}

void MmatrixK::SetNpoles(uint Npoles) {
  if (_mass.size() != 0) {std::cerr << "Error<MmatrixK::SetNpoles>: try to change nomber of poles!" << std::endl; return;}
  _mass.resize(Npoles);
  _coupling.resize(_Nch*Npoles);

  // Fill parameters map
  for (uint ipole = 0; ipole < Npoles; ipole++) {
    // generate name
    std::ostringstream mname; mname << "m" << ipole;
    _mass[ipole] = MParKeeper::gI()->add(mname.str(), 1.702, 1.001, 3.001);
    for (uint jch = 0; jch < _Nch; jch++) {
      std::ostringstream gname; gname << char(103+ipole + ((ipole >= 7) ? 1 : 0))
                                      << jch;
      _coupling[ipole*_Nch+jch] = MParKeeper::gI()->add(gname.str(), 0.0, -10., 10.);
    }
  }
  std::cout << "parameters for " << Npoles << " poles are allocated!" << std::endl;
}

void MmatrixK::addPole(const std::string &mass_name, const std::string &par_name) {
  // add index to MParKeeper for mass and couplings
  _mass.push_back(MParKeeper::gI()->add(mass_name, 1.702, 1.001, 3.001));
  for (uint jch = 0; jch < _Nch; jch++) {
    std::ostringstream gname; gname << par_name << jch;
    _coupling.push_back(MParKeeper::gI()->add(gname.str(), 0.0, -10., 10.));
  }
  std::cout << "parameters " << mass_name << " and " << par_name << "0.." << _Nch
            << " are allocated!" << std::endl;
}

void MmatrixK::addBackground(const std::string &bmass_name, const std::string &par_name) {
  _bmass.push_back(MParKeeper::gI()->add(bmass_name, 1, 0, 3));
  for (uint jch = 0; jch < _Nch; jch++) {
    std::ostringstream gname; gname << par_name << jch;
    _bcs.push_back(MParKeeper::gI()->add(gname.str(), 0.0, -10., 10.));
  }
  std::cout << "parameters " << bmass_name << " and " << par_name << "0.." << _Nch
            << " are allocated!" << std::endl;
}


template<typename sType>
void MmatrixK::tmpl_calculate(sType s) {
  // clear K
  b::symmetric_matrix<cd, b::upper> K = getK(s);
  // ph.sp.matrix
  b::symmetric_matrix<cd, b::upper> mrho(_Nch);
  // possibly add some background
  for (uint i = 0; i < _Nch; i++)
    for (uint j = 0; j < _Nch; j++) {
      mrho(i, j) = (i == j) ? _iso[i]->rholtilde(s) : 0.0;
    }

  b::matrix<cd> irhoK = cd(0., 0.5)*prod(mrho, K);
  b::matrix<cd> din = b::identity_matrix<cd>(_Nch) - irhoK;

  bool sing = false;
  b::matrix<cd> din_inv = gjinverse(din, sing);
  // std::cout << din_inv << std::endl;
  if (sing) {std::cerr << "ERROR: SINGULAR" << std::endl;}  // exit(); }

  // finally
  _value = prod(K, din_inv);
}

b::matrix<cd> MmatrixK::getK(cd s) {
  // clear K
  b::symmetric_matrix<cd, b::upper> K = b::zero_matrix<cd>(_Nch);  // (_Nch) K*=0.;
  // add poles terms
  const uint Npoles = _mass.size();
  for (uint i = 0; i < Npoles; i++) {
    double gpart[_Nch];
    for (uint j = 0; j < _Nch; j++) gpart[j] = MParKeeper::gI()->get(_coupling[i*_Nch+j]);
    double mass = MParKeeper::gI()->get(_mass[i]);
    for (uint j = 0; j < _Nch; j++)
      for (uint t = j; t < _Nch; t++)
        K(j, t) += gpart[j]*gpart[t]/(mass*mass-s);
  }
  // add background terms
  for (uint i = 0; i < _bmass.size(); i++) {
    double gpart[_Nch];
    for (uint j = 0; j < _Nch; j++) gpart[j] = MParKeeper::gI()->get(_bcs[i*_Nch+j]);
    double bmass = MParKeeper::gI()->get(_bmass[i]);
    for (uint j = 0; j < _Nch; j++)
      for (uint t = j; t < _Nch; t++)
        K(j, t) += gpart[j]*gpart[t]/(POW2(bmass)+s);
  }
  return K;
}

b::matrix<cd> MmatrixK::getSSInverseValue(cd s) {
  std::cerr << "Error<MmatrixK::getSSInverseValue>:The function is not supported anylonger, use \"getSSdemoninator(s)\"!\n";
  return b::identity_matrix<cd>(_Nch);
}

b::matrix<cd> MmatrixK::getSSdenominator(cd s) {
  // value on the first sheet
  b::matrix<cd> DI = getFSdenominator(s);
  // ph.sp.matrix
  b::symmetric_matrix<cd, b::upper> mrho(_Nch);
  for (uint i = 0; i < _Nch; i++)
    for (uint j = 0; j < _Nch; j++)
      mrho(i, j) = (i == j) ? _iso[i]->rho(s)*_iso[i]->DumpC(s) : 0.0;

  b::matrix<cd> K = getK(s);
  b::matrix<cd> i2rhoK = cd(0, 1)*prod(mrho, K);
  b::matrix<cd> DII = DI + i2rhoK;
  return DII;
}

b::matrix<cd> MmatrixK::getFSdenominator(cd s) {
  b::matrix<cd> K = getK(s);
  // ph.sp.matrix
  b::symmetric_matrix<cd, b::upper> mrho(_Nch);
  // possibly add some background
  for (uint i = 0; i < _Nch; i++)
    for (uint j = 0; j < _Nch; j++) {
      mrho(i, j) = (i == j) ? _iso[i]->rholtilde(s) : 0.0;
    }
  b::matrix<cd> irhoK = cd(0., 0.5)*prod(mrho, K);
  b::matrix<cd> DI = b::identity_matrix<cd>(_Nch) - irhoK;
  return DI;
}


void MmatrixK::Print() {
  std::cout << "------------ " << _Nch << "-channel K-matrix --------------" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Number of poles: " << _mass.size() << std::endl;
  std::cout << "Number of background terms: " << _bmass.size() << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
}


//b::zero_matrix<cd>(_Nch);  // (_Nch) K*=0.;
//  // add poles terms
//  const uint Npoles = _mass.size();
//  for (uint i = 0; i < Npoles; i++) {
//    double gpart[_Nch];
//    for (uint j = 0; j < _Nch; j++) gpart[j] = MParKeeper::gI()->get(_coupling[i*_Nch+j]);
//    double mass = MParKeeper::gI()->get(_mass[i]);
//    for (uint j = 0; j < _Nch; j++)
//      for (uint t = j; t < _Nch; t++)
//        K(j, t) += gpart[j]*gpart[t]/(mass*mass-s);
//  }
//  // add background terms
//  for (uint i = 0; i < _bmass.size(); i++) {
//    double gpart[_Nch];
//    for (uint j = 0; j < _Nch; j++) gpart[j] = MParKeeper::gI()->get(_bcs[i*_Nch+j]);
//    double bmass = MParKeeper::gI()->get(_bmass[i]);
//    for (uint j = 0; j < _Nch; j++)
//      for (uint t = j; t < _Nch; t++)
//        K(j, t) += gpart[j]*gpart[t]/(POW2(bmass)+s);
//  }
