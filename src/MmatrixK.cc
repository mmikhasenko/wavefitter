// Copyright [2016] Mikhail Mikhasenko

#include <MmatrixK.h>
#include <MatrixInverse.h>

#include <deflib.h>

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::vector;
using boost::numeric::ublas::identity_matrix;
using boost::numeric::ublas::zero_matrix;
using boost::numeric::ublas::symmetric_matrix;
using boost::numeric::ublas::upper;
using boost::numeric::ublas::lower;

MmatrixK::MmatrixK(int Nch,
                   int Npoles,
                   std::vector<int> amap,
                   matrix<cd> (*background)(double s, double *pars) ) :
  // private stuff
  _Nch(Nch), _fph(Nch), _Np(Npoles), _bkgr(background), _amap(amap),
  // matrixes and vectors
  _T(Nch, Nch), _alpha(Nch), _A(Nch), last_s(-11.11),
  // parameters
  _mass(Npoles), _coupling(Nch*Npoles),
  // production mechanism
  _Aprod(Nch),
  // minor stuff
  need_to_recalculate_A(true), need_to_recalculate_T(true) {
  for (int i = 0; i < _Nch; i++) _fph[i]   = [](double s) -> double {return 1./(8*M_PI);};
  for (int i = 0; i < _Nch; i++) _Aprod[i] = [](double s, const double*) -> cd {return 1.;};

  int tot = 0; std::for_each(_amap.begin(), _amap.end(), [&tot](int i) {tot+=i;});
  _Na = tot;
  _apar = std::vector<double>(_Na);

  // Fill parameters map
  for (int ipole = 0; ipole < Npoles; ipole++) {
    // generate name
    std::ostringstream mname;
    mname << "m" << ipole + 1;
    map_to_pars[mname.str()] = &(_mass.data()[ipole]);
    for (int jch = 0; jch < _Nch; jch++) {
      std::ostringstream gname; gname << char(103+ipole + ((ipole >= 7) ? 1 : 0))
                                      << jch+1;
      map_to_pars[gname.str()] = &(_coupling.data()[_Nch*ipole+jch]);
    }
  }
  std::cout << "MmatrixK instance is created!" << std::endl;
  std::cout << "map of parameters is created(" << map_to_pars.size() << ")!" << std::endl;
}

void MmatrixK::setPars(const double *pars) {
  for (int i = 0; i < _Np; i++)
    if (_mass[i] != pars[i]) {
      need_to_recalculate_T = true;
      _mass[i] = pars[i];
    }

  for (int i = 0; i < _Np; i++) for (int j = 0; j < _Nch; j++)
      if (_coupling[_Nch*i+j] != pars[_Np + (_Nch*i+j)]) {
        need_to_recalculate_T = true;
        _coupling[_Nch*i+j] = pars[_Np + (_Nch*i+j)];
      }

  // chech alphas
  for (int i = 0; i < _Na; i++)
    if (_apar[i] != pars[i+_Np*(_Nch+1)]) {
      need_to_recalculate_A = true;
      _apar[i] = pars[i+_Np*(_Nch+1)];
    }
}

void MmatrixK::calculateT(double s) {

  // clear K
  symmetric_matrix<cd, upper> K = zero_matrix<cd>(_Nch);  // (_Nch) K*=0.;
  for (int i = 0; i < _Np; i++) {
    symmetric_matrix<cd, upper> km(_Nch);
    double *gpart = &_coupling[i*_Nch];
    // for(int j=0;j<_Nch;j++) std::cout << gpart[j] << " "; std::cout << std::endl;
    for (int j = 0; j < _Nch; j++) for (int t = j; t < _Nch; t++)
      km(j, t) = gpart[j]*gpart[t]/(_mass[i]*_mass[i]-s);
    K += km;
  }
  // ph.sp.matrix
  symmetric_matrix<cd, upper> mrho(_Nch);
  // possibly add some background
  for (int i = 0; i < _Nch; i++)
    for (int j = 0; j < _Nch; j++)
      mrho(i, j) = (i == j) ? _fph[i](s) : 0.0;

  matrix<cd> irhoK = cd(0, 1)*prod(mrho, K);
  matrix<cd> din = identity_matrix<cd>(_Nch) - irhoK;

  bool sing = false;
  matrix<cd> din_inv = gjinverse(din, sing);
  // std::cout << din_inv << std::endl;
  if (sing) {std::cerr << "ERROR: SINGULAR" << std::endl;}  // exit(); }

  // finally
  _T = prod(K, din_inv);
  // change the flag
  last_s = s; need_to_recalculate_T = false;
}

void MmatrixK::calculateA(double s) {
  // calculate T only if it is neccesary
  if (s != last_s || need_to_recalculate_T) calculateT(s);
  for (int i = 0;  i < _Nch; i++) _alpha(i) = _Aprod[i](s, &_apar[ (i == 0) ? 0 : _amap[i-1] ]);
  // std::cout << "-----------!!!!!!----------------"  << std::endl;
  _A = prod(_T, _alpha);
  // change the flag
  last_s = s; need_to_recalculate_A = false;
}

const vector<cd>& MmatrixK::getA(double s, const double *pars) {
  setPars(pars);
  if (s != last_s || need_to_recalculate_A || need_to_recalculate_T) calculateA(s);
  return _A;
}
const vector<cd>& MmatrixK::getA(double s) {
  if (s != last_s || need_to_recalculate_A || need_to_recalculate_T) calculateA(s);
  return _A;
}

cd MmatrixK::getA(int i, double s, const double *pars) {
  if (i >= _Nch || i < 0) {std::cerr << "Error <MmatrixK::getA>: i > _Nch || i < 0!" << std::endl; return 0.;}
  return getA(s, pars)(i);
}

cd MmatrixK::getA(int i, double s) {
  if (i >= _Nch || i < 0) {std::cerr << "Error <MmatrixK::getA>: i > _Nch || i < 0!" << std::endl; return 0.;}
  return getA(s)(i);
}

void MmatrixK::setPhSp(int i, double (&ph)(double s)) {
  if (i >= _Nch || i < 0) {std::cerr << "Error <MmatrixK::setPhSp>: i > _Nch || i < 0!" << std::endl; return;}
  _fph[i] = &ph;
}

void MmatrixK::setProd(int i, cd (*prod)(double, const double *)) {
  if (i >= _Nch || i < 0) {std::cerr << "Error <MmatrixK::setPhSp>: i > _Nch || i < 0!" << std::endl; return;}
  _Aprod[i] = prod;
}


void MmatrixK::Print() {
  std::cout << "------------ " << _Nch << "-channel K-matrix --------------" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Number of poles: " << _Np << std::endl;
  std::cout << "Number of production parameters: " << _Na << std::endl;
  std::cout << "Total number of parameters: " << getNpar() << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
}
