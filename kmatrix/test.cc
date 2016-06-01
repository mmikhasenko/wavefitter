// Copyright [2016] Mikhail Mikhasenko

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
typedef std::complex<double> cd;

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace bnu = boost::numeric::ublas;

bnu::matrix<cd> calculate_matrix(cd s, double *par);
template <typename Type> Type rhoF2Pi(Type s);
template <typename Type> Type rhoRhoPi(Type s);

template<class T> bool InvertMatrix(const bnu::matrix<T>& input, bnu::matrix<T>& inverse);
template<class T> bool InvertMatrix(const bnu::symmetric_matrix<T,bnu::upper>& input, bnu::matrix<T>& inverse);

#define PI_MASS   0.139570
#define F2_MASS   1.27
#define RHO_MASS  0.7755

#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>

int main() {
  double initual_pars[] = {1.6, 1.88,
                           1.8, 0.9, 0.2,
                           0.15, 0.62, 2.1,
                           0.0001};
  const int Npoints = 100;
  const double sDown = (PI_MASS+F2_MASS)*(PI_MASS+F2_MASS)+0.05;
  const double sUp = 6.25;

  const cd alpha[3] = {cd(1., -0.4), cd(1, 0.1), cd(4., -0.5)};
  cd amps[3][Npoints];
  double mX[Npoints];
  for (int i=0; i < Npoints; i++) {
    double s = sDown + (sUp-sDown)/(Npoints-1)*i; mX[i] = sqrt(s);
    bnu::matrix<cd> m = calculate_matrix(s, initual_pars);
    for (int j=0; j < 3; j++) {
      amps[j][i] = 0;
      for (int t=0; t < 3; t++) amps[j][i] += alpha[t]*m(j, t);
    }
  }
  TCanvas canva("can", "Canva", 0, 0, 1000, 800);
  canva.Divide(3, 3);
  // intensities
  TGraph *intens[3];
  for (int i=0; i < 3; i++) {
    double Y[Npoints]; for (int t=0; t < Npoints; t++) Y[t] = pow(abs(amps[i][t]), 2);
    intens[i] = new TGraph(Npoints, mX, Y);
    intens[i]->SetName(TString::Format("intens_%d", i));
    intens[i]->SetTitle(TString::Format("intens_%d", i));
  }
  // interference
  TGraph *interf[3];
  std::pair<int, int> comb[3];
  comb[0] = std::make_pair(0, 1);
  comb[1] = std::make_pair(0, 2);
  comb[2] = std::make_pair(1, 2);
  for (int i=0; i < 3; i++) {
    double Y[Npoints]; for (int t=0; t < Npoints; t++) Y[t] = arg(amps[comb[i].first][t]*conj(amps[comb[i].second][t]));
    interf[i] = new TGraph(Npoints, mX, Y);
    interf[i]->SetName(TString::Format("interf_%d", i));
    interf[i]->SetTitle(TString::Format("interf_%d", i));
  }  
  canva.cd(1); intens[0]->Draw();
  canva.cd(5); intens[1]->Draw();
  canva.cd(9); intens[2]->Draw();
  canva.cd(2); interf[0]->Draw();
  canva.cd(3); interf[1]->Draw();
  canva.cd(6); interf[2]->Draw();
  canva.SaveAs("c1.pdf");

  return 0;
}

bnu::matrix<cd> calculate_matrix(cd s, double *par) {
  double m1 = par[0];
  double m2 = par[1];
  // couplings to I
  double g1[3] = {par[2], par[3], par[4]};
  // couplings to II
  double g2[3] = {par[5], par[6], par[7]};
  // background
  double bgrd = par[8];

  // std::cout << "(m1*m1-s) = " << (m1*m1-s) << std::endl;
  // fill the firsr matrix
  bnu::symmetric_matrix<cd, bnu::upper> km1(3, 3);
  for (unsigned i = 0; i < km1.size1 (); ++i)
    for (unsigned j = i; j < km1.size2 (); ++j)
      km1(i, j) = g1[i]*g1[j];
  km1 *= 1./(m1*m1-s);
  // std::cout << "km1, " << km1 << std::endl;
  // fill the second matrix
  bnu::symmetric_matrix<cd, bnu::upper> km2(3, 3);
  for (unsigned i = 0; i < km2.size1 (); ++i)
    for (unsigned j = i; j < km2.size2 (); ++j)
      km2(i, j) = g2[i]*g2[j];
  km2 *= 1./(m2*m2-s);
  // std::cout << "km2, " << km2 << std::endl;
  // add to km
  bnu::symmetric_matrix<cd, bnu::upper> km = km1+km2;
  // std::cout << "km, " << km << std::endl;
  // possibly add some background
  for (unsigned i = 0; i < km.size1 (); ++i)
    for (unsigned j = i; j < km.size2 (); ++j)
      km(i, j) += bgrd;
  // std::cout << km << std::endl;
  // ph.sp.matrix
  bnu::symmetric_matrix<cd, bnu::upper> mrho(3, 3);
  // possibly add some background
  for (unsigned i = 0; i < mrho.size1 (); ++i)
    for (unsigned j = i; j < mrho.size2 (); ++j)
      mrho(i, j) = 0.0;
  mrho(0, 0) = rhoF2Pi(s);
  mrho(1, 1) = rhoF2Pi(s);
  mrho(2, 2) = rhoRhoPi(s);
  // std::cout << "\nrho " << mrho << std::endl;
  // inverse and and phas
  bnu::identity_matrix<cd> uni(3);
  bnu::matrix<cd> irhoK = cd(0, 1)*prod(mrho, km);
  // std::cout << "\nirhoK " << irhoK << std::endl;
  bnu::matrix<cd> din = uni - irhoK;
  // std::cout << "\nmatrix " << din << std::endl;
  bnu::matrix<cd> din_inv(3, 3); InvertMatrix(din, din_inv);
  // std::cout << "inv matrix " << din_inv << std::endl;
  bnu::matrix<cd> T = prod(km, din_inv);
  return T;
  return km;
}

template <typename Type> Type rhoF2Pi(Type s) {
  double  mPlus2 = (PI_MASS+F2_MASS)*(PI_MASS+F2_MASS);
  double mMinus2 = (PI_MASS-F2_MASS)*(PI_MASS-F2_MASS);
  return 1./(8*M_PI)*sqrt((1.-mPlus2/s)*(1.-mMinus2/s));
}

template <typename Type> Type rhoRhoPi(Type s) {
  double  mPlus2 = (PI_MASS+RHO_MASS)*(PI_MASS+RHO_MASS);
  double mMinus2 = (PI_MASS-RHO_MASS)*(PI_MASS-RHO_MASS);
  return 1./(8*M_PI)*sqrt((1.-mPlus2/s)*(1.-mMinus2/s));
}

template<class T> bool InvertMatrix(const bnu::matrix<T>& input, bnu::matrix<T>& inverse) {
  typedef bnu::permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  bnu::matrix<T> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());
  // perform LU-factorization
  int res = lu_factorize(A, pm);
  if (res != 0) return false;
  // create identity matrix of "inverse"
  inverse.assign(bnu::identity_matrix<T> (A.size1()));
  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);
  return true;
}

template<class T> bool InvertMatrix(const bnu::symmetric_matrix<T,bnu::upper>& input, 
                                    bnu::matrix<T>& inverse) {
  typedef bnu::permutation_matrix<std::size_t> pmatrix;
  bnu::matrix<T> A(input);
  pmatrix pm(A.size1());
  int res = lu_factorize(A, pm);
  if (res != 0) return false;
  inverse.assign(bnu::identity_matrix<T> (A.size1()));
  lu_substitute(A, pm, inverse);
  return true;
}
