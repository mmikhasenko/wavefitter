// Copyright [14.07.2015] Misha Mikhasenko

#include "MHelicityDeck.h"
#include "mintegrate.h"

#include "TWigner.h"
#include "MIsobar.h"

#include "TWigner.h"
#include "TMath.h"
#include "Math/SpecFuncMathMore.h"

#define NUMERICAL_LIMIT_IM 1e-6

double MHelicityDeck::getDeck(double m1sq, double m2sq, double m3sq, double m4sq, double mtsq,
                      double s, double z,
                      uint Sp, int lamP,
                      uint S , int lamS,
                      double R) {
  if (s == POW2(sqrt(m3sq)+sqrt(m4sq))) s+=1e-5;
  if (s  < POW2(sqrt(m3sq)+sqrt(m4sq))) { std::cerr << "Warning<MHelicityDeck::getDeck> s = "
                                                    << s << " < sth = "
                                                    << POW2(sqrt(m3sq)+sqrt(m4sq))
                                                    << ", diff is " << s-POW2(sqrt(m3sq)+sqrt(m4sq))
                                                    << "\n"; return 0; }
  // std::cout << "~!!!" << s << "\n";
  // calculate dependent variables
  double t = m1sq + m3sq - ((s + m3sq - m4sq)*(s + m1sq - m2sq))/(2.*s) +
    sqrt(LAMBDA(s, m1sq, m2sq)*LAMBDA(s, m3sq, m4sq))/(2.*s)*z;
  cd CosChi2 = ((s + m2sq - m1sq)*(t + m2sq - m4sq) - 2.*m2sq*(m2sq + m3sq - m1sq - m4sq)) /
    sqrt(LAMBDA(s, m1sq, m2sq)*LAMBDA(t, m2sq, m4sq));
  cd CosChi3 = ((s + m3sq - m4sq)*(t + m3sq - m1sq) - 2.*m3sq*(m2sq + m3sq - m1sq - m4sq)) /
    sqrt(LAMBDA(s, m3sq, m4sq)*LAMBDA(t, m1sq, m3sq));
  cd Chi2 = acos(CosChi2);
  cd Chi3 = acos(CosChi3);
  double ClebschP = pow(2., Sp/2.)*TMath::Factorial(Sp)/sqrt(TMath::Factorial(2*Sp));
  double ClebschS = pow(2., S/2.)*TMath::Factorial(S)/sqrt(TMath::Factorial(2*S));

  // (* main fucnction *)
  cd f = ClebschP *
    // pow(cd(POW2(R)*POW2(R)*LAMBDA(t, m1sq, m3sq))/POW2(1. - POW2(R)*t), S/2.) *
    pow(cd(LAMBDA(t, m1sq, m3sq))/(LAMBDA(t, m1sq, m3sq) - 4.*t/POW2(R)), S/2.) *
    // pow(cd(LAMBDA(t, m1sq, m3sq))/(POW2(R)*LAMBDA(t, m1sq, m3sq) + 1.), S/2.) *
    1./(mtsq - t) *
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // temperary fix with real part
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Math::WignerD(2*Sp, 0, 2*lamP, real(Chi2)) *
      Math::WignerD(2*S,  0, 2*lamS, real(Chi3)) *
    // pow(cd(POW2(R)*POW2(R)*LAMBDA(t, m2sq, m4sq))/POW2(1. - t*POW2(R)), Sp/2.) /
    pow(cd(LAMBDA(t, m2sq, m4sq))/(LAMBDA(t, m2sq, m4sq) - 4.*t/POW2(R)), Sp/2.) /
    // pow(cd(LAMBDA(t, m2sq, m4sq))/(POW2(R)*LAMBDA(t, m2sq, m4sq) + 1.), Sp/2.) /
    pow(cd(m2sq), lamP/2.) * /* to suppress complexity */
    ClebschS;

  if (fabs(imag(f)) > NUMERICAL_LIMIT_IM) std::cerr << "Warning<MHelicityDeck::getDeck> imag part [b(s = "
                                                    << s << ", z = "
                                                    << z << ")] = " << f << "!= 0" << std::endl;
  return real(f);
}

double MHelicityDeck::getProjectedDeck(double m1sq, double m2sq, double m3sq, double m4sq, double mtsq,
                               double s,
                               uint J, uint L,
                               uint Sp, int lamP,
                               uint S, double R) {
  // integrate
  int minSJ = (S < J) ? S : J;
  double int_val =
    integrate([&, m1sq, m2sq, m3sq, m4sq, mtsq, s, J, L, Sp, lamP, S, minSJ](double z)->double{
      double val = 0;
      for (int lamS = -minSJ; lamS <= minSJ; lamS++)
        val +=
          ROOT::Math::wigner_3j(2*L, 2*S, 2*J, 2*0, 2*lamS, -2*lamS) *
          getDeck(m1sq, m2sq, m3sq, m4sq, mtsq, s, z, Sp, lamP, S, -lamS, R) *
          Math::WignerD(2*J, 2*lamP, 2*lamS, acos(z));
      return val;
      }, -1., 1.);
  int_val *= sqrt(2*J+1) * /*because of Clebsch and 3j*/
    sqrt(M_PI*(2*L+1)); /*because of normalisation*/

  return int_val;
}

void MHelicityDeck::setLookupTable(const std::vector<std::pair<double, double> > &ltable) {
  _ltable.resize(ltable.size());
  std::copy(ltable.begin(), ltable.end(), _ltable.begin());
}

void MHelicityDeck::makeLookupTable(const MIsobar &iso, double m3,
                            double from, double to, uint Npoints) {
  // first create integrand
  double sr;
  std::function<double(double)> integrand = [&](double s3)->double {
    double pD = getProjection(sr, s3);
    return iso.U(s3)*POW2(pD)*1./(2*M_PI) * 1./(8*M_PI)*sqrt(LAMBDA(sr, s3, POW2(m3)))/sr;
  };
  std::function<double(double)> drho = [&](double s3)->double {
    return iso.U(s3)*1./(2*M_PI) * 1./(8*M_PI)*sqrt(LAMBDA(sr, s3, POW2(m3)))/sr;
  };
  // then create and fill table
  _ltable.resize(Npoints);
  for (uint i = 0; i < Npoints; i++) {
    sr = from + (to-from)/(Npoints-1)*i;
    if (sr <= POW2(sqrt(iso.sth()) + m3)) { _ltable[i] = std::make_pair(sr, 0); continue; }
    double intD = integrate(integrand, iso.sth(), POW2(sqrt(sr)-m3));
    double intRho = integrate(drho, iso.sth(), POW2(sqrt(sr)-m3));
    _ltable[i] = std::make_pair(sr, sqrt(intD/intRho));
    std::cout << "--> " << sr << ", " << _ltable[i].second << std::endl;
  }
}

double MHelicityDeck::getPrecalculated(double s) const {
  if (!_ltable.size()) { std::cerr << "Warning<MHelicityDeck::getPrecalculated> TABLE is empty\n"; return 0; }
  if (s < _ltable[0].first) { std::cerr << "Warning<MHelicityDeck::getPrecalculated> s < FIRST\n"; return 0; }
  if (s < _ltable[_ltable.size()-1].first) return getvalue(s, _ltable);
  return (_ltable[_ltable.size()-1].second * pow(_ltable[_ltable.size()-1].first, 1+abs(_lamP))) / pow(s, 1+abs(_lamP));
}

void MHelicityDeck::Print() const {
  std::cout << " instance of Deck class: J = " << _J
            << ", Sp = " << _Sp
            << ", lamP = " << _lamP
            << ", S = " << _S
            << ", L = " << _L
            << ", R = " << _R << "\n";
}
