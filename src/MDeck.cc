// Copyright [14.07.2015] Misha Mikhasenko

#include "MDeck.h"
#include "mintegrate.hh"
#include "deflib.h"

#include "MIsobar.h"

#define NUMERICAL_LIMIT_IM 1e-6

double MDeck::getDeck(double m1sq, double m2sq, double m3sq, double m4sq, double mtsq,
                      double s, double z,
                      uint Sp, int lamP,
                      uint S , int lamS,
                      double R) {
  if (s == POW2(sqrt(m3sq)+sqrt(m4sq))) s+=1e-5;
  if (s  < POW2(sqrt(m3sq)+sqrt(m4sq))) { std::cerr << "Warning<MDeck::getDeck> s<sth\n"; return 0; }
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
    pow(cd(POW2(R)*POW2(R)*LAMBDA(t, m1sq, m3sq))/POW2(1. - POW2(R)*t), S/2.) *
    // pow(cd(LAMBDA(t, m1sq, m3sq))/(LAMBDA(t, m1sq, m3sq) - 4.*t/POW2(R)), S/2.) *
    // pow(cd(LAMBDA(t, m1sq, m3sq))/(POW2(R)*LAMBDA(t, m1sq, m3sq) + 1.), S/2.) *
    1./(mtsq - t) *
    rpwa::dFunction<cd>(2*Sp, 0, 2*lamP, Chi2) *
    rpwa::dFunction<cd>(2*S,  0, 2*lamS, Chi3) *
    pow(cd(POW2(R)*POW2(R)*LAMBDA(t, m2sq, m4sq))/POW2(1. - t*POW2(R)), Sp/2.) /
    // pow(cd(LAMBDA(t, m2sq, m4sq))/(LAMBDA(t, m2sq, m4sq) - 4.*t/POW2(R)), Sp/2.) /
    // pow(cd(LAMBDA(t, m2sq, m4sq))/(POW2(R)*LAMBDA(t, m2sq, m4sq) + 1.), Sp/2.) /
    pow(cd(m2sq), lamP/2.) * /* to suppress complexity */
    ClebschS;

  if (fabs(imag(f)) > NUMERICAL_LIMIT_IM) std::cerr << "Warning<MDeck::getDeck> imag part [b] = " << f << "!= 0" << std::endl;
  return real(f);
}

double MDeck::getProjectedDeck(double m1sq, double m2sq, double m3sq, double m4sq, double mtsq,
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
          rpwa::dFunction<double> (2*J, 2*lamP, 2*lamS, acos(z));
      return val;
      }, -1., 1.);
  int_val *= sqrt(2*J+1) * /*because of Clebsch and 3j*/
    sqrt(M_PI*(2*L+1)); /*because of normalisation*/

  return int_val;
}

void MDeck::makeLookupTable(const MIsobar &iso, double m3,
                              double from, double to, uint Npoints) {
  // first create integrand
  double sr;
  std::function<double(double)> integrand = [&](double s3)->double {
    double pD = getProjection(sr, s3);
    return iso.U(s3)*POW2(pD)*1./(2*M_PI);
  };
  // then create and fill table
  _ltable.resize(Npoints);
  for (uint i = 0; i < Npoints; i++) {
    sr = from + (to-from)/(Npoints-1)*i;
    _ltable[i] = std::make_pair(sr,
                                integrate(integrand, iso.sth(), POW2(sqrt(sr)-m3)));
    std::cout << "--> " << sr << ", " << _ltable[i].second << std::endl;
  }
}

double MDeck::getPrecalculated(double s) {
  if (!_ltable.size()) { std::cerr << "Warning<MDeck::getPrecalculated> TABLE is empty\n"; return 0; }
  if (s < _ltable[0].first) { std::cerr << "Warning<MDeck::getPrecalculated> s < FIRST\n"; return 0; }
  if (s < _ltable[_ltable.size()-1].first) return getvalue(s, _ltable);
  return (_ltable[_ltable.size()-1].second * pow(_ltable[_ltable.size()-1].first, 1+abs(_lamP))) / pow(s, 1+abs(_lamP));
}
