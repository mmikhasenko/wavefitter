// Copyright [2016] Mikhail Mikhasenko

#include "deflib.h"

cd ChewMandelstam(cd s, double m1sq, double m2sq) {
  double sth  = POW2(sqrt(m1sq)+sqrt(m2sq));
  double sth2 = POW2(sqrt(m1sq)-sqrt(m2sq));
  cd value = -cd(0., 1.)*(-2./M_PI) * (-sqrt(sth-s)*sqrt(sth2-s)/s*log((sqrt(sth-s)+sqrt(sth2-s))/(2.*sqrt(sqrt(m1sq*m2sq)))) +
                                      (m1sq-m2sq)/(4.0*s)*log(m1sq/m2sq) -
                                      (m1sq+m2sq)/(4.0)*
                                      (
                                       (sth2 < 1e-3) ? 1./m1sq : log(m1sq/m2sq)/(m1sq-m2sq)
                                       )
                                       - 1./2.);
  return value;
}
