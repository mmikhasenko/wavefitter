// Copyright [2016] Mikhail Mikhasenko

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include "MIsobarChannel.h"
#include "deflib.h"
#include "mintegrate.hh"

#define SCALEX_FOR_CUT 2.
#define SCALEY_FOR_CUT 10.0

MIsobarChannel::MIsobarChannel(const MIsobar &iso,
                               double m3,
                               int iL, double iR) :
  MChannel(iL, iR), _iso(iso), ltable(0) {;}


void MIsobarChannel::makeLookupTable(double from, double to, uint Npoints) {
  ltable.resize(Npoints);
  for (uint i = 0; i < Npoints; i++) {
    double s = from + (to-from)/(Npoints-1)*i;
    ltable[i].first = s;
    ltable[i].second = CalculateQuasiTwoBody(s);
  }
}

double MIsobarChannel::CalculateQuasiTwoBody(double s) const {
  if (s < POW2(sqrt(_iso.sth())+m3)) return 0;
  std::function<double(double)> drho = [&](double s12)->double {
    return 1./(2*M_PI)*_iso.U(s12)*RHO(s, s12, POW2(m3));
  };
  return integrate(drho, _iso.sth(), POW2(sqrt(s)+m3));
}

cd MIsobarChannel::CalculateQuasiTwoBodyStright(cd s) const {
  std::function<cd(double)> drho = [&](double t)->cd {
    cd s12 = _iso.sth()+(POW2(sqrt(s)+m3)-_iso.sth())*t;
    return 1./(2*M_PI)*_iso.U(s12)*RHO_PI(s, s12, m3*m3);
  };
  return (POW2(sqrt(s)+m3)-_iso.sth())*cintegrate(drho, 0, 1);
}

cd MIsobarChannel::CalculateQuasiTwoBodyEdge(cd s) const {
  std::function<cd(double)> drho = [&](double t)->cd {
    double th1 = _iso.sth();
    cd th2 = POW2(sqrt(s)-m3);
    cd thM(th1 + (real(th2)-th1)/SCALEX_FOR_CUT,
          imag(th2)/SCALEY_FOR_CUT);
    if (t < 1.) {
      cd s12 = th1+(thM-th1)*t;
      return 1./(2*M_PI)*_iso.U(s12)*RHO_PI(s, s12, m3*m3)  *  (thM-th1);
    } else {
      cd s12 = thM+(th2-thM)*(t-1);
      return 1./(2*M_PI)*_iso.U(s12)*RHO_PI(s, s12, m3*m3)  *  (th2-thM);
    }
  };
  return cintegrate(drho, 0, 2);
}

double MIsobarChannel::InterpolateQuasiTwoBody(double s) const {
  if (!ltable.size()) return -1;
  if (s >= ltable[ltable.size()-1].first) return RHO(s, POW2(_iso.GetM()), m3*m3);
  return getvalue(s, ltable);
}
