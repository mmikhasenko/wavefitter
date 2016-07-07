// Copyright [2016] Mikhail Mikhasenko

#include "MChannel.h"
#include "mintegrate.h"

void MChannel::makeDisperseLookupTable(double from, double to, uint Npoints) {
  dtable.resize(Npoints);
  for (uint i = 0; i < Npoints; i++) {
    double s = from + (to-from)/(Npoints-1)*i;
    dtable[i].first = s;
    dtable[i].second = DisperceRhoLtilda(s);
  }
}

cd MChannel::DisperceRhoLtilda(double s) const {
  double RH = rho(s);
  double DC = DumpC(s);
  std::function<cd(double)> drhoLtilda = [&](double up)->cd{
    double sp = 1./up;
    return (rho(sp)*DumpC(sp)-RH*DC)/(sp*(sp-s-cd(0., 1.e-6))) * 1./(up*up);  // minus in taken in integration range
  };
  cd integral = s*cintegrate(drhoLtilda, 0.0, 1./sth());
  cd add_term = RH*DC*cd(-log(fabs(s/sth()-1.)), M_PI);
  return cd(0, -1./(M_PI))*(integral+add_term);
}

cd MChannel::DisperceRhoLtilda(cd s) const {
  std::function<cd(double)> drhoLtilda = [&](double up)->cd{
    double sp = 1./up;
    return rho(sp)*DumpC(sp)/(sp*(sp-s-cd(0., 1.e-6))) * 1./(up*up);  // minus in taken in integration range
  };
  cd integral = cintegrate(drhoLtilda, 0.0, 1./sth());
  return s*cd(0, -1.)*integral/(M_PI);
}


cd MChannel::InterpolateRhoLtilda(double s) const {
  if (!dtable.size()) return -1.;
  if (s >= dtable[dtable.size()-1].first || s <= dtable[0].first) return 0.0;
  return getvalue(s, dtable);
}
