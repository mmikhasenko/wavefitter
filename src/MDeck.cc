// Copyright [14.07.2015] Misha Mikhasenko

#include "waves.h"
#include "deflib.h"

#include "MIsobar.h"
#include "MDeck.h"
#include "MAscoli.h"

#define NUMPRESICION 1e-8

cd MDeck::getAmplitude(double costheta, double phi,
                       double mS1sq, double RS1,
                       double costheta_pr, double phi_pr,
                       double wsq, double t,
                       double mtRsq,
                       double stot,
                       double mAsq, double mBsq, double mDsq,
                       double m1sq, double m2sq, double m3sq) {
  MIsobar rho(RHO_MASS, RHO_WIDTH, sqrt(m2sq), sqrt(m3sq), 1, RS1);
  MIsobar f2(RHO_MASS, RHO_WIDTH, sqrt(m2sq), sqrt(m3sq), 2, RS1);

  std::function<cd(double)> isobarT[3];
  isobarT[0] = [&](double s) -> cd { return waves::GKPY::T(s); };
  isobarT[1] = [&](double s) -> cd { return rho.T(s); };
  isobarT[2] = [&](double s) -> cd { return f2.T(s); };

  // loop over spin, helicity
  cd amp = 0.0;
  for (uint S = 0; S <= 2; S++) {
    cd deck_fixed_S_isobar_reduced = 0;
    for (int lamS = -S; lamS <= static_cast<int>(S); lamS++) {
      cd val = MAscoli::fullDeckTerm(costheta, phi,
                            mS1sq, S, lamS, RS1,
                            costheta_pr, phi_pr,
                            wsq, t,
                            mtRsq,
                            stot,
                            mAsq, mBsq, mDsq,
                            m1sq);
      deck_fixed_S_isobar_reduced += real(val);
    }
    if (imag(deck_fixed_S_isobar_reduced) > NUMPRESICION) {
      std::cerr << "Error<MDeck::getAmplitude> : imag = " << imag(deck_fixed_S_isobar_reduced) << " > 1e-8!!\n";
      return 0.0;
    }
    amp += deck_fixed_S_isobar_reduced * isobarT[S](mS1sq);
  }
  return amp;
}
