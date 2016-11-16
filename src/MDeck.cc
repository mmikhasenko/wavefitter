// Copyright [14.07.2015] Misha Mikhasenko

#include <functional>

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
  MIsobar f2(F2_MASS, F2_WIDTH, sqrt(m2sq), sqrt(m3sq), 2, RS1);

  std::function<cd(double)> isobarT[3];
  isobarT[0] = [&](double s) -> cd { return waves::GKPY::T(s) * sqrt(2.)/3.; };
  isobarT[1] = [&](double s) -> cd { return rho.T(s); };
  isobarT[2] = [&](double s) -> cd { return f2.T(s) * sqrt(2.)/3.; };

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
      deck_fixed_S_isobar_reduced += val;
    }
    if (imag(deck_fixed_S_isobar_reduced) > NUMPRESICION) {
      std::cerr << "Error<MDeck::getAmplitude> : imag = " << imag(deck_fixed_S_isobar_reduced) << " > 1e-8!!\n";
      return 0.0;
    }
    amp += (2.*S+1.) * deck_fixed_S_isobar_reduced * isobarT[S](mS1sq);
  }
  return amp;
}

cd MDeck::getAmplitude(double p1x, double p1y, double p1z,  // pi-
                       double p2x, double p2y, double p2z,  // pi+
                       double t, double mtRsq, double stot,
                       double mAsq, double mBsq, double mDsq,
                       double m1sq, double m2sq, double m3sq) {
  cd amp = 0.0;

  // calculate p2 as -p1-p3
  double p3x = -p1x-p2x; double p3y = -p1y-p2y; double p3z = -p1z-p2z;

  double e1 = sqrt(POW2(p1x)+POW2(p1y)+POW2(p1z)+m1sq);
  double e2 = sqrt(POW2(p2x)+POW2(p2y)+POW2(p2z)+m2sq);
  double e3 = sqrt(POW2(p3x)+POW2(p3y)+POW2(p3z)+m3sq);

  // invariants
  double s = POW2(e1+e2+e3)-POW2(p1x+p2x+p3x)-POW2(p1y+p2y+p3y)-POW2(p1z+p2z+p3z);
  double mS3sq = POW2(e2+e1)-POW2(p2x+p1x)-POW2(p2y+p1y)-POW2(p2z+p1z);
  double RS3 = 5;  // double RS1 = 5;
  double P1P3 = e1*e3-p1x*p3x-p1y*p3y-p1z*p3z;

  /* Isobar (12) */
  // cos theta, phi
  double costhetaS3 = (p2z+p1z)/sqrt(POW2(p2x+p1x)+POW2(p2y+p1y)+POW2(p2z+p1z));
  double phiS3 = atan2(p2y+p1y, p2x+p1x);

  // cos theta', phi
  double e3prS3 = (s-mS3sq-m3sq)/(2*sqrt(mS3sq));
  double p3prS3 = sqrt(LAMBDA(s, mS3sq, m3sq)/(4*mS3sq));
  double p1prS3 = sqrt(LAMBDA(mS3sq, m2sq, m1sq)/(4*mS3sq));
  double e1prS3 = (mS3sq+m1sq-m2sq)/(2*sqrt(mS3sq));
  double costheta_pr_S3 = (P1P3-e3prS3*e1prS3)/(p3prS3*p1prS3);
  double sinthetaS3 = sqrt(1-POW2(costhetaS3));
  double phi_pr_S3 = atan2(p1y*cos(phiS3)-p1x*sin(phiS3),
                           p1x*costhetaS3*cos(phiS3)-p1z*sinthetaS3+p1y*costhetaS3*sin(phiS3));
  // std::cout << "costhetaS3 = " << costheta_pr_S3 << ", phiS3 = " << phi_pr_S3 << "\n";

  amp += getAmplitude(costhetaS3, phiS3,
                      mS3sq, RS3,
                      costheta_pr_S3, phi_pr_S3,
                      s, t,
                      mtRsq,
                      stot,
                      mAsq, mBsq, mDsq,
                      m1sq, m2sq, m3sq);
  return amp;
}

#include "TLorentzVector.h"

// function is not finished
uint MDeck::fromLabToGJ(double pAx, double pAy, double pAz, double mAsq,
                        double pBx, double pBy, double pBz, double mBsq,
                        double p1x, double p1y, double p1z, double m1sq,  // pi-
                        double p2x, double p2y, double p2z, double m2sq,  // pi+
                        double p3x, double p3y, double p3z, double m3sq,  // pi-
                        double mDsq) {
  // function is not finished

  TLorentzVector pi1_lv, pi2_lv, pi3_lv, beam_lv, trgt_lv;
  pi1_lv.SetXYZM(p1x, p1y, p1z, sqrt(m1sq));
  pi2_lv.SetXYZM(p2x, p2y, p2z, sqrt(m2sq));
  pi3_lv.SetXYZM(p3x, p3y, p3z, sqrt(m3sq));
  // beam and target
  beam_lv.SetXYZM(pAx, pAy, pAz, sqrt(mAsq));
  trgt_lv.SetXYZM(pBx, pBy, pBz, sqrt(mBsq));
  TLorentzVector reso_lv = pi1_lv + pi2_lv + pi3_lv;
  TLorentzVector recl_lv = beam_lv + trgt_lv - reso_lv;
  if (mDsq != 0. && recl_lv.M2() != mDsq) {
    std::cerr << "Error <MDeck::fromLabToGJ> you provided inconsistent information, recl_lv.M2() != mDsq!";
    return 0;
  }

  /*************************** Transformation to GJ frame ********************************/
  // boost and rotation to GJ frame
  // boost
  TVector3 bv = -reso_lv.BoostVector();
  pi1_lv.Boost(bv);
  pi2_lv.Boost(bv);
  pi3_lv.Boost(bv);
  beam_lv.Boost(bv);
  recl_lv.Boost(bv);
  // to check
  reso_lv.Boost(bv);

  // rotation beam to z
  TVector3 oldZ = beam_lv.Vect().Unit();
  pi1_lv .RotateZ(-beam_lv.Phi());  pi1_lv .RotateY(-beam_lv.Theta());
  pi2_lv .RotateZ(-beam_lv.Phi());  pi2_lv .RotateY(-beam_lv.Theta());
  pi3_lv .RotateZ(-beam_lv.Phi());  pi3_lv .RotateY(-beam_lv.Theta());
  recl_lv.RotateZ(-beam_lv.Phi());  recl_lv.RotateY(-beam_lv.Theta());
  // to check
  beam_lv.RotateZ(-beam_lv.Phi());  beam_lv.RotateY(-beam_lv.Theta());

  // rotation xy
  double phi = M_PI-recl_lv.Phi();  // M_PI is really important to check!!
  pi1_lv.RotateZ(phi);
  pi2_lv.RotateZ(phi);
  pi3_lv.RotateZ(phi);
  beam_lv.RotateZ(phi);
  // to check
  recl_lv.RotateZ(phi);
  /*****************************************************************************************/

  return 0;
}
