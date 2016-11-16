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

// function is not finished
// [not finished] cd MDeck::getAmplitude(double pAx, double pAy, double pAz, double mAsq,
// [not finished]                        double pBx, double pBy, double pBz, double mBsq,
// [not finished]                        double p1x, double p1y, double p1z, double m1sq,  // pi-
// [not finished]                        double p2x, double p2y, double p2z, double m2sq,  // pi+
// [not finished]                        double p3x, double p3y, double p3z, double m3sq,  // pi-
// [not finished]                        double mDsq) {
// [not finished]   // function is not finished
// [not finished]
// [not finished]   double eA = sqrt(POW2(pAx)+POW2(pAy)+POW2(pAz)+mAsq);
// [not finished]   double eB = sqrt(POW2(pBx)+POW2(pBy)+POW2(pBz)+mBsq);
// [not finished]   double e1 = sqrt(POW2(p1x)+POW2(p1y)+POW2(p1z)+m1sq);
// [not finished]   double e2 = sqrt(POW2(p2x)+POW2(p2y)+POW2(p2z)+m2sq);
// [not finished]   double e3 = sqrt(POW2(p3x)+POW2(p3y)+POW2(p3z)+m3sq);
// [not finished]
// [not finished]   // invariants
// [not finished]   double stot = POW2(eA+eB)-POW2(pAx+pBx)-POW2(pAy+pBy)-POW2(pAz+pBz);
// [not finished]   double s = POW2(e1+e2+e3)-POW2(p1x+p2x+p3x)-POW2(p1y+p2y+p3y)-POW2(p1z+p2z+p3z);
// [not finished]   double t = POW2(e1+e2+e3-eA)-POW2(p1x+p2x+p3x-pAx)-POW2(p1y+p2y+p3y-pAy)-POW2(p1z+p2z+p3z-pAz);
// [not finished]
// [not finished]   // GJ-frame quantities
// [not finished]   double pAGJ = sqrt(LAMBDA(s, t, mAsq)/(4*s));
// [not finished]   double eAGJ = (s+mAsq-t)/(2.*sqrt(s));
// [not finished]
// [not finished]   double mS1sq = POW2(e2+e3)-POW2(p2x+p3x)-POW2(p2y+p3y)-POW2(p2z+p3z);
// [not finished]   double mS3sq = POW2(e2+e1)-POW2(p2x+p1x)-POW2(p2y+p1y)-POW2(p2z+p1z);
// [not finished]   double RS1 = 5; double RS3 = 5;
// [not finished]
// [not finished]   /* Isobar (23) */
// [not finished]   // cos theta, phi
// [not finished]   double P23PA = (e2+e3)*eA - (p2x+p3x)*pAx - (p2y+p3y)*pAy - (p2z+p3z)*pAz;
// [not finished]   double e23GJ = (s+mS1sq-m1sq)/(2*sqrt(s));
// [not finished]   double p23GJ = sqrt(LAMBDA(s, mS1sq, m1sq)/(4*s));
// [not finished]   double costhetaS1 = (eAGJ*e23GJ-P23PA) / (pAGJ*p23GJ);
// [not finished]   double phiS1 = 0.0;
// [not finished]   // cos theta', phi
// [not finished]   double costheta_pr_S1 = 0.0;
// [not finished]   double phi_pr_S1 = 0.0;
// [not finished]
// [not finished]   getAmplitude(costhetaS1, phiS1,
// [not finished]                mS1sq, RS1,
// [not finished]                costheta_pr_S1, phi_pr_S1,
// [not finished]                s, t,
// [not finished]                POW2(PI_MASS),
// [not finished]                stot,
// [not finished]                mAsq, mBsq, mDsq,
// [not finished]                m1sq, m2sq, m3sq);
// [not finished]
// [not finished]   /* Isobar (12) */
// [not finished]   // cos theta, phi
// [not finished]   double P23PA = (e2+e3)*eA - (p2x+p3x)*pAx - (p2y+p3y)*pAy - (p2z+p3z)*pAz;
// [not finished]   double e23GJ = (s+mS1sq-m1sq)/(2*sqrt(s));
// [not finished]   double p23GJ = sqrt(LAMBDA(s, mS1sq, m1sq)/(4*s));
// [not finished]   double costhetaS3 = (eAGJ*e23GJ-P23PA) / (pAGJ*p23GJ);
// [not finished]   double phiS3 = 0.0;
// [not finished]   // cos theta', phi
// [not finished]   double costheta_pr_S3 = 0.0;
// [not finished]   double phi_pr_S3 = 0.0;
// [not finished]   getAmplitude(costhetaS3, phiS3,
// [not finished]                mS3sq, RS3,
// [not finished]                costheta_pr_S3, phi_pr_S3,
// [not finished]                s, t,
// [not finished]                POW2(PI_MASS),
// [not finished]                stot,
// [not finished]                mAsq, mBsq, mDsq,
// [not finished]                m1sq, m2sq, m3sq);
// [not finished]   return 0.0;
// [not finished]
// [not finished] }
