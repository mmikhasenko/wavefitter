// Copyright [08.08.2016] Misha Mikhasenko

#ifndef SRC_MASCOLI_H_
#define SRC_MASCOLI_H_

#include "deflib.h"

#include <iostream>
#define _USE_MATH_DEFINES

namespace MAscoli {

double getReducedDeck(double costheta, double phi, double mS1sq, uint S1,
                      int lamS1, double RS1, double wsq, double t, double mtRsq,
                      double stot, double mAsq, double mBsq, double mDsq,
                      double m1sq);

double getReducedDeckProjectionM(int M, double costheta,
                                 double mS1sq, uint S1, int lamS1, double RS1,
                                 double wsq, double t, double mtRsq,
                                 double stot, double mAsq, double mBsq,
                                 double mDsq, double m1sq);

double getProjectedReducedDeck(uint J, int M, uint L, double mS1sq, uint S1,
                               double RS1, double wsq, double t, double mtRsq,
                               double stot, double mAsq, double mBsq,
                               double mDsq, double m1sq);

double getProjectedReducedDeck(uint J, int M, bool pos_refl, uint L,
                               double mS1sq, uint S1, double RS1, double wsq,
                               double t, double mtRsq, double stot, double mAsq,
                               double mBsq, double mDsq, double m1sq);

double getProjectionJMSlam(uint J, int M, int lam, double mS1sq, uint S1,
                           double RS1, double wsq, double t, double mtRsq,
                           double stot, double mAsq, double mBsq, double mDsq,
                           double m1sq);

double upperPart(double costheta, double mS1sq, uint S1, int lamS1, double RS1,
                 double wsq, double t, double mtRsq, double mAsq, double m1sq);

double sPionProton(double costheta, double phi, double mS1sq, double wsq,
                   double t, double stot, double mAsq, double mBsq, double mDsq,
                   double m1sq);
double tPionIsobar(double costheta, double mS1sq, double wsq, double t,
                   double mAsq, double m1sq);

cd fullDeckTerm(double costheta, double phi, double mS1sq, uint S1, int lamS1,
                double RS1, double costheta_pr, double phi_pr, double wsq,
                double t, double mtRsq, double stot, double mAsq, double mBsq,
                double mDsq, double m1sq);

double psi(double costheta, double mS1sq, double wsq, double t, double mAsq,
           double m1sq);
};

#endif // SRC_MASCOLI_H_
