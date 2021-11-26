// Copyright [2016] Misha Mikhasenko

// Coded from the paper
/*************************************************************************************/
/* @article{GarciaMartin:2011cn,                                                     */
/*       author         = "Garcia-Martin, R. and Kaminski, R. and Pelaez, J. R. and  */
/*                         Ruiz de Elvira, J. and Yndurain, F. J.",                  */
/*       title          = "{The Pion-pion scattering amplitude. IV: Improved         */
/*                         analysis with once subtracted Roy-like equations up to    */
/*                         1100 MeV}",                                               */
/*       journal        = "Phys. Rev.",                                              */
/*       volume         = "D83",                                                     */
/*       year           = "2011",                                                    */
/*       pages          = "074004",                                                  */
/*       doi            = "10.1103/PhysRevD.83.074004",                              */
/*       eprint         = "1102.2183",                                               */
/*       archivePrefix  = "arXiv",                                                   */
/*       primaryClass   = "hep-ph",                                                  */
/*       SLACcitation   = "%%CITATION = ARXIV:1102.2183;%%"                          */
/* }                                                                                 */
/*************************************************************************************/

// In addition,
// the amplitude above 1.4 GeV are exponentially smeared to 0
// see phi4 function.

#ifndef __GKPY_H__
#define __GKPY_H__

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>

#define uPOW2(v) ((v)*(v))
#define uKi(s, m) (sqrt((s)/4.-(m)*(m)))
#define uWM(s, s0) ((sqrt(s)-sqrt(s0-s))/(sqrt(s)+sqrt(s0-s)))

namespace waves {
  
  class GKPY {

   public:
    static double phi1(double s);
    static double phi2(double s);
    static double phi3(double s);
    static double phi4(double s);

    static double phi(double s);
    static std::complex<double> T(double s);
  
   private:
    static const double B0;
    static const double B1;
    static const double B2;
    static const double B3;
    static const double z0;
    static const double d0;
    static const double c;
    static const double B;
    static const double C;
    static const double D;

   private:
    static const double M0;
    static const double M2;

   private:
    static const double phi1_M0;
    static const double phi1_dM0;
    static const double phi3_M2;
    static const double phi3_dM2;

   private:
    static const double pi_mass;
    static const double k_mass;
    static const double eta_mass;
  };

};  // namespace waves

#endif  // __GKPY_H__
