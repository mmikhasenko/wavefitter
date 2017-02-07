// Copyright [2017] Misha Mikhasenko


#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "M3bodyAngularBasis.h"

#include "constants.h"
#include "mintegrate.h"
#include "deflib.h"


int main() {
  // set parameters
  // MIsobar rho(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.); rho.setIntU();
  MIsobarPiPiS pipiS;  // pipiS.setIntU();
  MIsobar *iso = &pipiS;

  uint J = 4, M = 0, L = 3, S = 1;

  uint Npoints = 21;
  uint NintPoints = 10000;
  std::pair<double, double> range(0.5, 2.5);
  for (uint i = 0; i < Npoints; i++) {
    double s = POW2(range.first + (range.second-range.first)/(Npoints-1)*i);
    /***********************************************************************************/
    // Phhase space
    double phsp = integrate([s](double s1)->double{
        return sqrt(LAMBDA(s, s1, POW2(PI_MASS))*LAMBDA(s1, POW2(PI_MASS), POW2(PI_MASS)))/s1;
      }, 4*POW2(PI_MASS), POW2(sqrt(s)-PI_MASS)) / (2*M_PI*POW2(8*M_PI)*s);

    std::cout << "sqrt(s), Ph.Sp., Symm, Non-Symm: " << sqrt(s) << "  " << phsp << "  ";
    /***********************************************************************************/
    // Symmetrized
    std::function<double(double, double, double, double, double,
                         double, double, double, double, double)>
      Msq_symm = [&, J, M, L, S](double s1, double theta1, double phi1, double theta23, double phi23,
                                 double s3, double theta3, double phi3, double theta12, double phi12)->double{
      cd Zf1 = Math::ZJMLS(J, M, L, S, theta1, phi1, theta23, phi23);
      cd Zf3 = Math::ZJMLS(J, M, L, S, theta3, phi3, theta12, phi12);
      // cd amp = Zf3;
      cd amp = (iso->T(s1)*Zf1 + iso->T(s3)*Zf3) / sqrt(2.);
      return norm(amp);};
    double int_symm = Math::integrate3bphs(Msq_symm, NintPoints, s, POW2(PI_MASS), POW2(PI_MASS), POW2(PI_MASS));
    std::cout << int_symm << "  ";

    /***********************************************************************************/
    // Non symmetrized
    std::function<double(double, double, double, double, double,
                         double, double, double, double, double)>
      Msq_nonsymm = [&, J, M, L, S](double s1, double theta1, double phi1, double theta23, double phi23,
                                    double s3, double theta3, double phi3, double theta12, double phi12)->double{
      cd Zf1 = Math::ZJMLS(J, M, L, S, theta1, phi1, theta23, phi23);
      // cd amp = Zf3;
      cd amp = iso->T(s1)*Zf1;
      return norm(amp);};
    double int_nonsymm = Math::integrate3bphs(Msq_nonsymm, NintPoints, s, POW2(PI_MASS), POW2(PI_MASS), POW2(PI_MASS));
    std::cout << int_nonsymm << "\n";
  }

  return 0;
}
