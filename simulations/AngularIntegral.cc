// Copyright [2017] Misha Mikhasenko

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "constants.h"
#include "deflib.h"

#include "TRandom.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TF1.h"

#include "M3bodyAngularBasis.h"

double fs1distr(double *x, double *par) {
  double s1 = x[0];
  double m1sq = par[0];
  double m2sq = par[1];
  double m3sq = par[2];
  double s = par[3];
  if ((s1 < (sqrt(m3sq)+sqrt(m2sq))*(sqrt(m3sq)+sqrt(m2sq))) ||
      (s1 > (sqrt(s)-sqrt(m1sq))*(sqrt(s)-sqrt(m1sq)))) return 0.;
  return sqrt(LAMBDA(s, s1, m1sq)*LAMBDA(s1, m2sq, m3sq))/s1;
}

#define EnergySq(p) ((p)[0]*(p)[0]+(p)[1]*(p)[1]+(p)[2]*(p)[2]+(p)[4]*(p)[4])

#define zBoostedEnergy(p, beta) (1./sqrt(1.-beta*beta)*((p)[3]+beta*(p)[2]))
#define zBoostedPz(p, beta) (1./sqrt(1.-beta*beta)*(beta*(p)[3]+(p)[2]))

#define invMassSq(p1, p2) ((p1[3]+p2[3])*(p1[3]+p2[3])-(p1[0]+p2[0])*(p1[0]+p2[0])-(p1[1]+p2[1])*(p1[1]+p2[1])-(p1[2]+p2[2])*(p1[2]+p2[2]))
#define pVecSq(p) ((p)[0]*(p)[0]+(p)[1]*(p)[1]+(p)[2]*(p)[2])

double find_phi_in_rotated_frame(double *p, double costheta, double phi) {
  return atan2(-sin(phi)*p[0]+cos(phi)*p[1],
               costheta*(cos(phi)*p[0]+sin(phi)*p[1])-sqrt(1-costheta*costheta)*p[2]);
}

void rotateZ(double *p, double phi) {
  double px = p[0], py = p[1];
  p[0] = cos(phi)*px - sin(phi)*py;
  p[1] = sin(phi)*px + cos(phi)*py;
}
void rotateY(double *p, double theta) {
  double px = p[0], pz = p[2];
  p[2] = cos(theta)*pz - sin(theta)*px;
  p[0] = sin(theta)*pz + cos(theta)*px;
}

int main() {
  // constants
  double m1sq = POW2(PI_MASS);
  double m2sq = POW2(PI_MASS);
  double m3sq = POW2(PI_MASS);
  double s = POW2(0.6);

  // random generator
  TF1 tfs1distr("t1", &fs1distr, POW2(sqrt(m2sq)+sqrt(m3sq)), POW2(sqrt(s)-sqrt(m1sq)), 4);
  double tfpars[] = {m1sq, m2sq, m3sq, s}; tfs1distr.SetParameters(tfpars);

  // plot for check
  TH2D hdalitz("hdalitz", "Dalitz-Plot;s1(GeV);s3(GeV)", 
               100, POW2(sqrt(m2sq)+sqrt(m3sq)), POW2(sqrt(s)-sqrt(m1sq)),
               100, POW2(sqrt(m1sq)+sqrt(m2sq)), POW2(sqrt(s)-sqrt(m3sq)));

  cd intgr(0., 0.);
  uint nLoop = 1e6;
  for (uint i = 0; i < nLoop; i++) {
    /****************************************************************/
    // generation
    /****************************************************************/
    // decay of main mass
    double costheta1 = 2.*gRandom->Rndm()-1.;
    double sintheta1 = sqrt(1-POW2(costheta1));
    double phi1 = M_PI*(2.*gRandom->Rndm()-1.);
    // decay of isobar
    double costheta23 = 2*gRandom->Rndm()-1.;
    double sintheta23 = sqrt(1-POW2(costheta23));
    double phi23 = M_PI*(2.*gRandom->Rndm()-1.);
    // isobar mass
    double s1 = tfs1distr.GetRandom();

    double p1 = sqrt(LAMBDA(s, s1, m1sq)/(4.*s));
    double p2pr = sqrt(LAMBDA(s1, m1sq, m2sq)/(4.*s1));
    double p1_v[] = {-p1*sintheta1*cos(phi1), -p1*sintheta1*sin(phi1), -p1*costheta1, 0.0, sqrt(m1sq)};
    p1_v[3] = sqrt(EnergySq(p1_v));

    double p2pr_v[] = { p2pr*sintheta23*cos(phi23),  p2pr*sintheta23*sin(phi23),  p2pr*costheta23, 0.0, sqrt(m2sq)};
    p2pr_v[3] = sqrt(EnergySq(p2pr_v));
    double p3pr_v[] = {-p2pr*sintheta23*cos(phi23), -p2pr*sintheta23*sin(phi23), -p2pr*costheta23, 0.0, sqrt(m3sq)};
    p3pr_v[3] = sqrt(EnergySq(p3pr_v));

    double p2_v[5], p3_v[5]; for (uint t = 0; t < 5; t++) { p2_v[t] = p2pr_v[t]; p3_v[t] = p3pr_v[t]; }
    // inverse boost and rotation
    double beta = p1/sqrt(p1*p1+s1);
    // p2 back to CMS
    p2_v[3] = zBoostedEnergy(p2pr_v, beta);
    p2_v[2] = zBoostedPz(p2pr_v, beta);
    rotateY(p2_v, acos(costheta1));
    rotateZ(p2_v, phi1);
    // p3 back to CMS
    p3_v[3] = zBoostedEnergy(p3pr_v, beta);
    p3_v[2] = zBoostedPz(p3pr_v, beta);
    rotateY(p3_v, acos(costheta1));
    rotateZ(p3_v, phi1);
    /****************************************************************/
    // check of the correctness
    /****************************************************************/
    double m_s3 = invMassSq(p1_v, p2_v);
    double m_s1 = invMassSq(p2_v, p3_v);
    double m_s2 = invMassSq(p3_v, p1_v);
    if (fabs(s1-m_s1) > 1e-3) std::cerr << "Error: something is wrong with s1!\n";
    if (fabs(p1_v[2]+p2_v[2]+p3_v[2]) > 1e-3) std::cerr << "Error: something is wrong with sum p [2]!\n";
    hdalitz.Fill(m_s1, m_s3);
    /****************************************************************/
    // calculations
    /****************************************************************/
    // variables to use
    // 23 system
    double m_phi1 = atan2(-p1_v[1], -p1_v[0]);
    double m_costheta1 = -p1_v[2]/sqrt(pVecSq(p1_v));
    double m_p2pr = sqrt(LAMBDA(m_s1, m2sq, m3sq)/(4.*m_s1));;
    double m_p1_rf23 = sqrt(LAMBDA(s, m_s1, m1sq)/(4.*m_s1));;
    double m_costheta23 = (m_s3-m1sq-m2sq-2*(m_s1+m2sq-m3sq)*(s-m_s1-m1sq)/(4.*m_s1))/(2*m_p2pr*m_p1_rf23);
    double m_phi23 = find_phi_in_rotated_frame(p2_v, m_costheta1, m_phi1);
    // std::cout << "--------- Mesurements(23) ----------\n";
    // std::cout << "m_costheta1: " << m_costheta1 << ", " << costheta1 << "\n";
    // std::cout << "m_phi1: " << m_phi1 << ", " << phi1 << "\n";
    // std::cout << "m_costheta23: " << m_costheta23 << ", " << costheta23 << "\n";
    // std::cout << "m_phi23: " << m_phi23 << ", " << phi23 << "\n";

    // 12 system
    double m_phi3 = atan2(-p3_v[1], -p3_v[0]);
    double m_costheta3 = -p3_v[2]/sqrt(pVecSq(p3_v));
    double m_p1pr = sqrt(LAMBDA(m_s3, m1sq, m2sq)/(4.*m_s3));;
    double m_p3_rf12 = sqrt(LAMBDA(s, m_s3, m3sq)/(4.*m_s3));;
    double m_costheta12 = (m_s2-m3sq-m1sq-2*(m_s3+m1sq-m2sq)*(s-m_s3-m3sq)/(4.*m_s3))/(2*m_p1pr*m_p3_rf12);
    double m_phi12 = find_phi_in_rotated_frame(p1_v, m_costheta3, m_phi3);
    // std::cout << "--------- Mesurements(12) ----------\n";
    // std::cout << "m_costheta3: " << m_costheta3 << "\n";
    // std::cout << "m_phi3: " << m_phi3 << "\n";
    // std::cout << "m_costheta12: " << m_costheta12 << "\n";
    // std::cout << "m_phi12: " << m_phi12 << "\n";
    // std::cout << "------------------------------------\n\n";
    
    double m_theta1 = acos(m_costheta1);
    double m_theta23 = acos(m_costheta23);
    double m_theta3 = acos(m_costheta3);
    double m_theta12 = acos(m_costheta12);

    // D-functions
    uint J = 4, M = 0, L = 3, S = 1;
    // FunctionZ(z)
    cd Zf1 = Math::ZJMLS(J, M, L, S, m_theta1, m_phi1, m_theta23, m_phi23);
    cd Zf3 = Math::ZJMLS(J, M, L, S, m_theta3, m_phi3, m_theta12, m_phi12);
    // cd amp = Zf3;
    cd amp = (Zf1 + Zf3) / sqrt(2.);
    intgr += amp*conj(amp);
  }
  intgr /= nLoop;
  intgr *= 4.*M_PI * 4.*M_PI;
  std::cout << "Result of integration is " << intgr << "\n";

  // using lib-function
  uint J = 4, M = 0, L = 3, S = 1;

  std::function<double(double, double, double, double, double,
                       double, double, double, double, double)>
    myf = [&, J, M, L, S](double s1, double theta1, double phi1, double theta23, double phi23,
                          double s3, double theta3, double phi3, double theta12, double phi12)->double{                   
    cd Zf1 = Math::ZJMLS(J, M, L, S, theta1, phi1, theta23, phi23);
    cd Zf3 = Math::ZJMLS(J, M, L, S, theta3, phi3, theta12, phi12);
    // cd amp = Zf3;
    cd amp = (Zf1 + Zf3) / sqrt(2.);
    return 4.*M_PI * 4.*M_PI * norm(amp);};
  double vV = Math::integrate3bphs(myf,
                                   10000,
                                   POW2(0.6), POW2(PI_MASS), POW2(PI_MASS), POW2(PI_MASS));
  std::cout << "Result of integration with library is " << vV << "\n";

  TCanvas can("can");
  hdalitz.Draw("col");
  can.Print("/tmp/AngularIntegral.pdf");

  return 0;
}


