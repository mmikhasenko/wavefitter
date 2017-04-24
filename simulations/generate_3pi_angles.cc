// Copyright [2017] Misha Mikhasenko
// A simple program to generate 3 particles phase space decay angles

#include <iostream>

#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TF2.h"
// #include "TH1D.h"

#include "MDeck.h"

#include "constants.h"
#include "mintegrate.h"
#include "deflib.h"

double fs1(double *x, double *par) {
  // double s = par[0];
  double m1sq = par[0];
  double m2sq = par[1];
  double m3sq = par[2];
  double s1 = x[0];
  double s = x[1];
  if ((s1 < POW2(sqrt(m2sq)+sqrt(m3sq)) || (s1 > POW2(sqrt(s)-sqrt(m1sq))))) return 0.;
  return sqrt(LAMBDA(s, s1, m1sq))/s*sqrt(LAMBDA(s1, m2sq, m3sq))/s1;
}

int main(int ac, char **av) {

  const char *fout_name = (ac >= 2) ? av[1] : "/tmp/test_generate_3pi_angles.root";
  const uint nEvents = (ac >= 3) ? atoi(av[2]) : 1e5;
  const double S_MIN = (ac >= 5)  ? atof(av[3])*atof(av[3]) : (9*PI_MASS*PI_MASS);
  const double S_MAX = (ac >= 5)  ? atof(av[4])*atof(av[4]) : 9.;
  const uint seed = (ac >= 6)  ? atoi(av[5]) : 0;
  std::cout << "Output: " << fout_name << "\n";
  std::cout << "Number of events: " << nEvents << "\n";
  std::cout << "M_{3pi} range: [ " << sqrt(S_MIN) << ", " << sqrt(S_MAX) << " ]\n";
  std::cout << "Seed is " << seed << "\n";

  const double m1sq = POW2(PI_MASS);
  const double m2sq = POW2(PI_MASS);
  const double m3sq = POW2(PI_MASS);

  gRandom->SetSeed(seed);

  TF2 tfs1("tfs1", fs1,
           POW2(2*PI_MASS), POW2(sqrt(S_MAX)-PI_MASS),
           S_MIN, S_MAX, 3);
  tfs1.SetParameter(0, m1sq);
  tfs1.SetParameter(1, m2sq);
  tfs1.SetParameter(2, m3sq);
  tfs1.SetNpx(500);
  tfs1.SetNpy(500);

  TFile fout(fout_name, "RECREATE");
  TTree tout("angles", "angles");
  // (23)-frame
  double s1, costheta1, phi1, costheta23, phi23;
  tout.Branch("s1", &s1);
  tout.Branch("costheta1", &costheta1);
  tout.Branch("phi1", &phi1);
  tout.Branch("costheta23", &costheta23);
  tout.Branch("phi23", &phi23);

  // (12)-frame
  double s3, costheta3, phi3, costheta12, phi12;
  tout.Branch("s3", &s3);
  tout.Branch("costheta3", &costheta3);
  tout.Branch("phi3", &phi3);
  tout.Branch("costheta12", &costheta12);
  tout.Branch("phi12", &phi12);

  // general
  double E_BEAM = 190;  // GeV
  double s0 = 2*E_BEAM*PROT_MASS + POW2(PROT_MASS) + POW2(PI_MASS);
  double t;
  tout.Branch("s0", &s0);
  tout.Branch("t", &t);

  double R = 5;

  // Deck-(23), Deck-(12)
  cd decklike1, decklike3;
  double decklike1_real, decklike1_imag, decklike3_real, decklike3_imag;
  tout.Branch("decklike1_real", &decklike1_real);
  tout.Branch("decklike1_imag", &decklike1_imag);
  tout.Branch("decklike3_real", &decklike3_real);
  tout.Branch("decklike3_imag", &decklike3_imag);

  double s;
  tout.Branch("s", &s);

  for (uint e = 0; e < nEvents; e++) {
    // from the function
    tfs1.GetRandom2(s1, s);  // BRANCH

    costheta1 = 2*gRandom->Rndm()-1.;  // BRANCH
    phi1 = M_PI*(2*gRandom->Rndm()-1.);  // BRANCH
    costheta23 = 2*gRandom->Rndm()-1.;  // BRANCH
    phi23 = M_PI*(2*gRandom->Rndm()-1.);  // BRANCH
    
    // calculate s3 in (23) frame
    s3 = m1sq + m2sq +  // BRANCH
      // 2*(  E2 * E1 -
      2*( (s1 + m2sq - m3sq)/(2*sqrt(s1)) * (s-m1sq-s1)/(2*sqrt(s1)) -
          // |p2|*cos(theta23) * (-|p1|)
          sqrt(LAMBDA(s1, m2sq, m3sq)/(4*s1))*costheta23 * (-sqrt(LAMBDA(s, m1sq, s1)/(4*s1))) );
    if (s3 != s3) continue;  // nan check

    // caluclate p3 in (23) frame
    double p23bu = sqrt(LAMBDA(s1, m2sq, m3sq)/(4*s1));
    double p3_in23[] = {(s1 + m3sq - m2sq)/sqrt(4*s1),
                        -p23bu*sqrt(1-POW2(costheta23))*cos(phi23),
                        -p23bu*sqrt(1-POW2(costheta23))*sin(phi23),
                        -p23bu*costheta23};
    double gamma1 = (s + s1 - m1sq)/sqrt(4*s*s1);
    double beta1 = sqrt(1.-1./POW2(gamma1));
    double p3_b[] = {gamma1*(p3_in23[0]+beta1*p3_in23[3]),
                     p3_in23[1],
                     p3_in23[2],
                     gamma1*(beta1*p3_in23[0]+p3_in23[3])};
    double p3_rot[4];
    { double ct = costheta1, st = sqrt(1.-POW2(costheta1)), cp = cos(phi1), sp = sin(phi1);
      // Rz(phi1) * Ry(theta1) * p3_boost
      p3_rot[0] = p3_b[0];
      p3_rot[1] = cp*ct* p3_b[1] + (-sp)* p3_b[2] + cp*st* p3_b[3];
      p3_rot[2] = sp*ct* p3_b[1] +   cp * p3_b[2] + sp*st* p3_b[3];
      p3_rot[3] =  - st* p3_b[1] +    0 * p3_b[2] +    ct* p3_b[3];
    }
    // calculate Omega12
    // BRANCH
    costheta3 = - p3_rot[3]/sqrt(POW2(p3_rot[1])+POW2(p3_rot[2])+POW2(p3_rot[3]));
    phi3 = atan2(-p3_rot[2], -p3_rot[1]);  // BRANCH

    // BRANCH
    costheta12 = (s1 - m2sq - m3sq - 2* (s3+m2sq-m1sq)/(2*sqrt(s3)) * (s-m3sq-s3)/(2*sqrt(s3))) /
      (2 * sqrt(LAMBDA(s3, m1sq, m2sq)/(4*s3)) * sqrt(LAMBDA(s, m3sq, s3)/(4*s3)) );
    if (costheta12 != costheta12) continue;  // nan check

    double n1[] = {0.,
                   -sqrt(1.-POW2(costheta1))*cos(phi1),
                   -sqrt(1.-POW2(costheta1))*sin(phi1),
                   -costheta1};
    // Ry(-theta12) * Rz(-phi12) * n1;
    double n1_rot[4];
    { double ct = costheta3, st = sqrt(1.-POW2(costheta3)), cp = cos(phi3), sp = sin(phi3);
      // Rz(phi23) * Ry(theta23) * p3_boost
      n1_rot[0] = n1[0];
      n1_rot[1] = cp*ct* n1[1] + sp*ct* n1[2] + (-st)* n1[3];
      n1_rot[2] = (-sp)* n1[1] +   cp * n1[2] +    0 * n1[3];
      n1_rot[3] = cp*st* n1[1] + sp*st* n1[2] +    ct* n1[3];
    }
    phi12 = atan2(n1_rot[2], n1_rot[1]);  // BRANCH

    /* deck-(23) */
    t = -0.1-gRandom->Exp(1./12.);
    decklike1 = MDeck::getAmplitude(costheta1, phi1,
                                    s1, R,
                                    costheta23 , phi23,
                                    s, t,
                                    POW2(PI_MASS),
                                    s0,
                                    POW2(PI_MASS), POW2(PROT_MASS), POW2(PROT_MASS),
                                    POW2(PI_MASS), POW2(PI_MASS), POW2(PI_MASS));
    // if (!(e%1000))std::cout << decklike1_real << "\n";
    decklike1_real = real(decklike1);  // BRANCH
    decklike1_imag = imag(decklike1);  // BRANCH
    /* deck-(12) */
    decklike3 = MDeck::getAmplitude(costheta3, phi3,
                                    s3, R,
                                    costheta12 , phi12,
                                    s, t,
                                    POW2(PI_MASS),
                                    s0,
                                    POW2(PI_MASS), POW2(PROT_MASS), POW2(PROT_MASS),
                                    POW2(PI_MASS), POW2(PI_MASS), POW2(PI_MASS));
    decklike3_real = real(decklike3);  // BRANCH
    decklike3_imag = imag(decklike3);  // BRANCH

    tout.Fill();
  }
  tout.Write();
  fout.Close();
  
  return 0;
}
