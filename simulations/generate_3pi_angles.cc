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

#include "M3bodyAngularBasis.h"


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
    
    if (!Math::changeAngularBasis( s1,  costheta1,  phi1,  costheta23,  phi23,
                                  &s3, &costheta3, &phi3, &costheta12, &phi12,
                                  m1sq, m2sq, m3sq, s)) continue;

    /**********************************************************************************************/
    /* Calculation of some quantities based on the angles *****************************************/
    /**********************************************************************************************/

    /**********************************************************************************************/
    /* deck-(23) */
    t = -0.1-gRandom->Exp(1./12.);  // BRANCH
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

    /**********************************************************************************************/
    /* Partial Waves */

    tout.Fill();
  }
  tout.Write();
  fout.Close();
  
  return 0;
}
