// Copyright [2017] Misha Mikhasenko
// A simple program to generate 3 particles phase space decay angles

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TF2.h"
// #include "TH1D.h"

#include "MDeck.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "waves.h"
#include "M3bodyAngularBasis.h"
#include "TText.h"


#include "constants.h"
#include "mintegrate.h"
#include "deflib.h"

typedef struct {
  uint index;
  uint J;
  bool parity;
  uint M;
  bool pos_refl;
  int S;
  uint L;
  double threshold;
  std::string title;
} wave;

void fill_wavepull(const char* wave_fname, std::vector<wave> *waves);


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

  /* prepare for the calculation of Deck-like amplitude */

  // Deck-(23), Deck-(12)
  cd decklike1, decklike3;
  double decklike1_real, decklike1_imag, decklike3_real, decklike3_imag;
  tout.Branch("decklike1_real", &decklike1_real);
  tout.Branch("decklike1_imag", &decklike1_imag);
  tout.Branch("decklike3_real", &decklike3_real);
  tout.Branch("decklike3_imag", &decklike3_imag);

  double s;
  tout.Branch("s", &s);

  /* prepare for the calculation of PWs */
  // read wavepull
  std::vector<wave> waves;
  fill_wavepull("/localhome/mikhasenko/results/pwa_3pi/wavelist_formated.txt", &waves);
  uint nWaves = waves.size();
  std::cout << "waves.size() = " << waves.size() << "\n";
  for (auto & w : waves)
    std::cout << w.index << ": " << w.title << " "
              << w.J << " " << (w.parity ? "+" : "-") << " " << w.M << " " << (w.pos_refl ? "+" : "-")
              << " " << w.S << " " << w.L << "\n";
  // branches
  cd amp[nWaves][2];
  double amp_real[nWaves][2], amp_imag[nWaves][2];
  for (uint w = 0; w < nWaves; w++) {
    tout.Branch(TString::Format("amp%d_frame1_real", waves[w].index), &amp_real[w][0]);
    tout.Branch(TString::Format("amp%d_frame1_imag", waves[w].index), &amp_imag[w][0]);
    tout.Branch(TString::Format("amp%d_frame3_real", waves[w].index), &amp_real[w][1]);
    tout.Branch(TString::Format("amp%d_frame3_imag", waves[w].index), &amp_imag[w][1]);
  }

  // Create isobars
  MIsobar rho(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.); rho.setIntU();
  MIsobar  f2(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.); f2.setIntU();
  MIsobar rho3(1.69, 0.16, PI_MASS, PI_MASS, 1, 5.); rho3.setIntU();
  MIsobarPiPiS pipiS; pipiS.setIntU();
  MIsobar *iso[] = {&pipiS, &rho, &f2, &rho3};
  MIsobar f980(0.99, 0.04, PI_MASS, PI_MASS, 0); f980.setIntU();
  MIsobar f1500(1.504, 0.11, PI_MASS, PI_MASS, 0); f1500.setIntU();
  MIsobar *iso_scalars[] = {&pipiS, &f980, &f1500};

  // main loop
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

    for (uint bose = 0; bose < 2; bose++) {
      /*************************************/
      double sI        = bose?         s3 :         s1;
      double costhetaI = bose?  costheta3 :  costheta1;
      double phiI      = bose?       phi3 :       phi1;
      double costheta  = bose? costheta12 : costheta23;
      double phi       = bose?      phi12 :      phi23;

      double thetaI = acos(costhetaI);
      double theta  = acos(costheta);
      // loop over waves
      for (uint w = 0; w < nWaves; w++) {
        cd iso_shape = waves[w].S == -7 ? 1. :
          (waves[w].S > 0 ?
           iso        [ waves[w].S]->ToneVertex(sI) :
           iso_scalars[-waves[w].S]->ToneVertex(sI));

        double R = 5.;
        double qsq = LAMBDA(s, sI, POW2(PI_MASS))/(4*s);
        double BlattWeisskopf = pow(R*R*qsq/(1.+R*R*qsq), waves[w].L/2.);

        amp[w][bose] = Math::ZJMLS_refl(waves[w].J, waves[w].M,
                                        (waves[w].pos_refl == waves[w].parity),  // (-1)*(-1) = (+1)*(+1) = true, otherwise is false
                                        waves[w].L, (waves[w].S > 0 ? waves[w].S : 0),
                                        thetaI, phiI, theta, phi) * iso_shape * BlattWeisskopf;
        amp_real[w][bose] = real(amp[w][bose]);
        amp_imag[w][bose] = imag(amp[w][bose]);
      }
    }

    tout.Fill();
  }
  tout.Write();
  for (uint w = 0; w < nWaves; w++) {
    TText tx(0., 0., waves[w].title.c_str()); tx.SetName(TString::Format("t%d", waves[w].index));
    tx.Write();
  }
  fout.Close();
  
  return 0;
}



void fill_wavepull(const char* wave_fname, std::vector<wave> *waves) {
  std::ifstream fin(wave_fname);
  if (!fin.is_open()) {
    std::cerr << "Error: can not find the file!\n";
    return;
  }
  std::string line;
  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    wave n;
    iss >> n.index;
    iss >> n.title;
    if (n.title == "FLAT") {
      n.J = 0; n.M = 0; n.S = -7; n.L = 0;
      n.parity = false;
      n.pos_refl = true;
      waves->push_back(n);
      continue;
    }
    iss >> n.J;
    std::string parity;
    iss >> parity;
    n.parity = (parity == "+");
    iss >> n.M;
    std::string epsilon;
    iss >> epsilon;
    n.pos_refl = (epsilon == "+");
    iss >> n.S;
    iss >> n.L;
    waves->push_back(n);
  }
}
