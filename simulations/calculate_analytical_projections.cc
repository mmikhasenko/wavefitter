// Copyright [2017] Misha Mikhasenko
// Discription:
//   The program calculate projections of deck analytically, meaning all but one
// integral
//   are performed analytically, the latest one is calculated in the
// GetReducedDeck function

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
// #include "TText.h"
#include "TH1D.h"

#include "MDeck.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "waves.h"
#include "mintegrate.h"
#include "M3bodyAngularBasis.h"
#include "MAscoli.h"

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

#define SHIFT 7

void fill_wavepull(const char *wave_fname, std::vector<wave> *waves);

int main(int ac, char *av[]) {
  if (ac < 2) {
    std::cerr << "Usage: ./calculate_analytical_projections [FOUT_NAME]";
    return 0;
  }
  const char *fout_name = av[1];
  std::vector<wave> waves;
  fill_wavepull("/localhome/mikhasenko/results/pwa_3pi/wavelist_formated.txt",
                &waves);
  uint Nwaves = waves.size();
  std::cout << "waves.size() = " << waves.size() << "\n";
  for (auto &w : waves)
    std::cout << w.index << ": " << w.title << " " << w.J << " "
              << (w.parity ? "+" : "-") << " " << w.M << " "
              << (w.pos_refl ? "+" : "-") << " " << w.S << " " << w.L << "\n";
  // isobars
  std::vector<std::pair<uint, double> > isobars(SHIFT + 3 + 1);
  isobars[SHIFT + 3] = std::make_pair(3, 1.69);
  isobars[SHIFT + 2] = std::make_pair(2, F2_MASS);
  isobars[SHIFT + 1] = std::make_pair(1, RHO_MASS);
  isobars[SHIFT + 0] = std::make_pair(0, 0.6);  // sigma
  isobars[SHIFT - 1] = std::make_pair(0, 0.99); // f0(980)
  isobars[SHIFT - 2] = std::make_pair(0, 1.5);  // f0(1500)
  isobars[0] = std::make_pair(0, 0.3);          // FLAT
  // isospin clebsch
  double BrF2pipi = 0.845; // branching ratio f2->pipi
  std::vector<double> clebsch = { sqrt(2. / 3), 1, sqrt(2. / 3 * BrF2pipi), 1 };

  // for out mode;
  MIsobar rho(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
  rho.setIntU();
  MIsobar f2(F2_MASS, F2_WIDTH, PI_MASS, PI_MASS, 2, 5.);
  f2.setIntU();
  MIsobar rho3(1.69, 0.16, PI_MASS, PI_MASS, 1, 5.);
  rho3.setIntU();
  MIsobarPiPiS pipiS;
  pipiS.setIntU();
  // MIsobar *iso[] = {&pipiS, &rho, &f2, &rho3};
  MIsobar f980(0.99, 0.04, PI_MASS, PI_MASS, 0);
  f980.setIntU();
  MIsobar f1500(1.504, 0.11, PI_MASS, PI_MASS, 0);
  f1500.setIntU();
  // MIsobar *iso_scalars[] = {&pipiS, &f980, &f1500};

  std::vector<double> nS(SHIFT + 3 + 1);
  nS[SHIFT + 3] = rho3.IntU();
  nS[SHIFT + 2] = f2.IntU();
  nS[SHIFT + 1] = rho.IntU();
  nS[SHIFT + 0] = pipiS.IntU();
  nS[SHIFT - 1] = f980.IntU();
  nS[SHIFT - 2] = f1500.IntU();
  nS[0] = 1.;

  // result
  TH1D *hresr[Nwaves], *hresi[Nwaves], *hint[Nwaves];
  for (uint i = 0; i < Nwaves; i++) {
    hresr[i] = new TH1D(TString::Format("hr%d", i + 1),
                        TString::Format("Expansion coeff, real part, %s",
                                        waves[i].title.c_str()),
                        100, 0.5, 2.5);
    hresi[i] = new TH1D(TString::Format("hi%d", i + 1),
                        TString::Format("Expansion coeff, imag part, %s",
                                        waves[i].title.c_str()),
                        100, 0.5, 2.5);
    hint[i] = new TH1D(TString::Format("h%d", i + 1),
                       TString::Format("Intensity %s", waves[i].title.c_str()),
                       100, 0.5, 2.5);
  }
  double E_BEAM_LAB = 190.;
  double mAsq = POW2(PI_MASS), mBsq = POW2(PROT_MASS), mDsq = POW2(PROT_MASS);
  double stot = POW2(PROT_MASS) + POW2(PI_MASS) + 2 * PROT_MASS * E_BEAM_LAB;
  double m1sq = POW2(PI_MASS), R = 5., t = -0.1, mtRsq = POW2(PI_MASS);
  const uint Nbins = 100;
  for (uint b = 0; b < Nbins; b++) {
    double en = 0.5 + (2.5 - 0.5) / Nbins * (b + 0.5);
    double s = POW2(en);

    // loop over waves
    for (uint w = 0; w < Nwaves; w++) {
      if (waves[w].S == -7)
        continue;
      uint S1 = isobars[waves[w].S + SHIFT].first;
      double m23 = isobars[waves[w].S + SHIFT].second;
      double cl = clebsch[S1];
      double val = 0., phsp = 0., hL = 0.;
      if (sqrt(s) > m23 + sqrt(m1sq)) {
        // calculate phase space
        double psq = LAMBDA(s, POW2(m23), POW2(PI_MASS)) / (4 * s);
        phsp = 1. / (8 * M_PI) * sqrt(4 * psq / s);
        hL = pow(R * R * psq / (1 + R * R * psq), waves[w].L / 2.);
        psq *= hL;
        // calculate amplitude
        // std::cout << "w.index = " << waves[w].index << ", " << b << "\n";
        // std::cout << waves[w].J << ", " << waves[w].M << ", " << waves[w].L
        // << ", "
        //           << POW2(m23) << ", " << S1 << ", " <<  R << ", " << s << ",
        // "
        //           << t << ", " << mtRsq << ", " << stot << ", " << mAsq << ",
        // "
        //           << mBsq << ", " << mDsq << ", " << POW2(PI_MASS) << "\n";
        val = MAscoli::getProjectedReducedDeck(
            waves[w].J, waves[w].M, waves[w].L, POW2(m23), S1, R, s, t, mtRsq,
            stot, mAsq, mBsq, mDsq, POW2(PI_MASS));
        // std::cout << "val = " << val << "\n";
        // multiply to production clebsch, Blatt-Weisskopf
        double qsq =
            LAMBDA(POW2(m23), POW2(PI_MASS), POW2(PI_MASS)) / (4 * POW2(m23));
        double hS = pow(R * R * qsq / (R * R * qsq + 1.), S1 / 2.);
        val *= POW2(cl) / sqrt(2.) * hS / hL * nS[SHIFT + waves[w].S];
      }
      // set bin content
      hresr[w]->SetBinContent(b + 1, val);
      hresi[w]->SetBinContent(b + 1, 0);
      hint[w]->SetBinContent(b + 1, POW2(val) * phsp * hL);
    }
    /*************************************/
  }
  TFile *fout =
      new TFile(fout_name, "RECREATE");
  for (uint w = 0; w < Nwaves; w++) {
    hresr[w]->Write();
    hresi[w]->Write();
    hint[w]->Write();
  }
  std::cout << "File " << fout->GetName() << " have been completed!\n";
  fout->Close();
  return 0;
}

void fill_wavepull(const char *wave_fname, std::vector<wave> *waves) {
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
      n.J = 0;
      n.M = 0;
      n.S = -7;
      n.L = 0;
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
