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
#include "TH2D.h"

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
  // waves.resize(2);
  uint Nwaves = waves.size();
  std::cout << "waves.size() = " << waves.size() << "\n";
  for (auto &w : waves)
    std::cout << w.index << ": " << w.title << " " << w.J << " "
              << (w.parity ? "+" : "-") << " " << w.M << " "
              << (w.pos_refl ? "+" : "-") << " " << w.S << " " << w.L << "\n";
  // isobars
  double m1sq = POW2(PI_MASS), m2sq = POW2(PI_MASS), m3sq = POW2(PI_MASS);
  double R = 5.;
  MIsobar rho(RHO_MASS, RHO_WIDTH, sqrt(m2sq), sqrt(m3sq), 1, R);
  MIsobar f2(F2_MASS, F2_WIDTH, sqrt(m2sq), sqrt(m3sq), 2, R);
  // functions
  std::function<cd(double)> isobarT[4];
  isobarT[0] = [&](double s)->cd {
    return waves::GKPY::T(s) * sqrt(2.) / 3.;
  };
  isobarT[1] = [&](double s)->cd {
    return rho.T(s) * 1. / sqrt(2.);
  };
  double BrF2pipi = 0.845;
  isobarT[2] = [&](double s)->cd {
    return f2.T(s) * sqrt(2.) / 3. * BrF2pipi;
  };
  isobarT[3] = [&](double s)->cd {
    return 0.;
  };

  // result
  TH2D *hr[Nwaves], *hi[Nwaves], *hint[Nwaves], *hphi[Nwaves];
  const uint NbinsX = 100;
  const uint NbinsY = 100;
  const std::pair<double, double> rangeX = std::make_pair(0.5, 2.5);
  const std::pair<double, double> rangeY = std::make_pair(2 * PI_MASS, 2.4);
  for (uint i = 0; i < Nwaves; i++) {
    hr[i] = new TH2D(TString::Format("h2r%d", i + 1),
                     TString::Format("Expansion coeff, real part, %s",
                                     waves[i].title.c_str()),
                     NbinsX, rangeX.first, rangeX.second, NbinsY, rangeY.first,
                     rangeY.second);
    hi[i] = new TH2D(TString::Format("h2i%d", i + 1),
                     TString::Format("Expansion coeff, imag part, %s",
                                     waves[i].title.c_str()),
                     NbinsX, rangeX.first, rangeX.second, NbinsY, rangeY.first,
                     rangeY.second);
    hint[i] = new TH2D(
        TString::Format("h2int%d", i + 1),
        TString::Format("Expansion coeff, absolute part x diff ph.sp., %s",
                        waves[i].title.c_str()),
        NbinsX, rangeX.first, rangeX.second, NbinsY, rangeY.first,
        rangeY.second);
    hphi[i] = new TH2D(
        TString::Format("h2phi%d", i + 1),
        TString::Format("Expansion coeff., phase, %s", waves[i].title.c_str()),
        NbinsX, rangeX.first, rangeX.second, NbinsY, rangeY.first,
        rangeY.second);
  }
  double E_BEAM_LAB = 190.;
  double mAsq = POW2(PI_MASS), mBsq = POW2(PROT_MASS), mDsq = POW2(PROT_MASS);
  double stot = POW2(PROT_MASS) + POW2(PI_MASS) + 2 * PROT_MASS * E_BEAM_LAB;
  double t = -0.1, mtRsq = POW2(PI_MASS);
  for (uint bx = 0; bx < NbinsX; bx++) {
    double en =
        rangeX.first + (rangeX.second - rangeX.first) / NbinsX * (bx + 0.5);
    double s = POW2(en);
    std::cout << "---> bin no." << bx + 1 << ", e = " << en << "\n";
    for (uint by = 0; by < NbinsY; by++) {
      double m23 =
          rangeY.first + (rangeY.second - rangeY.first) / NbinsY * (by + 0.5);
      double s1 = POW2(m23);
      // phase space jacobian
      double phsp =
          (en > m23 + PI_MASS)
              ? 1. / (8 * M_PI) * sqrt(LAMBDA(s, s1, POW2(PI_MASS))) / s
              : 0.0;
      double phsp1 = (m23 > 2 * PI_MASS)
                         ? 1. / (8 * M_PI) *
                               sqrt(LAMBDA(s1, POW2(PI_MASS), POW2(PI_MASS))) /
                               s1
                         : 0.0;
      double jac = 1. / (2 * M_PI) * phsp * phsp1 * 1. / POW2(4 * M_PI);

      // loop over waves
      for (uint w = 0; w < Nwaves; w++) {
        if (waves[w].S == -7)
          continue;
        uint S1 = (waves[w].S > 0) ? waves[w].S : 0;
        cd val(0., 0.);
        if (en > m23 + sqrt(m1sq)) {
          // calculate phase space
          val = MAscoli::getProjectedReducedDeck(
              waves[w].J, waves[w].M, (waves[w].pos_refl == waves[w].parity),
              waves[w].L, s1, S1, R, s, t, mtRsq, stot, mAsq, mBsq, mDsq,
              POW2(PI_MASS));
          // multiply to production clebsch, Blatt-Weisskopf
          val *= isobarT[S1](s1);
          // reflectivity basis
        }
        // set bin content
        hr[w]->SetBinContent(bx + 1, by + 1, real(val));
        hi[w]->SetBinContent(bx + 1, by + 1, imag(val));
        hint[w]->SetBinContent(bx + 1, by + 1, std::norm(val) * jac);
        hphi[w]->SetBinContent(bx + 1, by + 1, arg(val));
      }
    }
  }
  TFile *fout = new TFile(fout_name, "RECREATE");
  for (uint w = 0; w < Nwaves; w++) {
    hr[w]->Write();
    hi[w]->Write();
    hint[w]->Write();
    hphi[w]->Write();
  }
  std::cout << "File " << fout->GetName() << " have been completed!\n";
  fout->Close();
  return 0;
}

void fill_wavepull(const char *wave_fname, std::vector<wave> *waves) {
  uint index = 0;
  for (uint J = 0; J <= 4; J++) {
    for (uint M = 0; M <= J; M++) {
      for (int refl = 1; refl >= -1; refl -= 2) {
        for (uint S = 0; S <= 2; S++) {
          uint Lmin = (J > S) ? J - S : S - J;
          uint Lmax = S + J;
          for (uint L = Lmin; L <= Lmax; L++) {
            wave w;
            w.J = J;
            w.M = M;
            w.pos_refl = (refl == 1) ? true : false;
            w.S = S;
            w.L = L;
            w.parity =
                (S + L + 1) % 2 == 0; // true = positive, false = negarive
            if (M == 0) {
              int factor = 1 - (w.pos_refl ? 1 : -1) * (J % 2 == 0 ? 1 : -1) *
                               (w.parity ? 1 : -1);
              if (factor == 0)
                continue;
            }
            w.index = ++index;
            w.title = TString::Format("1-(%d%c+)%d%c(S=%d)#pi %d-wave", J,
                                      (w.parity ? '+' : '-'), M, (w.pos_refl ? '+' : '-'),
                                      S, L);
            waves->push_back(w);
          }
        }
      }
    }
  }
}
