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
  uint S;
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
    std::cout << w.index << ": " << w.title << " " << w.J << " " << w.S << " "
              << "\n";
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
  TH2D *hint[Nwaves];
  const uint NbinsX = 100;
  const uint NbinsY = 100;
  const std::pair<double, double> rangeX = std::make_pair(0.5, 2.5);
  const std::pair<double, double> rangeY = std::make_pair(2 * PI_MASS, 2.4);
  for (uint i = 0; i < Nwaves; i++) {
    hint[i] = new TH2D(TString::Format("h2int%d", i + 1),
                       TString::Format("Intensity, %s", waves[i].title.c_str()),
                       NbinsX, rangeX.first, rangeX.second, NbinsY,
                       rangeY.first, rangeY.second);
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
        double valSq = 0.0;
        for (int M = -waves[w].J; M <= static_cast<int>(waves[w].J); M++) {
          uint minJS = (waves[w].S < waves[w].J) ? waves[w].S : waves[w].J;
          for (int lam = -minJS; lam <= static_cast<int>(minJS); lam++) {
            cd val(0., 0.);
            if (en > m23 + sqrt(m1sq)) {
              // calculate phase space
              val = MAscoli::getProjectionJMSlam(
                  waves[w].J, M, lam, s1, waves[w].S, R, s, t, mtRsq, stot,
                  mAsq, mBsq, mDsq, POW2(PI_MASS));
              // multiply to production clebsch, Blatt-Weisskopf
              val *= isobarT[waves[w].S](s1);
            }
            valSq += std::norm(val);
          }
        }
        // set bin content
        hint[w]->SetBinContent(bx + 1, by + 1, valSq * jac);
      }
    }
  }
  TFile *fout = new TFile(fout_name, "RECREATE");
  for (uint w = 0; w < Nwaves; w++)
    hint[w]->Write();
  std::cout << "File " << fout->GetName() << " have been completed!\n";
  fout->Close();
  return 0;
}

void fill_wavepull(const char *wave_fname, std::vector<wave> *waves) {
  uint index = 1;
  // only positive reflectivity
  for (uint J = 0; J < 5; J++) {
    for (uint S = 0; S < 3; S++) {
      wave w;
      w.index = index++;
      w.J = J;
      w.S = S;
      w.title = TString::Format(
          "sum{#lambda,M} |1-(J^{PC}=%d^{?+})M+(S=%d)pi|^{2}", J, S);
      waves->push_back(w);
    }
  }
}
