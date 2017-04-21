// Copyright [2016] Misha Mikhasenko
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "MDeck.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "waves.h"
#include "mintegrate.h"
#include "M3bodyAngularBasis.h"

typedef struct {
  uint index;
  uint J;
  uint M;
  bool neg_refl;
  int S;
  uint L;
  double threshold;
  std::string title;
} wave;

// wave make_wave(uint _J, uint _M, uint _S, uint _L, bool _neg_refl = true, double _threshold = 0.5,
//                const std::string &_title = "") {
//   wave n; n.J = _J; n.M = _M; n.S = _S; n.L = _L;
//   n.neg_refl = _neg_refl; n.title = _title; n.threshold = _threshold;
//   return n;
// }

void fill_wavepull(const char* wave_fname, std::vector<wave> *waves);

int main(int argc, char *argv[]) {


  std::vector<wave> waves;
  fill_wavepull("/localhome/mikhasenko/results/pwa_3pi/wavelist_formated.txt", &waves);
  std::cout << "waves.size() = " << waves.size() << "\n";
  for (auto & w : waves)
    std::cout << w.index << ": " << w.title << " " << w.J << " " << w.M << " " << w.S << " " << w.L << "\n";


  TH1D *hsymm[waves.size()];
  TH1D *hnsym[waves.size()];
  for (uint w = 0; w < waves.size(); w++) {
    hsymm[w] = new TH1D(TString::Format("h%d", waves[w].index),
                     waves[w].title.c_str(),
                     100, 0.5, 2.5);
    hnsym[w] = new TH1D(TString::Format("h%d_nsym", waves[w].index),
                     waves[w].title.c_str(),
                     100, 0.5, 2.5);
  }

  // for out mode;
  MIsobar rho(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.); rho.setIntU();
  MIsobar  f2(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.); f2.setIntU();
  MIsobar rho3(1.69, 0.16, PI_MASS, PI_MASS, 1, 5.); rho3.setIntU();
  MIsobarPiPiS pipiS; pipiS.setIntU();
  MIsobar *iso[] = {&pipiS, &rho, &f2, &rho3};
  MIsobar f980(0.99, 0.04, PI_MASS, PI_MASS, 0); f980.setIntU();
  MIsobar f1500(1.504, 0.11, PI_MASS, PI_MASS, 0); f1500.setIntU();
  MIsobar *iso_scalars[] = {&pipiS, &f980, &f1500};

  for (uint e = 0; e < 100; e++) {
    std::cout << "---> File #" << e << "\n";
    TString fin_name = TString::Format("/mnt/data/compass/2008/phase_space_MC/%d.root", e);

    // open file and check
    TFile *f = TFile::Open(fin_name);
    if (!f) {std::cout << "Error: no file" << std::endl; return 0;}
    TTree *tin = 0; gDirectory->GetObject("angles", tin);
    if (!tin) {std::cout << "Error: no tree" << std::endl; return 0;}

    double s;
    tin->SetBranchAddress("s", &s);

    // (23)-frame
    double s1, costheta1, phi1, costheta23, phi23;
    tin->SetBranchAddress("s1", &s1);
    tin->SetBranchAddress("costheta1", &costheta1);
    tin->SetBranchAddress("phi1", &phi1);
    tin->SetBranchAddress("costheta23", &costheta23);
    tin->SetBranchAddress("phi23", &phi23);

    // (12)-frame
    double s3, costheta3, phi3, costheta12, phi12;
    tin->SetBranchAddress("s3", &s3);
    tin->SetBranchAddress("costheta3", &costheta3);
    tin->SetBranchAddress("phi3", &phi3);
    tin->SetBranchAddress("costheta12", &costheta12);
    tin->SetBranchAddress("phi12", &phi12);

    // Phhase space
    tin->GetEntry(0);
    double phsp = integrate([s](double _s1)->double{
        return sqrt(LAMBDA(s, _s1, POW2(PI_MASS))*LAMBDA(_s1, POW2(PI_MASS), POW2(PI_MASS)))/_s1;
      }, 4*POW2(PI_MASS), POW2(sqrt(s)-PI_MASS)) / (2*M_PI*POW2(8*M_PI)*s);
    // with isobar
    // double quasi_two_body_phsp = integrate([&, s](double s1)->double{
    //     return sqrt(LAMBDA(s, s1, POW2(PI_MASS))*LAMBDA(s1, POW2(PI_MASS), POW2(PI_MASS)))/s1 *
    //       norm(iso[1]->ToneVertex(s1));
    //   }, 4*POW2(PI_MASS), POW2(sqrt(s)-PI_MASS)) / (2*M_PI*POW2(8*M_PI)*s);

    // select a set of non-zero waves
    uint nWaves = waves.size();
    // create and clean integral variables
    double integral[nWaves], integral_symm[nWaves];
    for (uint w = 0; w < nWaves; w++) { integral[w] = 0; integral_symm[w] = 0; }

    const int Nentries = tin->GetEntries();
    for (int i = 0; i < Nentries; i++) {
      if (i%1000000 == 0 && i != 0) std::cout << "Processing entry " << i << "\n";
      tin->GetEntry(i);

      cd amp[nWaves][2];
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
          amp[w][bose] = Math::ZJMLS_refl(waves[w].J, waves[w].M, waves[w].L, waves[w].S > 0 ? waves[w].S : 0,
                                          waves[w].neg_refl,
                                          thetaI, phiI, theta, phi) * iso_shape;
        }
        /*************************************/
      }
      for (uint w = 0; w < nWaves; w++) {
        integral[w] += norm(amp[w][0]);
        integral_symm[w] += norm(amp[w][0]+amp[w][1])/2.;
      }
    }
    f->Close();

    for (uint w = 0; w < nWaves; w++) {
      hsymm[w]->SetBinContent(e+1, integral_symm[w]*phsp*POW2(4*M_PI)/Nentries * (8*M_PI));
      hnsym[w]->SetBinContent(e+1, integral[w]*phsp*POW2(4*M_PI)/Nentries * (8*M_PI));
    }
  }

  TCanvas c1("c1", "title", 1000, 1000);
  c1.DivideSquare(waves.size());
  for (uint w = 0; w < waves.size(); w++) {
    c1.cd(w+1);
    hsymm[w]->Draw();
    hnsym[w]->SetLineColor(kRed); hnsym[w]->Draw("same");
  }
  c1.Print("/tmp/ph.sp.PhoPiS.pdf");
  c1.SaveAs("/tmp/ph.sp.PhoPiS.pdf");

  TFile fout("/tmp/waves.calculate_phase_space.root", "RECREATE");
  for (uint w = 0; w < waves.size(); w++) { hsymm[w]->Write(); hnsym[w]->Write(); }
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
      n.neg_refl = true;
      waves->push_back(n);
      continue;
    }
    iss >> n.J;
    iss >> n.M;
    std::string epsilon;
    iss >> epsilon;
    if (epsilon == "-") n.neg_refl = true; else n.neg_refl = false;
    iss >> n.S;
    iss >> n.L;
    waves->push_back(n);
  }  
}
