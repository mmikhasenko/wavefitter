// Copyright [2016] Misha Mikhasenko
// Discription:
//    It calculates Deck projections based on event by event phase space MonteCarlo

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
#include "TRandom.h"
#include "TText.h"

#include "MDeck.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "waves.h"
#include "mintegrate.h"
#include "M3bodyAngularBasis.h"

typedef struct {
  uint index;
  std::string title;
} wave;

int main(int argc, char *argv[]) {

  const char *fin_tmpl = "/mnt/data/compass/2008/phase_space_MC/_with_deck_and_PWs_%d_large1e6.root";  // _large1e6
  TFile *f = TFile::Open(TString::Format(fin_tmpl, 0));
  if (!f) {std::cout << "Error: no file 0!\n"; return 0;}
  // check and load
  double nWaves = 88;
  std::vector<wave> waves(nWaves);
  for (uint w = 0; w < nWaves; w++) {
    waves[w].index = w+1;
    TText *t;
    gDirectory->GetObject(TString::Format("t%d", w+1), t);
    if (!t) { std::cerr << "Text field t" << w+1 << " is not found. You gonna get empty title!\n"; }
    else waves[w].title = std::string(t->GetTitle());
  }

  TH1D *hreal[waves.size()], *himag[waves.size()];
  for (uint iw = 0; iw < waves.size(); iw++) {
    hreal[iw] = new TH1D(TString::Format("r%d", waves[iw].index),
                         waves[iw].title.c_str(),
                         100, 0.5, 2.5);
    himag[iw] = new TH1D(TString::Format("i%d", waves[iw].index),
                         waves[iw].title.c_str(),
                         100, 0.5, 2.5);
  }

  for (uint e = 0; e < 100; e++) {
    std::cout << "---> File #" << e << "\n";
    TString fin_name = TString::Format(fin_tmpl, e);

    // open file and check
    TFile *f = TFile::Open(fin_name);
    if (!f) {std::cout << "Error: no file" << std::endl; return 0;}
    TTree *tin = 0; gDirectory->GetObject("angles", tin);
    if (!tin) {std::cout << "Error: no tree" << std::endl; return 0;}

    double s;
    tin->SetBranchAddress("s", &s);


    // general
    double s0, t;
    tin->SetBranchAddress("s0", &s0);
    tin->SetBranchAddress("t", &t);

    // deck
    double decklike_real[2], decklike_imag[2];
    tin->SetBranchAddress("decklike1_real", &decklike_real[0]);
    tin->SetBranchAddress("decklike1_imag", &decklike_imag[0]);
    tin->SetBranchAddress("decklike3_real", &decklike_real[1]);
    tin->SetBranchAddress("decklike3_imag", &decklike_imag[1]);

    // add a branch for every wave
    uint nWaves = waves.size();
    double amp_real[nWaves][2], amp_imag[nWaves][2];
    for (uint w = 0; w < nWaves; w++) {
      tin->SetBranchAddress(TString::Format("amp%d_frame1_real", waves[w].index), &amp_real[w][0]);
      tin->SetBranchAddress(TString::Format("amp%d_frame1_imag", waves[w].index), &amp_imag[w][0]);
      tin->SetBranchAddress(TString::Format("amp%d_frame3_real", waves[w].index), &amp_real[w][1]);
      tin->SetBranchAddress(TString::Format("amp%d_frame3_imag", waves[w].index), &amp_imag[w][1]);
    }

    // Phhase space
    tin->GetEntry(0);  // to get s
    double phsp = integrate([s](double _s1)->double{
        return sqrt(LAMBDA(s, _s1, POW2(PI_MASS))*LAMBDA(_s1, POW2(PI_MASS), POW2(PI_MASS)))/_s1;
      }, 4*POW2(PI_MASS), POW2(sqrt(s)-PI_MASS)) / (2*M_PI*POW2(8*M_PI)*s);

    // create and clean integral variables
    cd integrals[nWaves];
    for (uint iw = 0; iw < waves.size(); iw++) integrals[iw] = 0.;

    // integration loop
    const int Nentries = tin->GetEntries();
    for (int i = 0; i < Nentries; i++) {
      if (i%1000000 == 0 && i != 0) std::cout << "Processing entry " << i << "\n";
      tin->GetEntry(i);

      cd deck[2];
      cd amp[nWaves][2];
      for (uint bose = 0; bose < 2; bose++) {
        for (uint w = 0; w < nWaves; w++) amp[w][bose] = cd(amp_real[w][bose], amp_imag[w][bose]);
        deck[bose] = cd(decklike_real[bose], decklike_imag[bose]);
      }
      for (uint iw = 0; iw < waves.size(); iw++)
        integrals[iw] +=
          // conj(amp[iw][0])*(deck[0]);
          conj(amp[iw][0]+amp[iw][1])*(deck[0]+deck[1])/2.;
    }
    f->Close();

    for (uint iw = 0; iw < waves.size(); iw++) {
      hreal[iw]->SetBinContent(e+1, real(integrals[iw])*phsp*POW2(4*M_PI)/Nentries * (8*M_PI));
      himag[iw]->SetBinContent(e+1, imag(integrals[iw])*phsp*POW2(4*M_PI)/Nentries * (8*M_PI));
    }
  }

  TFile fout("/tmp/waves.calculate_Deck_integrals_symm.root", "RECREATE");
  for (uint iw = 0; iw < waves.size(); iw++) {
    hreal[iw]->Write();
    himag[iw]->Write();
  }

  fout.Close();
  std::cout << "File " << fout.GetName() << " has been created\n";
  return 0;
}
