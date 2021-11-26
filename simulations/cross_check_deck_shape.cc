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

int main(int argc, char *argv[]) {

  TH1D *h = new TH1D("h1", "h1", 100, 0.5, 2.5);

  const uint Nbins = 100;
  for (uint e = 0; e < Nbins; e++) {
    std::cout << "--- Process " << e << " file \n";
    double m3pi = 0.5 + (2.5-0.5)/Nbins*e;
    double s = POW2(m3pi);
    // Phhase space
    double phsp = integrate([s](double _s1)->double{
        return sqrt(LAMBDA(s, _s1, POW2(PI_MASS))*LAMBDA(_s1, POW2(PI_MASS), POW2(PI_MASS)))/_s1;
      }, 4*POW2(PI_MASS), POW2(sqrt(s)-PI_MASS)) / (2*M_PI*POW2(8*M_PI)*s);

    std::cout << "sqrt(s) = " << m3pi << ", phase space = " << phsp << "\n";
    TFile *file = TFile::Open(TString::Format("/mnt/data/compass/2008/phase_space_MC/direct_%d_large1e6.root", e));
    TTree *tin; gDirectory->GetObject("angles", tin); if (!tin) return 1;
    tin->Draw("sqrt(s)>>+h1", "(pow(decklike1_real+decklike3_real,2)+pow(decklike1_imag+decklike3_imag,2))");
    // double decklike1_real, decklike3_real, decklike1_imag, decklike3_imag;
    // tin->SetBranchAddress("decklike1_real", &decklike1_real);
    // tin->SetBranchAddress("decklike1_imag", &decklike1_imag);
    // tin->SetBranchAddress("decklike3_real", &decklike3_real);
    // tin->SetBranchAddress("decklike3_imag", &decklike3_imag);
    // double dw = 0;
    // const uint nEvents = tin->GetEntries();
    // for (uint j = 0; j < nEvents; j++) {
    //   tin->GetEntry(j);
    //   double dr = decklike1_real + decklike3_real;
    //   double di = decklike1_imag + decklike3_imag;
    //   dw += dr*dr+di*di;
    // }
    h->SetBinContent(e+1, phsp*h->GetBinContent(e+1)/tin->GetEntries());
    std::cout << h->GetBinContent(e+1) << "\n";
    file->Close();
  }
  TCanvas *c1 = new TCanvas("c1", "title", 1000, 1000);
  h->Draw();
  c1->Print("/tmp/c1.pdf");

  return 0;
}
