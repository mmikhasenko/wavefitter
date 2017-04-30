// Copyright [2016] Misha Mikhasenko
// Description:
//   The program calculates integrals numerically based on pregenerated MC data sample

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

  const char *fin_tmpl = "/mnt/data/compass/2008/phase_space_MC/_with_deck_and_PWs_%d_large1e6.root";
  TFile *f = TFile::Open(TString::Format(fin_tmpl, 0));
  if (!f) {std::cout << "Error: no file 0!\n"; return 0;}
  // check and load
  uint nWaves = 88;
  std::vector<wave> waves(nWaves);
  for (uint w = 0; w < nWaves; w++) {
    waves[w].index = w+1;
    TText *t;
    gDirectory->GetObject(TString::Format("t%d", w+1), t);
    if (!t) { std::cerr << "Text field t" << w+1 << " is not found. You gonna get empty title!\n"; }
    else waves[w].title = std::string(t->GetTitle());
  }

  TH1D *hdiag[waves.size()];
  TH1D *hinterf[waves.size()][waves.size()];
  for (uint iw = 0; iw < waves.size(); iw++) {
    for (uint jw = 0; jw < waves.size(); jw++) {
      if (jw > iw)
        hinterf[iw][jw] = new TH1D(TString::Format("h1%03d%03d", waves[iw].index, waves[jw].index),
                                   TString::Format("Re[ %s* %s ]", waves[iw].title.c_str(), waves[jw].title.c_str()),
                                   100, 0.5, 2.5);
      if (jw < iw)
        hinterf[iw][jw] = new TH1D(TString::Format("h2%03d%03d", waves[jw].index, waves[iw].index),
                                   TString::Format("Im[ %s* %s ]", waves[jw].title.c_str(), waves[iw].title.c_str()),
                                   100, 0.5, 2.5);
      if (jw == iw)
        hdiag[iw] = new TH1D(TString::Format("h%d", waves[iw].index),
                            waves[iw].title.c_str(),
                            100, 0.5, 2.5);
    }
  }

  for (uint e = 0; e < 100; e++) {
    std::cout << "---> File #" << e << "\n";
    TString fin_name = TString::Format(fin_tmpl, e);  //

    // open file and check
    TFile *f = TFile::Open(fin_name);
    if (!f) {std::cout << "Error: no file" << std::endl; return 0;}
    TTree *tin = 0; gDirectory->GetObject("angles", tin);
    if (!tin) {std::cout << "Error: no tree" << std::endl; return 0;}

    // add a branch for every wave
    double amp_real[nWaves][2], amp_imag[nWaves][2];
    for (uint w = 0; w < nWaves; w++) {
      tin->SetBranchAddress(TString::Format("amp%d_frame1_real", waves[w].index), &amp_real[w][0]);
      tin->SetBranchAddress(TString::Format("amp%d_frame1_imag", waves[w].index), &amp_imag[w][0]);
      tin->SetBranchAddress(TString::Format("amp%d_frame3_real", waves[w].index), &amp_real[w][1]);
      tin->SetBranchAddress(TString::Format("amp%d_frame3_imag", waves[w].index), &amp_imag[w][1]);
    }

    // Phhase space
    double s; tin->SetBranchAddress("s", &s); tin->GetEntry(0);  // to use it at the next line
    double phsp = integrate([s](double _s1)->double{
        return sqrt(LAMBDA(s, _s1, POW2(PI_MASS))*LAMBDA(_s1, POW2(PI_MASS), POW2(PI_MASS)))/_s1;
      }, 4*POW2(PI_MASS), POW2(sqrt(s)-PI_MASS)) / (2*M_PI*POW2(8*M_PI)*s);

    // select a set of non-zero waves

    // create and clean integral variables
    cd integrals[nWaves][nWaves];
    for (uint iw = 0; iw < nWaves; iw++)
      for (uint jw = 0; jw < nWaves; jw++)
        integrals[iw][jw] = 0.;

    // integration loop
    const int Nentries = tin->GetEntries();
    for (int i = 0; i < Nentries; i++) {
      if (i%1000000 == 0 && i != 0) std::cout << "Processing entry " << i << "\n";
      tin->GetEntry(i);

      cd amp[nWaves][2];
      for (uint w = 0; w < nWaves; w++)
        for (uint bose = 0; bose < 2; bose++)
          amp[w][bose] = cd(amp_real[w][bose], amp_imag[w][bose]);

      for (uint iw = 0; iw < nWaves; iw++)
        for (uint jw = iw; jw < nWaves; jw++)
          integrals[iw][jw] +=
            // conj(amp[iw][0])*(amp[jw][0]);  // non-symmetrized
            conj(amp[iw][0]+amp[iw][1])*(amp[jw][0]+amp[jw][1])/2.;  // symmetrized
    }
    f->Close();

    for (uint iw = 0; iw < waves.size(); iw++) {
      for (uint jw = 0; jw < waves.size(); jw++) {
        if (jw > iw)
          hinterf[iw][jw]->SetBinContent(e+1, real(integrals[iw][jw])*phsp*POW2(4*M_PI)/Nentries * (8*M_PI));
        if (jw < iw)
          hinterf[iw][jw]->SetBinContent(e+1, imag(integrals[jw][iw])*phsp*POW2(4*M_PI)/Nentries * (8*M_PI));
        if (jw == iw) {
          hdiag[iw]->SetBinContent(e+1,  abs(integrals[iw][jw])*phsp*POW2(4*M_PI)/Nentries * (8*M_PI));
        }
      }
    }
  }

  TFile fout("/tmp/waves.calculate_phase_space_symm.root", "RECREATE");
  for (uint iw = 0; iw < waves.size(); iw++)
    for (uint jw = 0; jw < waves.size(); jw++)
      if (jw == iw) hdiag[iw]->Write();
      else hinterf[iw][jw]->Write();

  fout.Close();
  std::cout << "File " << fout.GetName() << " has been created\n";
  return 0;
}

/***************************************************************************************/
/**** script to plot matrix ************************************************************/
/***************************************************************************************/
/*
TCanvas c1
c1.Divide(4,4)
for(int i=0; i<4; i++) {
for (int j=0; j<4; j++) {
c1.cd(1+i*4+j);
if(i==j) gDirectory->Get(TString::Format("h%d",i+1))->Draw();
if(i<j) gDirectory->Get(TString::Format("h1%03d%03d",i+1,j+1))->Draw();
if(i>j) gDirectory->Get(TString::Format("h2%03d%03d",j+1,i+1))->Draw();
}
}
*/
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/


/***************************************************************************************/
/**** script to plot in the matrix form and colors *************************************/
/***************************************************************************************/
/*
TMatrixD m(88,88);
for(int i=0; i<88;i++) {
  for (int j=0; j<88; j++) {
    TH1D *h = (TH1D*)gDirectory->Get((i==j) ? TString::Format("h%d",i+1) : ((i<j) ? TString::Format("h1%03d%03d",i+1,j+1) : TString::Format("h1%03d%03d",j+1,i+1)));
    m(i,j) = h->GetBinContent(81);
  }
}
TCanvas *c1 = new TCanvas("c1", "title", 1000, 1000);
m.Draw("colz");
TH2D *h = (TH2D*)gROOT->FindObject("TMatrixDBase");
if (h) {
  for(int i=0; i<88;i++) {
    TH1D *hdiag = (TH1D*)gROOT->FindObject(TString::Format("h%d", i+1));
    if (!hdiag) break;
    h->GetXaxis()->SetBinLabel(i+1, hdiag->GetTitle());
    h->GetYaxis()->SetBinLabel(i+1, hdiag->GetTitle());
  }
}
h->SetStats(kFALSE);
h->GetXaxis()->SetLabelSize(0.02);
h->GetYaxis()->SetLabelSize(0.02);
h->GetZaxis()->SetRangeUser(0.01, 1.0);
h->GetXaxis()->SetTickSize(0);
h->GetYaxis()->SetTickSize(0);
h->SetTitle("PW integral matrix");
h->Draw("colz");
*/
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
