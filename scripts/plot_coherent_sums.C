// Copyright [2017] Misha Mikhasenko
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>

#include "TH1D.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"

typedef std::complex<double> cd;

TH1D *plot_coherent_sums(const std::vector<uint> &_indeces,
                         const char *ph_sp, const char *inversions) {
  std::vector<uint> indeces(_indeces);
  std::sort(indeces.begin(), indeces.end());

  const uint nBins = 100;
  const uint nInd = indeces.size();
  // phase space
  TFile *fin_ph_sp = new TFile(ph_sp); if (!fin_ph_sp) { std::cerr << "No file 1\n"; return 0; }
  TFile *fin_invs  = new TFile(inversions); if (!fin_invs) { std::cerr << "No file 2\n"; return 0; }

  TString title("Coherent sum");  // for (auto & i : indeces) title += TString::Format("%d,", i); title += "]";
  TH1D *hint = new TH1D("hint", title, nBins, 0.5, 2.5);

  for (uint b = 0; b < nBins; b++) {
    double M = 0.51+0.02*b;
    fin_ph_sp->cd();
    cd amp[88][88];
    for (uint i = 0; i < nInd; i++) {
      for (uint j = 0; j < nInd; j++) {
        TH1D *h;
        if (i == j) {
          gDirectory->GetObject(TString::Format("h%d", indeces[i]), h); if (!h) {std::cerr << "No hist 1\n"; return 0;}
          amp[i][j] = cd(h->GetBinContent(b+1), 0);
        } else if (i < j) {
          gDirectory->GetObject(TString::Format("h1%03d%03d", indeces[i], indeces[j]), h); if (!h) {std::cerr << "No hist " << i << j << "\n"; return 0;}
          amp[i][j] += h->GetBinContent(b+1) - real(amp[i][j]);  // upper triangle
          amp[j][i] += h->GetBinContent(b+1) - real(amp[j][i]);  // conjugated
        } else {
          gDirectory->GetObject(TString::Format("h2%03d%03d", indeces[j], indeces[i]), h); if (!h) {std::cerr << "No hist 21\n"; return 0;}
          amp[j][i] += cd(0, h->GetBinContent(b+1) - imag(amp[j][i]));  // upper triangle
          amp[i][j] += cd(0, -h->GetBinContent(b+1) - imag(amp[i][j]));  // conjugated
        }
      }  // j
    }  // i
    // std::cout << "Integrals extruction is done\n";

    // Deck
    cd proj[88];
    fin_invs->cd();
    for (uint i = 0; i < nInd; i++) {
      TH1D *h1; gDirectory->GetObject(TString::Format("hr%d", indeces[i]), h1); if (!h1) {std::cerr << "No hist r1\n"; return 0;}
      TH1D *h2; gDirectory->GetObject(TString::Format("hi%d", indeces[i]), h2); if (!h2) {std::cerr << "No hist i1\n"; return 0;}
      proj[i] = cd(h1->GetBinContent(b+1), h2->GetBinContent(b+1));
    }
    // std::cout << "Projections extraction is done\n";

    // result
    cd sum = 0;
    for (uint i = 0; i < nInd; i++) {
      for (uint j = 0; j < nInd; j++) {
        sum += conj(proj[i])*proj[j] * M *
        (amp[i][j] * 1./(8*M_PI)) * // because of the normalization of integrals
        1./((4*M_PI)*(4*M_PI)); // additional factors from the phase space
        //
      }
    }
    hint->SetBinContent(b+1, real(sum));
    if (imag(sum) > 1e-6) {std::cerr << "Something is wrong!\n"; return 0;}
    // std::cout << "Some is done\n";
  }
  return hint;
}
