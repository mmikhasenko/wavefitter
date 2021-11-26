// Copyright [2017] Misha Mikhasenko

#include <utility>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TH2D.h"


uint plot_deck_2D(const char *fin_name = "/tmp/result.invert_matrix_non_symm.root") {
  TFile *fin = TFile::Open(fin_name);
  if(!fin) return 1;
  const uint nWaves = 88;
  TH1D *h[nWaves];
  for (uint w = 0; w < nWaves; w++) {
    h[w] = static_cast<TH1D*>(gDirectory->Get(TString::Format("h%d", w+1)));
    if (!h[w]) { std::cout << "No hist!\n"; return 1; }
  }
  std::vector<pair<uint, double> > w_max(nWaves-1);  // minus FLAT
  for (uint w = 1; w < nWaves; w++) {
    // w_max[w-1] = std::make_pair(w, h[w]->GetBinContent(h[w]->GetMaximumBin()));
    w_max[w-1] = std::make_pair(w, h[w]->Integral());
  }

  std::sort(w_max.begin(), w_max.end(),
            [](std::pair<double, double> p1, std::pair<double, double> p2)->bool {
      return (p1.second > p2.second);
    });
  for (auto && i : w_max) std::cout << i.second << " " << h[i.first]->GetTitle() << "\n";

  const uint nBins = 100;
  TH2D *h2 = new TH2D("Combined", "Combined distrs for Deck", 100, 0.5, 2.5, nWaves-1, 0, nWaves-1);
  for (uint w = 0; w < nWaves-1; w++) {
    for (uint e = 0; e < nBins; e++) {
      h2->SetBinContent(e+1, w+1, h[w_max[w].first]->GetBinContent(e+1));
    }
    h2->GetYaxis()->SetBinLabel(w+1, h[w_max[w].first]->GetTitle());
  }

  h2->SetStats(kFALSE);
  // h2->GetYaxis()->SetLabelSize(0.02);
  h2->Draw("colz");
  return 0;
}
