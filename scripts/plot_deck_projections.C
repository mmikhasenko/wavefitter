// Copyright [2017] Misha Mikhasenko

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TRandom.h"


uint plot_deck_projections(uint nPads = 4);
uint plot_deck_projections(uint nPads) {
  TFile *fin1 = new TFile("/mnt/data/compass/2008/Deck.88waves.linal.misha/inverted.Deck.large.root");
  if (!fin1) { std::cerr << "Error with file1!\n"; return 1; }
  TFile *fin2 = new TFile("/mnt/data/compass/2008/Deck.88waves.linal.misha/inverted.Deck.thresholded.large.root");
  if (!fin2) { std::cerr << "Error with file2!\n"; return 1; }

  TCanvas *c1 = new TCanvas("c1");
  c1->DivideSquare(nPads);
  for (uint p = 0; p < nPads; p++) {
    const uint index = 1+gRandom->Integer(81);

    fin1->cd(); TH1D *h1 = static_cast<TH1D*>(gDirectory->Get(TString::Format("h%d", index)));
    if (!h1) { std::cerr << "Error with hist h" << index << " in file 1!\n"; return 1; }
    // h1->SetStats(kFALSE);

    fin2->cd(); TH1D *h2 = static_cast<TH1D*>(gDirectory->Get(TString::Format("h%d", index)));
    if (!h2) { std::cerr << "Error with hist h" << index << " in file 2!\n"; return 1; }
    // h2->SetStats(kFALSE);

    c1->cd(p+1);
    h2->Draw();
    h1->SetLineColor(kGray); h1->SetFillColor(0); h1->Draw("same");
    // Alpha(kBlue, 0.
  }
  return 0;
}
