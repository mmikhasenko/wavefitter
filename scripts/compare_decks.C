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


uint compare_decks(uint nPads = 4);
uint compare_decks(uint nPads) {
  TFile *fin1 = new TFile("/mnt/data/compass/2008/Deck.88waves.linal.misha/waves.calculate_Deck_integrals.root");
  TFile *fin2 = new TFile("/mnt/data/compass/2008/Deck.88waves.linal.misha/waves.calculate_Deck_integrals_quicker.root");
  if (!fin1) { std::cerr << "Error with file1!\n"; return 1; }
  if (!fin2) { std::cerr << "Error with file2!\n"; return 1; }

  TCanvas *c1 = new TCanvas("c1");
  c1->DivideSquare(nPads);
  for (uint p = 0; p < nPads; p++) {
    const uint index = 1+gRandom->Integer(87);
    fin1->cd(); TH1D *h1; gDirectory->GetObject(TString::Format("r%d", index), h1);
    if (!h1) { std::cerr << "Error with hist1!\n"; return 1; }
    fin2->cd(); TH1D *h2; gDirectory->GetObject(TString::Format("r%d", index), h2);
    if (!h2) { std::cerr << "Error with hist2!\n"; return 1; }

    h1->SetStats(kFALSE);
    h2->SetStats(kFALSE);

//     h2->SetFillStyle(0);
    h2->SetLineColor(kGray);

    c1->cd(p+1);
    h1->Draw();
    h2->Draw("same");
  }
  return 0;
}
