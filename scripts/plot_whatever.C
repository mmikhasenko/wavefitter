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
#include "THStack.h"
#include "TPad.h"

// header
TCanvas *compare_whatever(std::vector<uint> set,
                          const char *file1_name,
                          const char *file2_name,
                          const char *tmpl);
TCanvas *compare_whatever(uint nPads,
                          const char *file1_name,
                          const char *file2_name,
                          uint MAX, const char *tmpl);

TCanvas *plot_whatever(uint nPads = 4,
                       const char *file_name = "/tmp/result.invert_matrix.root",
                       uint MAX = 88,
                       const char *tmpl = "h%d");

// functions
TCanvas *plot_whatever(uint nPads, const char *file_name, uint MAX, const char *tmpl) {
  TFile *fin2 = new TFile(file_name);
  if (!fin2) { std::cerr << "Error with file2!\n"; return 0; }

  TCanvas *c1 = new TCanvas("c1");
  c1->DivideSquare(nPads);
  for (uint p = 0; p < nPads; p++) {
    const uint index = 1+gRandom->Integer(MAX);
    fin2->cd(); TH1D *hh = static_cast<TH1D*>(gDirectory->Get(TString::Format(tmpl, index)));
    if (!hh) { std::cerr << "Error with hist h" << index << "!\n"; return 0; }
    hh->SetStats(kFALSE);

    c1->cd(p+1);
    hh->Draw();
  }
  return c1;
}

TCanvas *compare_whatever(uint nPads,
                      const char *file1_name,
                      const char *file2_name,
                      uint MAX, const char *tmpl) {
  std::vector<uint> set(nPads);
  for (auto & i : set) i = 1+gRandom->Integer(MAX);
  return compare_whatever(set, file1_name, file2_name, tmpl);
}

TCanvas *compare_whatever(std::vector<uint> set,
                      const char *file1_name,
                      const char *file2_name,
                      const char *tmpl) {
  TFile *fin1 = new TFile(file1_name);
  if (!fin1) { std::cerr << "Error with file1!\n"; return 0; }
  TFile *fin2 = new TFile(file2_name);
  if (!fin2) { std::cerr << "Error with file2!\n"; return 0; }

  TCanvas *c1 = new TCanvas("c1");
  uint nPads = set.size();
  c1->DivideSquare(nPads);
  for (uint p = 0; p < nPads; p++) {
    c1->cd(p+1);
    gPad->IncrementPaletteColor(2, "pfc");

    uint index = set[p];
    fin1->cd(); TH1D *h1 = static_cast<TH1D*>(gDirectory->Get(TString::Format(tmpl, index)));
    if (!h1) { std::cerr << "Error with hist1 h" << index << "!\n"; return 0; }
    THStack *hs = new THStack("hs", h1->GetTitle());
    h1->SetFillColorAlpha(gPad->NextPaletteColor(), 0.5);
    h1->SetStats(kFALSE); hs->Add(h1);
    fin2->cd(); TH1D *h2 = static_cast<TH1D*>(gDirectory->Get(TString::Format(tmpl, index)));
    if (!h2) { std::cerr << "Error with hist2 h" << index << "!\n"; return 0; }
    h2->SetFillColorAlpha(gPad->NextPaletteColor(), 0.5);
    h2->SetStats(kFALSE); hs->Add(h2);

    // h2->SetFillStyle(0); h2->SetLineColor(kMagenta);
    hs->Draw("nostack");
    
  }
  return c1;
}
