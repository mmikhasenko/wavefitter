// Copyright [2016] Misha Mikhasenko
// Discription

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "draw_distribution_from_tree.C"

TCanvas *compare_deck_expansion_and_pure(const char *quantity,
                                         const char *fproj_name,
                                         const char *fin_templ,
                                         uint bin,
                                         bool SYMM = false);

TCanvas *compare_deck_expansion_and_pure(const char *quantity,
                                         const char *fproj_name,
                                         const char *fin_templ,
                                         uint bin,
                                         bool SYMM) {
  std::vector< std::pair<std::string, std::vector<uint> > > waves;
  { std::vector<uint> v = {}; waves.push_back(std::make_pair("Coherent sum", v)); }
  { std::vector<uint> v = { 2,  5,  8,  9, 10, 18, 19, 20, 21, 22,
                           25, 33, 34, 35, 36, 37, 38, 39, 43, 45,
                           46, 47, 48, 51, 52, 53, 56, 57, 62, 63,
                           64, 67, 78, 76, 77, 78, 80, 82, 83, 84, 85};
    waves.push_back(std::make_pair("#rho", v)); }
  { std::vector<uint> v = { 7, 13, 14, 15, 23, 24, 26, 27, 28, 29,
                           30, 31, 32, 44, 49, 50, 58, 59, 61, 65,
                           66, 71, 72, 73, 79, 81, 86, 87, 88};
    waves.push_back(std::make_pair("f_{2}", v)); }
  { std::vector<uint> v = {6};
    waves.push_back(std::make_pair("f_{0}(1500)", v)); }
  { std::vector<uint> v = {3, 11, 12, 40, 41, 54, 55, 60, 69, 70, 74, 75};
    waves.push_back(std::make_pair("#(pipi)_{S}", v)); }
  { std::vector<uint> v = {4, 16, 17, 42};
    waves.push_back(std::make_pair("f_{0}(980)", v)); }
  { std::vector<uint> v = {18, 19, 25, 38, 39, 51, 52, 53, 64, 68, 78};
    waves.push_back(std::make_pair("#rho_{3}", v)); }

  // all waves
  THStack *hs = new THStack("hs_dex", TString::Format("d #sigma / d (%s) distribution, sqrt(s)=%2.2f GeV;%s}",
                                                      quantity, 0.5+(2.5-0.5)/100*bin, quantity));
  for (uint i = 0; i < waves.size(); i++) {
    TString hname = TString::Format("h%d", i);
    draw_distribution_from_tree(TString::Format("%s>>%s(300)", quantity, hname.Data()),
                                fproj_name,
                                fin_templ,
                                bin,
                                waves[i].second,
                                SYMM);
    TH1D *h = static_cast<TH1D*>(gROOT->FindObject(hname));
    if (!h) return 0;
    h->SetStats(kFALSE);
    h->SetTitle(waves[i].first.c_str());

    // if (i == 0) { h->SetLineColor(kBlack); h->SetFillStyle(0); }
    hs->Add(h);
  }

  TFile *fin = new TFile(TString::Format(fin_templ, bin-1));
  if (!fin) return 0;
  TTree *tin; gDirectory->GetObject("angles", tin);
  if (!tin) return 0;
  if (SYMM) {
    tin->Draw(TString::Format("%s>>h2(300)", quantity), "(decklike1_real+decklike3_real)**2 / 2 + (decklike1_imag+decklike3_imag)**2 / 2");
  } else {
    tin->Draw(TString::Format("%s>>h2(300)", quantity), "decklike1_real**2+decklike1_imag**2");
  }
  TH1D *h2 = static_cast<TH1D*>(gROOT->FindObject("h2"));
  if (!h2) return 0;

  h2->SetLineColor(kRed);
  h2->SetFillStyle(0);
  hs->Add(h2);
  h2->SetTitle("Analytic");

  TCanvas *can = new TCanvas("c1", "title");
  hs->Draw("hist pfc nostack");
  // set opacity

  // for (const auto && obj : *(hs->GetHists())) static_cast<TH1D*>(obj)->SetFillColorAlpha(static_cast<TH1D*>(obj)->GetFillColor(), 0.3);
    
  can->BuildLegend(0.55, 0.65, 0.9, 0.9);
  
  return can;
}
