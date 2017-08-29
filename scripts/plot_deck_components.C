// Copyright [2017] Misha Mikhasenko
// Discription:
//   the script plots projected Deck-like signal

#include <iostream>
#include <vector>
#include <string>
#include <utility>

#include "TROOT.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "THStack.h"

#include "plot_coherent_sums.C"

int plot_deck_components(const char *fin_proj_name, const char *fin_phsp_name,
                         const char *save_name = 0);
int plot_deck_components(const char *fin_proj_name, const char *fin_phsp_name,
                         const char *save_name) {

  std::vector< std::pair<std::string, std::vector<uint> > > waves;
  { std::vector<uint> v = {}; waves.push_back(std::make_pair("Coherent sum;M_{3#pi}", v)); }
  { std::vector<uint> v = {2, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}; waves.push_back(std::make_pair("1^{++}", v)); }
  { std::vector<uint> v = {6, 3, 4, 5, 7}; waves.push_back(std::make_pair("0^{-+}", v)); }
  { std::vector<uint> v = {26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 88}; waves.push_back(std::make_pair("2^{-+}", v)); }
  { std::vector<uint> v = {21, 22, 23, 24, 25, 85, 86, 87}; waves.push_back(std::make_pair("2^{++}", v)); }
  { std::vector<uint> v = {20, 82, 83}; waves.push_back(std::make_pair("1^{-+}", v)); }
  // { std::vector<uint> v = {43, 44}; waves.push_back(std::make_pair("3^{-+}", v)); }
  // { std::vector<uint> v = {45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55}; waves.push_back(std::make_pair("3^{++}", v)); }
  // { std::vector<uint> v = {56, 57, 58, 59, 60, 61}; waves.push_back(std::make_pair("4^{-+}", v)); }
  // { std::vector<uint> v = {62, 63, 64, 65, 66}; waves.push_back(std::make_pair("4^{++}", v)); }
  // { std::vector<uint> v = {67, 68, 69, 70, 71, 72, 73}; waves.push_back(std::make_pair("5^{++}", v)); }
  // { std::vector<uint> v = {74, 75, 76, 77, 78, 79}; waves.push_back(std::make_pair("6^{-+}", v)); }
  // { std::vector<uint> v = {80, 81}; waves.push_back(std::make_pair("6^{++}", v)); }
  // { std::vector<uint> v = {82, 83, 84, 85, 86, 87, 88}; waves.push_back(std::make_pair("#epsilon = (-)", v)); }

  waves[0].second.resize(88);
  for (uint i = 0; i < waves[0].second.size(); i++) waves[0].second[i] = i+1;

  TCanvas *c1 = new TCanvas("c1", "title", 1000, 1000);

  THStack *hs = new THStack("hs", "Intensities for the coherent sum of projections;M_{3#pi}");
  for (auto & i : waves) {
    TH1D *h =
      plot_coherent_sums(i.second,
                         fin_phsp_name,
                         fin_proj_name);
    h->SetTitle(i.first.c_str());
    h->SetStats(kFALSE);
    hs->Add(h);
  }
  hs->Draw("pfc nostack");

  c1->BuildLegend(0.7, 0.6, 0.9, 0.9);
  if (save_name != 0) c1->SaveAs(TString::Format("%s_all.pdf", save_name));

  TFile *fin = new TFile(fin_proj_name);
  for (auto & i : waves)
  std::sort(i.second.begin(), i.second.end(), [&](uint i1, uint i2)->bool{
      TH1D *h1 = 0; gDirectory->GetObject(TString::Format("h%d", i1), h1); if (!h1) return false;
      TH1D *h2 = 0; gDirectory->GetObject(TString::Format("h%d", i2), h2); if (!h2) return false;
      if (h1->GetBinContent(h1->GetMaximumBin()) > h2->GetBinContent(h2->GetMaximumBin())) return true;
      return false;
    });
  for (auto & i : waves) {
    THStack *hsi = new THStack(TString::Format("hs%s", i.first.c_str()),
                              TString::Format("Intensities for %s;M_{3#pi}", i.first.c_str()));
    fin->cd();
    for (auto & j : i.second) {
      TH1D *h = 0; gDirectory->GetObject(TString::Format("h%d", j), h);
      if (!h) return -2;
      if (j == 6) h->SetFillColorAlpha(kYellow, 0.1);
      hsi->Add(h);
    }
    TCanvas *can = new TCanvas(TString::Format("hs%s", i.first.c_str()),
                               i.first.c_str(), 1000, 1000);
    TH1D *h =
      plot_coherent_sums(i.second,
                         fin_phsp_name,
                         fin_proj_name);
    h->SetLineColor(kRed);
    h->SetFillStyle(0);
    hsi->Add(h);
    hsi->Draw("pfc nostack");
    // h->Draw("same");
    can->BuildLegend(0.55, 0.6, 0.9, 0.9);
    if (save_name != 0) can->SaveAs(TString::Format("%s_%s.pdf", save_name, i.first.c_str()));
  }

  return 0;
}
