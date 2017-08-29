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

int plot_onepp_fnot_sector(const char *fin_proj_name, const char *fin_phsp_name) {

  std::vector<uint> onepp = {2, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

  TCanvas *c1 = new TCanvas("c1", "title", 1000, 1000);

  THStack *hs = new THStack("hs", "Intensities of the partial waves in 1^{++} sector;M_{3#pi}");
  
  // gStyle->SetHatchesLineWidth(2);
  // all
  TH1D *hall =
    plot_coherent_sums(onepp,
                       fin_phsp_name,
                       fin_proj_name);
  hall->SetTitle("1^{++};M_{3#pi}");
  hall->SetStats(kFALSE);
  // hall->SetFillStyle(3019);
  hs->Add(hall);
  // [0++]
  TH1D *hfnots =
    plot_coherent_sums({11, 16},
                       fin_phsp_name,
                       fin_proj_name);
  hfnots->SetTitle("1^{++}0^{+} [#pi#pi]_{0^{++}}");
  hfnots->SetStats(kFALSE);
  // hfnots->SetFillStyle(3004);
  hs->Add(hfnots);
  // sigma(600)
  TH1D *h11 =
    plot_coherent_sums({11},
                       fin_phsp_name,
                       fin_proj_name);
  h11->SetTitle("1^{++}0^{+} (#pi#pi)_{S}");
  h11->SetStats(kFALSE);
  // h11->SetFillStyle(3005);
  hs->Add(h11);
  // f0(980)
  TH1D *h16 =
    plot_coherent_sums({16},
                       fin_phsp_name,
                       fin_proj_name);
  h16->SetTitle("1^{++}0^{+} f_{0}(980);M_{3#pi}");
  h16->SetStats(kFALSE);
  h16->SetFillStyle(3016);
  hs->Add(h16);

  hs->Draw("pfc nostack");
  for (const TObject *h : *(hs->GetHists())) ((TH1D*)h)->SetFillStyle(3023);
    // ((TH1D*)h)->SetFillColorAlpha(((TH1D*)h)->GetFillColor(), 0.35);
  c1->Update();
  
  THStack *hs2 = (THStack*)hs->Clone("hs2");

  c1->BuildLegend(0.7,0.75,0.88,0.88);
  TPad *p = new TPad("p1","p1",0.45,0.25,0.95,0.7);
  p->SetFillColorAlpha(kWhite, 0);
  p->Draw("same");
  p->cd();
  hs2->SetTitle("");
  hs2->Draw("p");
  hs2->GetHistogram()->GetYaxis()->SetLimits(0, 25000);
  hs2->GetHistogram()->GetYaxis()->SetRangeUser(0, 25000);
  hs2->GetHistogram()->SetMaximum(25000);
  hs2->GetYaxis()->SetLimits(0, 25000);
  hs2->GetYaxis()->SetRangeUser(0, 25000);
  hs2->Draw("pfc nostack");
  p->Update();
  
  
  return 0;
}
