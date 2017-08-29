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


TCanvas *compare_integrals(uint nPads,
                           const char *hist_file,
                           const char *gr_file = "/mnt/data/compass/2008/integrals_stefan/integrals.root");

TCanvas *compare_integrals(uint nPads,
                           const char *hist_file,
                           const char *gr_file) {
  TFile *fin1 = new TFile(gr_file);
  TFile *fin2 = new TFile(hist_file);
  if (!fin1) { std::cerr << "Error with file1!\n"; return 0; }
  if (!fin2) { std::cerr << "Error with file2!\n"; return 0; }

  TCanvas *c1 = new TCanvas("c1");
  c1->DivideSquare(nPads);
  for (uint p = 0; p < nPads; p++) {
    const uint index = 1+gRandom->Integer(87);
    fin1->cd(); TGraph *gr; gDirectory->GetObject(TString::Format("g%d", index), gr);
    if (!gr) { std::cerr << "Error with graph!\n"; return 0; }
    fin2->cd(); TH1D *hh; gDirectory->GetObject(TString::Format("h%d", index), hh);
    if (!hh) { std::cerr << "Error with hist!\n"; return 0; }
    hh->SetStats(kFALSE);

    
    double scale1 = 0; for (int i=0; i < 200; i++) scale1 += gr->GetY()[i]/200;
    double scale2 = 0; for (int i=1; i <= 100; i++) scale2 += hh->GetBinContent(i)/100.;
    for (uint i = 0; i < gr->GetN(); i++) gr->GetY()[i] *= scale2/scale1;

    c1->cd(p+1);
    hh->Draw();
    gr->Draw("same");
  }
  return c1;
}
