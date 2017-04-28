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


uint plot_whatever(uint nPads = 4,
                   const char *file_name = "/tmp/result.invert_matrix.root",
                   uint MAX = 88,
                   const char *tmpl = "h%d");
uint plot_whatever(uint nPads, const char *file_name, uint MAX, const char *tmpl) {
  // TFile *fin1 = new TFile("/mnt/data/compass/2008/integrals_stefan/integrals.root");
  // if (!fin1) { std::cerr << "Error with file1!\n"; return 1; }
  TFile *fin2 = new TFile(file_name);
  if (!fin2) { std::cerr << "Error with file2!\n"; return 1; }

  TCanvas *c1 = new TCanvas("c1");
  c1->DivideSquare(nPads);
  for (uint p = 0; p < nPads; p++) {
    const uint index = 1+gRandom->Integer(MAX);
    // fin1->cd(); TGraph *gr; gDirectory->GetObject(TString::Format("g%d", index), gr);
    // if (!gr) { std::cerr << "Error with graph!\n"; return 1; }
    fin2->cd(); TH1D *hh = static_cast<TH1D*>(gDirectory->Get(TString::Format(tmpl, index)));
    if (!hh) { std::cerr << "Error with hist h" << index << "!\n"; return 1; }
    hh->SetStats(kFALSE);

    // double scale1 = gr->GetY()[179];
    // double scale2 = hh->GetBinContent(90);
    // for (uint i = 0; i < gr->GetN(); i++) gr->GetY()[i] *= scale2/scale1;

    c1->cd(p+1);
    hh->Draw();
    // gr->Draw("same");
  }
  return 0;
}
