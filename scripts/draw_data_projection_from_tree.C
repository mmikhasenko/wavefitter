// Copyright [2016] Misha Mikhasenko
// Discription:
//  It draws Deck expansion

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TCanvas.h"

TH1D *draw_deck_projection_from_tree(const char *request,
                                     const char *fproj_name,
                                     const char *fin_templ,
                                     uint bin) {
  const uint Nwaves = 88;
  
  TFile *fproj = new TFile(fproj_name);
  TH1D *hproj[Nwaves][2];
  for (uint w = 1; w < Nwaves; w++) {
    hproj[w][0] = 0; hproj[w][1] = 0;
    gDirectory->GetObject(TString::Format("hr%d", w+1), hproj[w][0]); if (!hproj[w][0]) return 0;
    gDirectory->GetObject(TString::Format("hi%d", w+1), hproj[w][1]); if (!hproj[w][1]) return 0;
  }
  double bamp[Nwaves][2];
  for (uint w = 1; w < Nwaves; w++) {
    bamp[w][0] = hproj[w][0]->GetBinContent(bin);
    bamp[w][1] = hproj[w][1]->GetBinContent(bin);
  }
  fproj->Close();

  // work on tree
  TFile *fin = new TFile(TString::Format(fin_templ, bin-1));
  TTree *tin = 0; gDirectory->GetObject("angles", tin); if (!tin) return 0;
  TString weight("(");
  // real part
  for (uint w = 1; w < Nwaves; w++) weight += TString::Format("amp%d_frame1_real*(%f)+", w+1, bamp[w][0]);
  weight += "0)**2+(";
  for (uint w = 1; w < Nwaves; w++) weight += TString::Format("amp%d_frame1_imag*(%f)+", w+1, bamp[w][1]);
  weight += "0)**2";
  tin->Draw(TString::Format("%s>>hgen", request), weight);
  TH1D *hgen = static_cast<TH1D*>(gROOT->FindObject("hgen"));
  return hgen;
}
