// Copyright [2016] Misha Mikhasenko
// Discription:
//  It draws Deck expansion

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "v5/TFormula.h"

void draw_distribution_from_tree(const char *request,
                                 const char *fproj_name,
                                 const char *fin_templ,
                                 uint bin = 1,
                                 const std::vector<uint> _v = {},
                                 bool SYMM = false,
                                 const char *cut = "s>0");

void draw_distribution_from_tree(const char *request,
                                 const char *fproj_name,
                                 const char *fin_templ,
                                 uint bin,
                                 const std::vector<uint> _v,
                                 bool SYMM,
                                 const char *cut) {
  std::vector<uint> v(_v);
  if (v.size() == 0) {
    v.resize(88); for (uint i = 0; i < 88; i++) v[i] = i+1;
  }
  const uint Nwaves = v.size();

  TFile *fproj = new TFile(fproj_name);
  double bamp[Nwaves][2];
  for (uint w = 0; w < v.size(); w++) {
    TH1D *hr = 0; gDirectory->GetObject(TString::Format("hr%d", v[w]), hr); if (!hr) {std::cerr << "No hr" << v[w] << "\n!"; return;}
    TH1D *hi = 0; gDirectory->GetObject(TString::Format("hi%d", v[w]), hi); if (!hi) {std::cerr << "No hi" << v[w] << "\n!"; return;}
    bamp[w][0] = hr->GetBinContent(bin);
    bamp[w][1] = hi->GetBinContent(bin);
  }
  fproj->Close();

  // work on tree
  TFile *fin = new TFile(TString::Format(fin_templ, bin-1));
  TTree *tin = 0; gDirectory->GetObject("angles", tin); if (!tin) {std::cerr << "No tree!\n"; return;}
  TString weight(TString::Format("(%s)*(", cut));
  ROOT::v5::TFormula::SetMaxima(100000, 100000, 100000);
  //*** real part
  if (SYMM) {
    for (uint w = 0; w < v.size(); w++) weight += TString::Format("(amp%d_frame1_real+amp%d_frame3_real)*(%f)+", v[w], v[w], bamp[w][0]);
    for (uint w = 0; w < v.size(); w++) weight += TString::Format("(-1.)*(amp%d_frame1_imag+amp%d_frame3_imag)*(%f)+", v[w], v[w], bamp[w][1]);
    weight += "0)**2 / 2";
   } else {
    for (uint w = 0; w < v.size(); w++) weight += TString::Format("amp%d_frame1_real*(%f)+", v[w], bamp[w][0]);
    for (uint w = 0; w < v.size(); w++) weight += TString::Format("(-1.)*amp%d_frame1_imag*(%f)+", v[w], bamp[w][1]);
    weight += "0)**2";
  }
  gROOT->cd();
  tin->Draw(request, weight);
  // std::cout << "first request: "<< request << "\n";
  //**** imag part
  weight = TString::Format("(%s)*(", cut);
  if (SYMM) {
    for (uint w = 0; w < v.size(); w++) weight += TString::Format("(amp%d_frame1_real+amp%d_frame3_real)*(%f)+", v[w], v[w], bamp[w][1]);
    for (uint w = 0; w < v.size(); w++) weight += TString::Format("(amp%d_frame1_imag+amp%d_frame3_imag)*(%f)+", v[w], v[w], bamp[w][0]);
    weight += "0)**2 / 2.";
  } else {
    for (uint w = 0; w < v.size(); w++) weight += TString::Format("amp%d_frame1_real*(%f)+", v[w], bamp[w][1]);
    for (uint w = 0; w < v.size(); w++) weight += TString::Format("amp%d_frame1_imag*(%f)+", v[w], bamp[w][0]);
    weight += "0)**2";
  }
  // change the request "s1>>h3(2000)" to "s1>>+h3"
  TString srequest(request);
  if (srequest.Last('+') != srequest.First(">>")+2)
    srequest.Insert(srequest.First(">>")+2, "+");
  if (srequest.First(">>") < srequest.Last('('))
    srequest = srequest(0, srequest.Last('('));
  // std::cout << "second request: "<< srequest << "\n";
  tin->Draw(srequest, weight);
  fin->Close();
}
