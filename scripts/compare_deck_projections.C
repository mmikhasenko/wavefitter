#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TGraph.h"
#include <iostream>
#include <vector>

void compare_deck_projections() {

  TFile *_file0 = TFile::Open("~/cernbox/tmp/linal_hfit_0.100000-0.112853.root");
  if (!_file0) { std::cerr << "no file0\n"; return; }
  TFile *_file1 = TFile::Open("~/cernbox/tmp/pwa_hfit_0.100000-0.112853.root");
  if (!_file1) { std::cerr << "no file1\n"; return; }
  TFile *_file2 = TFile::Open("~/cernbox/tmp/anal.proj_0.100000-0.112853.root");
  if (!_file2) { std::cerr << "no file2\n"; return; }
  std::vector<uint> num = {2,33,26,45,49,56};
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->DivideSquare(num.size());
  std::cout << "cc " << num.size() << ", " << num[2] << "\n";
  for(uint i = 0; i < num.size(); i++) {
    THStack *hs = new THStack(TString::Format("hs%d",num[i]),"");
    _file0->cd();
    c1->cd(i+1);
    TH1F *h0 = (TH1F*)gDirectory->Get(TString::Format("h%d",num[i]));
    if (!h0) { std::cerr << "no h0\n"; return; }
    h0->Scale(1e6);
    hs->Add(h0);

    _file1->cd();
    c1->cd(i+1);
    TH1F *h1 = (TH1F*)gDirectory->Get(TString::Format("h%d",num[i]));
    if (!h1) { std::cerr << "no h1\n"; return; }
    hs->Add(h1);
    h1->Scale(1);
    hs->SetTitle(h0->GetTitle());

    _file2->cd();
    c1->cd(i+1);
    TGraph *h2 = (TGraph*)gDirectory->Get(TString::Format("h%d",num[i]));
    if (!h2) { std::cerr << "no h2\n"; return; }
    h2->SetFillColor(kRed-1);
    //h2->SetFillColorAlpha(kRed-1, 0.7);
    // hs->Add(h2);

    hs->Draw("nostack");
    h2->Draw("same");
  }
}
