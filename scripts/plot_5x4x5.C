#include "TDirectory.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

uint plot_5x4x4() {
  TFile *_file0 = TFile::Open("pwa_hfit_0.100000-0.112853.root");
  TCanvas *c1 = new TCanvas("c1","title", 1500, 1200);
  c1->Divide(5,4);
  for (p = 0; p < 4; p++){
    for (uint i = 1; i <= 20; i++) {c1->cd(i); gDirectory->Get(TString::Format("h%d",20*p+i))->Draw();}
    c1->Print(TString::Format("/tmp/pwa%d.pdf",p));
  }
  return 1;
}

uint plot_together_5x4x4() {
  TFile *_file0 = TFile::Open("pwa_hfit_0.100000-0.112853.root");
  TFile *_file1 = TFile::Open("result.invert_matrix.root.thresholded.symm");
  TCanvas *c1 = new TCanvas("c1","title", 1500, 1200);
  c1->Divide(5,4);
  for (p = 0; p < 4; p++){
    for (uint i = 1; i <= 20; i++) {
      c1->cd(i);
      std::cout << "h" << 20*p+i << ", ";
      _file0->cd();
      TH1D *h0 = (TH1D*)gDirectory->Get(TString::Format("h%d",20*p+i));
      h0->Draw();
      _file1->cd();
      TH1D *h1 = (TH1D*)gDirectory->Get(TString::Format("h%d",20*p+i));
      for (uint bn = 1; bn <= h1->GetNbinsX(); bn++) {
        h1->SetBinContent(bn, h1->GetBinContent(bn)*h1->GetXaxis()->GetBinCenter(bn) );
      }
      h1->SetLineColor(kRed);
      h1->Scale(15.);
      h1->Draw("same");
    }
    const char *b = (p != 0 && p != 3) ? "pdf" : ((p == 0) ? "pdf(" : "pdf)");
    c1->Print(TString::Format("/tmp/plot_5x4x4.%s",b));
  }
  return 1;
}
