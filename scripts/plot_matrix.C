
#include "TMatrixD.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TROOT.h"

TCanvas *plot_matrix(uint iBin) {
  TMatrixD m(88, 88);
  for (int i = 0; i < 88; i++) {
    for (int j = 0; j < 88; j++) {
      TH1D *h = static_cast<TH1D*>(gDirectory->Get((i == j) ?
                                       TString::Format("h%d", i+1) :
                                       ((i < j) ?
                                        TString::Format("h1%03d%03d", i+1, j+1) :
                                        TString::Format("h1%03d%03d", j+1, i+1))));
      m(i, j) = h->GetBinContent(iBin);
    }
  }
  TCanvas *c1 = new TCanvas("c1", "title", 1000, 1000);
  m.Draw("colz");
  TH2D *h = static_cast<TH2D*>(gROOT->FindObject("TMatrixDBase"));
  if (h) {
    for (int i = 0; i < 88; i++) {
      TH1D *hdiag = static_cast<TH1D*>(gROOT->FindObject(TString::Format("h%d", i+1)));
      if (!hdiag) break;
      h->GetXaxis()->SetBinLabel(i+1, hdiag->GetTitle());
      h->GetYaxis()->SetBinLabel(i+1, hdiag->GetTitle());
    }
  }
  h->SetStats(kFALSE);
  h->GetXaxis()->SetLabelSize(0.02);
  h->GetYaxis()->SetLabelSize(0.02);
  h->GetZaxis()->SetRangeUser(0.01, 1.0);
  h->GetXaxis()->SetTickSize(0);
  h->GetYaxis()->SetTickSize(0);
  h->SetTitle("PW integral matrix");
  // h->Draw("colz");
  return c1;
}
