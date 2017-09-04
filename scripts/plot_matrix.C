
#include "TMatrixD.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TROOT.h"

TCanvas *plot_matrix(uint iBin, bool norm = false);
TCanvas *plot_matrix(uint iBin, bool norm) {
  TMatrixD m(88, 88);
  for (int i = 0; i < 88; i++) {
    for (int j = 0; j < 88; j++) {
      TH1D *h = static_cast<TH1D*>(gDirectory->Get((i == j) ?
                                       TString::Format("h%d", i+1) :
                                       ((i < j) ?
                                        TString::Format("h2%03d%03d", i+1, j+1) :
                                        TString::Format("h1%03d%03d", j+1, i+1))));
      m(i, j) = h->GetBinContent(iBin);
    }
  }
  // notmalized
  if (norm) {
    for (int i = 0; i < 88; i++)
      for (int j = 0; j < 88; j++)
        if (i != j) m(i, j) = m(i, j)/sqrt(m(i, i)*m(j, j));
    for (int j = 0; j < 88; j++) m(j, j) = 1;
  }
  // Draw matrix
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
  if (norm) h->SetTitle("PW integral matrix normalized by the diagonal");
  c1->SetBottomMargin(0.14);
  c1->SetLeftMargin(0.14);
  // c1->SetGrid();

  TText *re = new TText(10, 80, "Re");
  re->SetTextAlign(22);
  re->SetTextFont(43);
  re->SetTextSize(80);
  re->Draw();
  TText *im = new TText(75, 10, "Im");
  im->SetTextAlign(22);
  im->SetTextFont(43);
  im->SetTextSize(80);
  im->Draw();

  double M3pi = 0.51+0.02*(iBin-1);
  TLatex *sqrtS = new TLatex(80, 90, Form("M_{3#pi} = %1.2f GeV", M3pi));
  sqrtS->SetTextAlign(22);
  sqrtS->SetTextFont(43);
  sqrtS->SetTextSize(20);
  sqrtS->Draw();
  // h->Draw("colz");
  return c1;
}
