// Copyright [2016] Misha Mikhasenko
// Discription

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "draw_distribution_from_tree.C"

TCanvas *draw_2pi_vs_3pi_for_mindf() {

  // extract binning
  const int NbinsX = 50;
  double binsX[NbinsX+1] = {0.5, 0.54, 0.58, 0.62, 0.66, 0.7, 0.74, 0.78, 0.82, 0.86, 0.9, 0.94, 0.98, 1.02, 1.06, 1.1, 1.14, 1.18, 1.22, 1.26, 1.3, 1.34, 1.38, 1.42, 1.46, 1.5, 1.54, 1.58, 1.62, 1.66, 1.7, 1.74, 1.78, 1.82, 1.86, 1.9, 1.94, 1.98, 2.02, 2.06, 2.1, 2.14, 2.18, 2.22, 2.26, 2.3, 2.34, 2.38, 2.42, 2.46, 2.5};
  const int NbinsY = 62;
  double binsY[NbinsY+1] = {0.27914, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.12, 1.16, 1.2, 1.24, 1.28, 1.32, 1.36, 1.4, 1.44, 1.48, 1.52, 1.56, 1.6, 1.64, 1.68, 1.72, 1.76, 1.8, 1.84, 1.88, 1.92, 1.96, 2, 2.04, 2.08, 2.12, 2.16, 2.2, 2.24, 2.28};
  TH2D *hh = new TH2D("hist3", "Intensity of 1^{++}0^{+} [#pi#pi]_{0^{++}}#pi P-wave;Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV);Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV)",
                      NbinsX, binsX, NbinsY, binsY);
  
  TCanvas *can = new TCanvas("c1", "can", 1000, 1000);

  const char *quantity = "sqrt(s1):sqrt(s)"; 
  const char *fproj_name = "~/cernbox/tmp/pwa_results/amplitudes/mfit_0.127471-0.144385.root"; //0.127471-0.144385
  // ~/Documents/phsp_with_COMPASS_basis/pwa_results.root
  const char *fin_templ = "/home/mikhasenko/Documents/phsp_with_COMPASS_basis/mc.e5/phsp.m%d.e5.root";
  // extract and draw
  for (uint bin = 1; bin <= 100; bin++) {  // 1--100
    draw_distribution_from_tree(TString::Format("%s>>+hist3", quantity),
                                fproj_name,
                                fin_templ,
                                bin,
                                {11, 16},
                                "false");
  }
  // devided by the bin width
  for (int i = 1; i <= NbinsX; i++)
    for (int j = 1; j <= NbinsY; j++)
      hh->SetBinContent(hh->GetBin(i,j),
                        hh->GetBinContent(hh->GetBin(i,j)) / hh->GetYaxis()->GetBinWidth(j)
                        );
  hh->SetStats(kFALSE);
  hh->Draw("colz");
  return can;
}
