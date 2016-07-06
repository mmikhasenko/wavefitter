// Copyright [2016] Mikhail Mikhasenko

#include <functional>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TGraph.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH1D.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "mstructures.hh"

TH1D *draw(const std::vector<double> &v) {
  double max = *std::max_element(v.begin(), v.end());  // , [](double v1, double v2)->bool{return (v1>v2);}
  double min = *std::min_element(v.begin(), v.end());
  TH1D *h = new TH1D("h1","title", 100, min, max);
  std::for_each(v.begin(), v.end(), [&](double x)->void{h->Fill(x);});
  return h;
}

int main(int argc, char *argv[]) {

  const uint Np = 4;
  double x[Np] = {-1, 0, 1, 2};
  double y[Np] = {1., 2.5, 0.3, 2.5};
  double dy[Np] = {0.1, 3., 0.1, 0.1};

  TGraphErrors ge(Np, x, y, 0, dy);
  
  TCanvas c1("c1");
  TF1 *f1 = new TF1("f1_with",    "[0]+[1]*x+[2]*x*x", -1.1, 2.1);
  TF1 *f2 = new TF1("f2_without", "[0]+[1]*x+[2]*x*x", -1.1, 2.1);

  ge.Fit(f2,"WN");
  ge.Fit(f1,"EN");
  SET2(&ge, SetMarkerStyle(8), SetTitle("Test of a fit with/without errors"))->Draw("ap");
  f1->SetLineColor(kRed); f1->Draw("same");
  f2->SetLineColor(kRed); f2->SetLineStyle(2); f2->Draw("same");
  c1.SaveAs("/tmp/a.test_Errors.pdf");

  const uint Natt = 10000;
  std::vector<double> pw[3];
  std::vector<double> pwo[3];
  for (uint e=0; e < Natt; e++) {
    // redistribute data
    for (uint i=0; i < Np; i++) ge.GetY()[i] = gRandom->Gaus(y[i], dy[i]);
    // fit again
    ge.Fit(f2, "WNQ");
    ge.Fit(f1, "ENQ");
    // push
    for (uint i = 0; i < 3; i++)  pw[i].push_back(f1->GetParameter(i));
    for (uint i = 0; i < 3; i++) pwo[i].push_back(f2->GetParameter(i));
  }
  // plot
  c1.Clear(); c1.Divide(2, 2);
  for (uint i=0; i < 3; i++) {
    c1.cd(i+1);
    TH1D *hw  = draw(pwo[i]); hw ->SetFillColor(kGray);   hw ->SetTitle(TString::Format("par[%d] distribution", i)); hw ->Draw();
    TH1D *hwo = draw( pw[i]); hwo->SetFillColor(kOrange); hwo->Draw("same");
  }
  c1.SaveAs("/tmp/b.test_Errors.pdf");
  
  return 0;
}
