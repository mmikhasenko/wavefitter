// Copyright [2016] Mikhail Mikhasenko

#include <list>
#include <string>
#include <vector>

#include "TMultiGraph.h"
#include "TGraphErrors.h"

#include "mstructures.hh"

using std::vector;

TGraphErrors *draw(const DP & data) {
  const uint nPoints = data.data.size();
  vector<double> x, y, dy;
  for (auto && dp : data.data) {
    x.push_back(dp.x);
    y.push_back(dp.y);
    dy.push_back(dp.dy);
  }
  TGraphErrors *lgr = new TGraphErrors(nPoints, x.data(), y.data(), 0, dy.data());
  lgr->SetTitle(data.title.c_str());
  return lgr;
}

TGraph *draw(std::function<double(double)> func, double lrange, double rrange, uint nPoints) {
  TGraph *lgr = new TGraph(nPoints);
  for (uint i = 0; i < nPoints; i++) {
    lgr->GetX()[i] = lrange + (rrange-lrange)/(nPoints-1)*i;
    lgr->GetY()[i] = func(lgr->GetX()[i]);
  }
  return lgr;
}

TGraph *style(TGraph *lgr, double color, double style) {
  lgr->SetLineColor(color);
  lgr->SetLineStyle(style);
  return lgr;
}

TMultiGraph *combine(std::list<TGraph*> grs) {
  TMultiGraph *lm = new TMultiGraph();
  for (auto && gr : grs) {
    lm->Add(gr, gr->GetDrawOption()); 
      std::cout << "----------" << gr->GetDrawOption() << "\n";
  }
  return lm;
}
TMultiGraph *combine(TGraph* g1) {
  return combine({g1});
}
TMultiGraph *combine(TGraph* g1, TGraph* g2) {
  return combine({g1, g2});
}
TMultiGraph *combine(TGraph* g1, TGraph* g2, TGraph* g3) {
  return combine({g1, g2, g3});
}
TMultiGraph *combine(TGraph* g1, TGraph* g2, TGraph* g3, TGraph* g4) {
  return combine({g1, g2, g3, g4});
}
