// Copyright [2016] Mikhail Mikhasenko

#include <list>
#include <string>
#include <vector>

#include "TMultiGraph.h"
#include "TGraphErrors.h"

#include "mstructures.hh"

// relation make_relation(const DP & data, double (*func)(double s, int flag)) {
//   relation lrel;
//   lrel.data = data;
//   lrel.func = func;
//   return lrel;
// }
using std::vector;


TGraphErrors *draw(const DP & data) {
  const uint nPoints = data.data.size();
  vector<double> x, y, dy;
  for (auto && dp : data.data)
    if (dp.x > data.lrange && dp.x < data.rrange) {
      x.push_back(dp.x);
      y.push_back(dp.y);
      dy.push_back(dp.dy);
    }
  TGraphErrors *lgr = new TGraphErrors(nPoints, x.data(), y.data(), 0, dy.data());
  lgr->SetTitle(data.title.c_str());
  return lgr;
}

TGraph *draw(double (*func)(double s), double lrange, double rrange, uint nPoints) {
  double x[nPoints], y[nPoints], dy[nPoints];
  for (int i = 0; i < nPoints; i++) {
    x[i] = (rrange-lrange)/(nPoints-1)*i;
    y[i] = func(x[i]);
  }
  TGraph *lgr = new TGraph(nPoints, x, y);
}

TGraph *style(TGraph *lgr, double color, double style) {
  lgr->SetLineColor(color);
  lgr->SetLineStyle(style);
  return lgr;
}

TMultiGraph *combine(std::list<TGraph*> grs) {
  TMultiGraph *lm = new TMultiGraph();
  for (auto && gr : grs) lm->Add(gr);
  return lm;
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

