// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MSTRUCTURES_H_
#define SRC_MSTRUCTURES_H_

#include <list>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

typedef struct {
  double x;
  double y;
  double dy;
} data_point;

typedef struct {
  std::list<data_point> data;
  std::string name;
  std::string title;
  double lrange;
  double rrange;
} data_points;

typedef data_points DP;

typedef struct {
  const DP & data;
  std::function<double(double)> func;
} relation;

// relation make_relation(const DP & data, double (*func)(double s, int flag));

TGraphErrors *draw(const DP & data);

TGraph *draw(std::function<double(double)> funct,
             double lrange, double rrange, uint nPoints = 100);

TGraph *draw(const std::vector<double> & xv, const std::vector<double> & yv);

TGraph *style(TGraph *lgr, double color, double style = 1);

#define SET1(f, a) ([&]()->TGraph* {TGraph *g = f; g->a; return g;})()
#define SET2(f, a, b) ([&]()->TGraph* {TGraph *g = f; g->a; g->b; return g;})()
#define SET3(f, a, b, c) ([&]()->TGraph* {TGraph *g = f; g->a; g->b; g->c; return g;})()
#define SET4(f, a, b, c, d) ([&]()->TGraph* {TGraph *g = f; g->a; g->b; g->c; g->d; return g;})()
#define SET5(f, a, b, c, d, e) ([&]()->TGraph* {TGraph *g = f; g->a; g->b; g->c; g->d; g->e; return g;})()

TMultiGraph *combine(std::list<TGraph*> grs);
TMultiGraph *combine(TGraph* g1, TGraph* g2);
TMultiGraph *combine(TGraph* g1, TGraph* g2, TGraph* g3);
TMultiGraph *combine(TGraph* g1, TGraph* g2, TGraph* g3, TGraph* g4);
TMultiGraph *combine(TGraph* g1, TGraph* g2, TGraph* g3, TGraph* g4, TGraph* g5);
TMultiGraph *combine(TGraph* g1, TGraph* g2, TGraph* g3, TGraph* g4, TGraph* g5, TGraph* g6);

template<class Array>
TGraph *draw(const Array &a) {
  const uint nPoints = a.size();
  TGraph *lgr = new TGraph(nPoints);
  for (int i = 0; i < nPoints; i++) {
    lgr->GetX()[i] = a[i].first;
    lgr->GetY()[i] = a[i].second;
  }
  lgr->SetTitle("");
  return lgr;
}

#endif  // SRC_MSTRUCTURES_H_

