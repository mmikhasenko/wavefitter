// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_MSTRUCTURES_H_
#define SRC_MSTRUCTURES_H_

#include <list>
#include <string>
#include <vector>
#include <functional>

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

TGraph *draw(double (*func)(double s),
             double lrange, double rrange, uint nPoints);

TGraph *style(TGraph *lgr, double color, double style = 1);

TMultiGraph *combine(std::list<TGraph*> grs);
TMultiGraph *combine(TGraph* g1, TGraph* g2);
TMultiGraph *combine(TGraph* g1, TGraph* g2, TGraph* g3);
TMultiGraph *combine(TGraph* g1, TGraph* g2, TGraph* g3, TGraph* g4);

#endif  // SRC_MSTRUCTURES_H_

