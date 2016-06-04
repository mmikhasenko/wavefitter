// Copyright [2016] Mikhail Mikhasenko

#include <mstructures.hh>

#include <functional>

#include <MmatrixK.h>
#include <MRelationHolder.h>

int main(int argc, char *argv[]) {
  MmatrixK k(2, 3, {1, 2, 1});
  const double pars[] = {1.1, 1.2, 1.4,  2, 5,  10, 11,  20, 22,  100,  7, 2,  9};
  k.setPars(pars);

  std::cout << k.CreatePool({"m1", "g1", "g2",
        "m2", "h1", "h2",
        "m3", "i1", "i2"
        }) << "\n";

  // k.PrintParameters();

  std::function<double(double)> f0 = [&](double s) -> double {return norm(k.getA(0, s));};

  MRelationHolder rels;
  DP points;
  // read points
  rels.AddRelation(points, f0);

  return 0;
}
