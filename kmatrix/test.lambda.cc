// Copyright [2016] Mikhail Mikhasenko

int f(double *x) {
  return x[0]+x[1]+x[2];
}

#include <vector>
#include <iostream>
#include <functional>
#include <algorithm>

#include "mintegrate.hh"

int ff(const std::function<double(double*)> &amp, double *x) {
  return (amp)(x);
}


int main(int argc, char **argv) {
  double factor = 5;
  auto func0 = [](double *x) -> int { return f(x); };
  auto func1 = [&factor](double *x) -> int { return factor*f(x); };

  auto fl2 = func1;
  int (*lf)(double *x) = func0;
  {
    lf = [](double *x) -> int { return x[1]+x[0]-x[2]; };
  }
  double arr[3] = {1, 2, 3};
  std::cout << lf(arr) << "\n";
  std::cout << ff([=](double *p)->double {return f(p);}, arr) << "\n";
  // test 2
  std::vector<double> aff = {3, 2, 4, 767, 4, 2, 5, 2, 7, 34, 2};

  struct ex { static bool f(int a, int b) { return a < b; } };  // create function to sort
  std::sort(aff.begin(), aff.end(), ex::f);

  std::for_each(aff.begin(), aff.end(), [&](int i) {std::cout << i << " ";});
  std::cout << "\n";

  // --------------------------------------------------------------------------------
  // Now we check what is the best way to speed up integration

  // test1: struct
  if (argc == 1) {
    for (int i=0; i < 1e7; i++) {
      struct com { double c; double operator() (double x) const { return c*x*x + 1.; } } ex1;  // create function to sort
      ex1.c = i;
      double val = tintegrate(ex1, 0., 1.);
      if (i == 0 || i == 10) std::cout << val << "\n";
    }
  } else if (argc == 2) {
    // test2: lambda
    for (int i=0; i < 1e7; i++) {
      std::function<double(double)> ff = [&](double x)->double{return i*x*x+1;};
      double val = integrate(ff, 0., 1.);
      if (i == 0 || i == 10) std::cout << val << "\n";
    }
  } else if (argc == 3) {
    // test3: direct lambda
    for (int i=0; i < 1e7; i++) {
      double val = integrate([&](double x)->double{return i*x*x+1;}, 0., 1.);
      if (i == 0 || i == 10) std::cout << val << "\n";
    }
  } else if (argc == 4) {
    // test4: lambda + template
    for (int i=0; i < 1e7; i++) {
      double val = tintegrate([&](double x)->double{return i*x*x+1;}, 0., 1.);
      if (i == 0 || i == 10) std::cout << val << "\n";
    }
  }
  return 0;
}
