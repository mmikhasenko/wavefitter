// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace b = boost::numeric::ublas;

int main(int argc, char *argv[]) {
  b::matrix<double> b1(5, 5);
  b::matrix<double> b2(5, 5);
  b::matrix<double> b3 = prod(b1, b2);
  std::cout << b3 << "\n";
  return 0;
}
