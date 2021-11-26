// Copyright [2016] Mikhail Mikhasenko

#ifndef SRC_DEFLIB_H_
#define SRC_DEFLIB_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <vector>
#include <utility>
#include <iostream>

#define LAMBDA(a, b, c) ((a)*(a)+(b)*(b)+(c)*(c)-2.*(a)*(b)-2.*(b)*(c)-2.*(c)*(a))
#define sqrtPi(x) sqrt(x*exp(-cd(0, 1.)*M_PI))*exp(cd(0, 1.)*M_PI/2.)
#define sqrtPhi(x, phi) sqrt(x*exp(-cd(0, 1.)*phi))*exp(cd(0, 1.)*phi)
#define POW2(a) ((a)*(a))
#define RHO(a, b, c) (1./(8*M_PI)*sqrt(LAMBDA(a, b, c))/(a))
#define RHO_PI(a, b, c) (1./(8*M_PI)*sqrtPi(LAMBDA(a, b, c))/(a))

typedef std::complex<double> cd;

cd ChewMandelstam(cd s, double m1sq, double m2sq);

// technical tools

/*************************************************/
/* The method was a bit slower then the next one */
/*************************************************/

// template <typename Type>
// Type getvalue(double M, const std::vector<std::pair<double, Type> > &table) {
//   // some super cases
//   if (M == table[0].first) return table[0].second;
//   if (M == table[table.size()-1].first) return table[table.size()-1].second;
//   // standard search
//   const uint N = table.size();
//   const double lft = table[0].first;
//   const double Mstep = table[1].first - lft;
//   int Nsteps = (M - lft)/Mstep;
//   if (Nsteps == static_cast<int>(N-1) && M < table[N-1].first) Nsteps--; /* presision fix */
//   if (Nsteps < 0 || Nsteps >= static_cast<int>(-1+N)) {
//     std::cerr << "Error!! in getvalue! M = " << M << " not in ["
//               << table[0].first << ", "
//               << table[N-1].first << "], diff to max = " << M - table[N-1].first
//               << std::endl;
//     return 0;
//   }
//   const Type value = table[Nsteps].second +
//     (table[Nsteps+1].second - table[Nsteps].second) /
//     (table[Nsteps+1].first  - table[Nsteps].first) * (M - table[Nsteps].first);
//   return value;
// }
// 
template <typename Type>
Type getvalue(double M, const std::pair<double, Type> *table, uint N) {
  // some super cases
  if (M == table[0].first) return table[0].second;
  if (M == table[N-1].first) return table[N-1].second;
  // standard search
  const double lft = table[0].first;
  const double Mstep = table[1].first - lft;
  int Nsteps = (M - lft)/Mstep;
  if (Nsteps == static_cast<int>(N-1) && M < table[N-1].first) Nsteps--; /* presision fix */
  if (Nsteps < 0 || Nsteps >= static_cast<int>(-1+N)) {
    std::cerr << "Error!! in getvalue! M = " << M << " not in ["
              << table[0].first << ", "
              << table[N-1].first << "], diff to max = " << M - table[N-1].first
              << std::endl;
    return 0;
  }
  const Type value = table[Nsteps].second +
    (table[Nsteps+1].second - table[Nsteps].second) /
    (table[Nsteps+1].first  - table[Nsteps].first) * (M - table[Nsteps].first);
  return value;
}

#endif  // SRC_DEFLIB_H_
