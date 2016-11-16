// Copyright [2016] Misha Mikhasenko

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "TWigner.h"
#include "dFunction.hpp"

#include "deflib.h"

#include "TRandom.h"

void test1();
void test2();

int main() {
  // Wigner d-functions
  std::cout << "\n-----> Test 1: Wigner d-functions: \n\n";
  test1();
  // Wigner D-functions
  std::cout << "\n-----> Test 2: Wigner D-functions: \n\n";
  test2();

  return 0.;
}


void test1() {
  const int Nent = 1000000;
  for (int e=0; e < Nent; e++) {
    int j = gRandom->Integer(10);
    int m1 = gRandom->Integer(2*j)-j;
    int m2 = gRandom->Integer(2*j)-j;
    double theta = M_PI*gRandom->Rndm();
    double rpwa = rpwa::dFunction<double>(2*j, 2*m1, 2*m2, theta);
    double wigd = WignerD(j, m1, m2, theta);
    if (fabs(rpwa - wigd) > 1e-10 || (e % (Nent/10)) == 0) {
      if ((e % (Nent/10)) == 0) std::cout << "Good! ";
      else std::cout << "Error<> : ";
      std::cout << "{" << j << m1 << m2 << "}, " << theta
                << " : " << rpwa << " " << wigd << "\n";
    }
  }
}

void test2() {
  const int Nent = 1000000;
  for (int e=0; e < Nent; e++) {
    int j = gRandom->Integer(10);
    int m1 = gRandom->Integer(2*j)-j;
    int m2 = gRandom->Integer(2*j)-j;
    double theta = M_PI*gRandom->Rndm();
    double phi   = 2*M_PI*gRandom->Rndm();
    double psi   = 2*M_PI*gRandom->Rndm();
    cd rpwa = rpwa::DFunction<cd>(2*j, 2*m1, 2*m2, phi, theta, psi);
    cd wigd = WignerD(j, m1, m2, phi, theta, psi);
    if (fabs(rpwa - wigd) > 1e-10 || (e % (Nent/10)) == 0) {
      if ((e % (Nent/10)) == 0) std::cout << "Good! ";
      else std::cout << "Error<> : ";
      std::cout << "{" << j << m1 << m2 << "}, (" << phi << "," << theta << "," << psi << ")"
                << " : " << rpwa << " " << wigd << "\n";
    }
  }
}
