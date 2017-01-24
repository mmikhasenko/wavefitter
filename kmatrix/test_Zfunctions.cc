// Copyright [2016] Misha Mikhasenko

#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "TCanvas.h"

#include "M3bodyAngularBasis.h"
#include "mstructures.h"

int main() {
  TCanvas can("canva");
  combine(
          SET1(
               draw([](double theta)->double{return real(Math::ZJMLS(3, 0, 2, 0, theta, 0., 1., 0.));},
                    0, M_PI, 100), SetLineColor(kBlack)),
          SET1(
               draw([](double theta)->double{return imag(Math::ZJMLS(3, 0, 2, 0, theta, 0., 1., 0.));},
                    0, M_PI, 100), SetLineColor(kRed)) )->Draw("al");
  can.Print("/tmp/test_Zfunctions.pdf");
  return 0;
}

