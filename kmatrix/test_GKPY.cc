// Copyright [2016] Misha Mikhasenko

#include <iostream>
#include <vector>
#include <fstream>

#include "TGraph.h"
#include "TCanvas.h"
#include "waves.h"
#include "mstructures.h"
#include "constants.h"

#include "MIsobar.h"

int main(int ac, char** av) {

  /* generate range */
  const uint Npoints = 100;
  std::vector<double> sVls(Npoints);
  std::pair<double, double> ranges(2*PI_MASS+1e-6, 2.2);
  for (uint i=0; i < sVls.size(); i++) sVls[i] = ranges.first + (ranges.second-ranges.first)/(Npoints-1)*i;

  /* data arrays */
  std::vector<cd> amp(Npoints);
  std::vector<double> reA(Npoints);
  std::vector<double> imA(Npoints);
  std::vector<double> deltaA(Npoints);

  /* fill rho amplitude */
  /* rho */
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd { return waves::GKPY::T(e*e); });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });
  std::transform(sVls.begin(), sVls.end(), deltaA.begin(), [&](double e)->double { return waves::GKPY::phi(e*e); });

  std::ofstream myfile;
  myfile.open("/tmp/GKPY.misha.txt");
  for (uint i=0; i < sVls.size(); i++) myfile << sVls[i] << " " << real(amp[i]) << " "  << imag(amp[i]) << "\n";
  myfile.close();

  TCanvas c1("isobar_shape", "isobar_shape", 1000, 500);
  c1.Divide(2, 1);
  c1.cd(1);
  combine(
          SET3(draw(sVls, reA),
               SetTitle("Real part"),
               SetFillStyle(0),
               SetLineColor(kBlue) ),
          SET3(draw(sVls, imA),
               SetTitle("Imag part"),
               SetFillStyle(0),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.cd(1)->BuildLegend();
  c1.cd(2);
  SET3(draw(sVls, deltaA),
       SetTitle("Phase"),
       SetFillStyle(0),
       SetLineColor(kRed) )->Draw("al");
            
  c1.SaveAs("/tmp/test_GKPY.png");

  return 0.;
}
