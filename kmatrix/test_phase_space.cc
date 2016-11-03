// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <vector>
#include <fstream>
// #include <>

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "MIsobarChannel.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {
  // stable
  MIsobar rho_iso(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
  MIsobar  f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.);
  MIsobarPiPiS  pipiS_iso;

  /* generate range */
  const uint Npoints = 100;
  std::vector<double> sVls(Npoints);
  std::pair<double, double> ranges(2*PI_MASS, 3.0);
  for (uint i=0; i < sVls.size(); i++) sVls[i] = ranges.first + (ranges.second-ranges.first)/(Npoints-1)*i;

  /* data arrays */
  std::vector<cd> amp(Npoints);
  std::vector<double> reA(Npoints);
  std::vector<double> imA(Npoints);

  /* fill rho amplitude */
  /* rho */
  MIsobarChannel RhoPiS(rho_iso, PI_MASS, 0);
  RhoPiS.makeLookupTable(RhoPiS.sth(), 10., 100);
  RhoPiS.makeDisperseLookupTable(0.01, 10., 100);
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd { return RhoPiS.rholtilde(e*e); });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  std::ofstream myfile;
  myfile.open("/tmp/RhoPiS.misha.txt");
  for (uint i=0; i < sVls.size(); i++) myfile << sVls[i] << " " << real(amp[i]) << " "  << imag(amp[i]) << "\n";
  myfile.close();

  TCanvas c1("isobar_shape");
  combine(
          SET3(draw(sVls, reA),
               SetTitle("Real part"),
               SetFillStyle(0),
               SetLineColor(kBlue) ),
          SET3(draw(sVls, imA),
               SetTitle("Imag part"),
               SetFillStyle(0),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_phase_space.pdf(");

  /* f2 */
  MIsobarChannel F2PiS(f2_iso, PI_MASS, 0);
  F2PiS.makeLookupTable(F2PiS.sth(), 10., 100);
  F2PiS.makeDisperseLookupTable(0.01, 10., 100);
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd { return F2PiS.rholtilde(e*e); });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  myfile.open("/tmp/F2PiS.misha.txt");
  for (uint i=0; i < sVls.size(); i++) myfile << sVls[i] << " " << real(amp[i]) << " "  << imag(amp[i]) << "\n";
  myfile.close();

  combine(
          SET3(draw(sVls, reA),
               SetTitle("Real part"),
               SetFillStyle(0),
               SetLineColor(kBlue) ),
          SET3(draw(sVls, imA),
               SetTitle("Imag part"),
               SetFillStyle(0),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_phase_space.pdf");

  /* pipiS */
  MIsobarChannel PiPiSPiS(pipiS_iso, PI_MASS, 0);
  PiPiSPiS.makeLookupTable(PiPiSPiS.sth(), 10., 100);
  PiPiSPiS.makeDisperseLookupTable(0.01, 10., 100);
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd { return PiPiSPiS.rholtilde(e*e); });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  myfile.open("/tmp/PiPiSPiS.misha.txt");
  for (uint i=0; i < sVls.size(); i++) myfile << sVls[i] << " " << real(amp[i]) << " "  << imag(amp[i]) << "\n";
  myfile.close();

  combine(
          SET3(draw(sVls, reA),
               SetTitle("Real part"),
               SetFillStyle(0),
               SetLineColor(kBlue) ),
          SET3(draw(sVls, imA),
               SetTitle("Imag part"),
               SetFillStyle(0),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.BuildLegend();
  c1.SaveAs("/tmp/test_phase_space.pdf)");

  return 0;
}
