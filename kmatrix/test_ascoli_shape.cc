// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <vector>
#include <fstream>
// #include <>

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "MAscoli.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {

  // Deck parameters
  uint J = 2;
  int M = 0;
  double mAsq = POW2(PI_MASS);
  double mBsq = POW2(PROT_MASS);
  double mDsq = POW2(PROT_MASS);
  double mtRsq = POW2(PI_MASS);
  double stot = 2*190*PROT_MASS;
  double t = -0.01;
  double R = 5.;

  /* generate range */
  const uint Npoints = 100;
  std::vector<double> sVls(Npoints);
  std::pair<double, double> ranges;

  /* data arrays */
  std::vector<cd> amp(Npoints);
  std::vector<double> reA(Npoints);
  std::vector<double> imA(Npoints);

  /************************************************************************/
  /* rho pi P - wave */
  double iso_mass = RHO_MASS;
  uint S = 1;
  uint L = 1;
  // fill X array
  ranges.first  = iso_mass + PI_MASS;
  ranges.second = 2.5;
  for (uint i=0; i < sVls.size(); i++) sVls[i] = ranges.first + (ranges.second-ranges.first)/(Npoints-1)*i;
  // calculate value of the function
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd {
      return MAscoli::getProjectedReducedDeck(J, M, L,
                                              POW2(iso_mass), S, R,
                                              e*e, t,
                                              mtRsq,
                                              stot,
                                              mAsq, mBsq, mDsq,
                                              POW2(PI_MASS));
    });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  std::ofstream myfile;
  myfile.open("/tmp/ascoli.RhoPiP.misha.txt");
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
  c1.SaveAs("/tmp/test_ascoli_shape.pdf(");

  /************************************************************************/
  /* rho pi F - wave */
  iso_mass = RHO_MASS;
  S = 1;
  L = 3;
  // fill X array
  ranges.first  = iso_mass + PI_MASS;
  ranges.second = 2.5;
  for (uint i=0; i < sVls.size(); i++) sVls[i] = ranges.first + (ranges.second-ranges.first)/(Npoints-1)*i;
  // calculate value of the function
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd {
      return MAscoli::getProjectedReducedDeck(J, M, L,
                                              POW2(iso_mass), S, R,
                                              e*e, t,
                                              mtRsq,
                                              stot,
                                              mAsq, mBsq, mDsq,
                                              POW2(PI_MASS));
    });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  myfile.open("/tmp/ascoli.RhoPiF.misha.txt");
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
  c1.SaveAs("/tmp/test_ascoli_shape.pdf");

  /************************************************************************/
  /* f2 pi S - wave */
  iso_mass = F2_MASS;
  S = 2;
  L = 0;
  // fill X array
  ranges.first  = iso_mass + PI_MASS;
  ranges.second = 2.5;
  for (uint i=0; i < sVls.size(); i++) sVls[i] = ranges.first + (ranges.second-ranges.first)/(Npoints-1)*i;
  // calculate value of the function
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd {
      return MAscoli::getProjectedReducedDeck(J, M, L,
                                              POW2(iso_mass), S, R,
                                              e*e, t,
                                              mtRsq,
                                              stot,
                                              mAsq, mBsq, mDsq,
                                              POW2(PI_MASS));
    });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  myfile.open("/tmp/ascoli.F2PiS.misha.txt");
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
  c1.SaveAs("/tmp/test_ascoli_shape.pdf");

  /************************************************************************/
  /* f2 pi D - wave */
  iso_mass = F2_MASS;
  S = 2;
  L = 2;
  // fill X array
  ranges.first  = iso_mass + PI_MASS;
  ranges.second = 2.5;
  for (uint i=0; i < sVls.size(); i++) sVls[i] = ranges.first + (ranges.second-ranges.first)/(Npoints-1)*i;
  // calculate value of the function
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd {
      return MAscoli::getProjectedReducedDeck(J, M, L,
                                              POW2(iso_mass), S, R,
                                              e*e, t,
                                              mtRsq,
                                              stot,
                                              mAsq, mBsq, mDsq,
                                              POW2(PI_MASS));
    });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  myfile.open("/tmp/ascoli.F2PiD.misha.txt");
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
  c1.SaveAs("/tmp/test_ascoli_shape.pdf");

  /************************************************************************/
  /* pi pi D - wave */
  iso_mass = 0.6;
  S = 0;
  L = 2;
  // fill X array
  ranges.first  = iso_mass + PI_MASS;
  ranges.second = 2.5;
  for (uint i=0; i < sVls.size(); i++) sVls[i] = ranges.first + (ranges.second-ranges.first)/(Npoints-1)*i;
  // calculate value of the function
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd {
      return MAscoli::getProjectedReducedDeck(J, M, L,
                                              POW2(iso_mass), S, R,
                                              e*e, t,
                                              mtRsq,
                                              stot,
                                              mAsq, mBsq, mDsq,
                                              POW2(PI_MASS));
    });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  myfile.open("/tmp/ascoli.PiPiSPiD.misha.txt");
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
  c1.SaveAs("/tmp/test_ascoli_shape.pdf)");

  return 0;
}
