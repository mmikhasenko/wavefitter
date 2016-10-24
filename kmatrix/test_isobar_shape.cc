// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <vector>
#include <fstream>
// #include <>

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "mstructures.h"

int main(int argc, char *argv[]) {
  // stable
  MIsobar rho_iso(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
  MIsobar  f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.);
  MIsobarPiPiS  pipiS_iso;

  /* generate range */
  const uint Npoints = 100;
  std::vector<double> sVls(Npoints);
  std::pair<double, double> ranges(2*PI_MASS, 2.2);
  for (uint i=0; i < sVls.size(); i++) sVls[i] = ranges.first + (ranges.second-ranges.first)/(Npoints-1)*i;

  /* data arrays */
  std::vector<cd> amp(Npoints);
  std::vector<double> reA(Npoints);
  std::vector<double> imA(Npoints);

  /* fill rho amplitude */
  /* rho */
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd { /*std::cout << e << " ";*/ return rho_iso.T(e*e); });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  std::ofstream myfile;
  myfile.open("/tmp/rho.misha.txt");
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
  c1.SaveAs("/tmp/test_isobar_shape.pdf(");

  /* f2 */
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd { /*std::cout << e << " ";*/ return f2_iso.T(e*e); });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  myfile.open("/tmp/f2.misha.txt");
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
  c1.SaveAs("/tmp/test_isobar_shape.pdf");

  /* pipiS */
  std::transform(sVls.begin(), sVls.end(), amp.begin(), [&](double e)->cd { /*std::cout << e << " ";*/ return pipiS_iso.T(e*e); });
  std::transform(amp.begin(), amp.end(), reA.begin(), [](cd A)->double { return real(A); });
  std::transform(amp.begin(), amp.end(), imA.begin(), [](cd A)->double { return imag(A); });

  myfile.open("/tmp/pipiS.misha.txt");
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
  c1.SaveAs("/tmp/test_isobar_shape.pdf)");

  return 0;
}
// =======
// #include "TH2D.h"
// #include "TGraph.h"
// #include "TCanvas.h"
// #include "MIsobar.h"
// #include "MIsobarChannel.h"
// #include "MIsobarPiPiS.h"
// #include "MTwoBodyChannel.h"
// #include "mstructures.h"
// 
// int main(int argc, char *argv[]) {
//   
//   MIsobarPiPiS  pipiS_iso;
//   MIsobar rho_iso(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
//   MIsobar  f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.);
// 
//   TCanvas c1("c1");
//   TMultiGraph *m;  
//   /* RHO */
//   m = combine(
//           SET3(
//                draw([&](double m)->double{return real(pipiS_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Real part"),
//                SetFillColor(0),
//                SetLineColor(kBlue)),
//           SET3(
//                draw([&](double m)->double{return imag(pipiS_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Imag part"),
//                SetFillColor(0),
//                SetLineColor(kOrange) ),
//           SET3(
//                draw([&](double m)->double{return abs(pipiS_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Abs part"),
//                SetFillColor(0),
//                SetLineColor(kGreen) ) );
//   m->SetTitle("Amplitude pipiS(s)");
//   m->Draw("al");
//   c1.BuildLegend();
//   c1.SaveAs("/tmp/test_isobar_shale.pdf(");
//   c1.Clear();
// 
//   /* RHO */
//   m = combine(
//           SET3(
//                draw([&](double m)->double{return real(rho_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Real part"),
//                SetFillColor(0),
//                SetLineColor(kBlue)),
//           SET3(
//                draw([&](double m)->double{return imag(rho_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Imag part"),
//                SetFillColor(0),
//                SetLineColor(kOrange) ),
//           SET3(
//                draw([&](double m)->double{return abs(rho_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Abs part"),
//                SetFillColor(0),
//                SetLineColor(kGreen) ) );
//   m->SetTitle("Amplitude rho(s)");
//   m->Draw("al");
//   c1.BuildLegend();
//   c1.SaveAs("/tmp/test_isobar_shale.pdf");
//   c1.Clear();
// 
//   /* F2 */
//   m = combine(
//           SET3(
//                draw([&](double m)->double{return real(f2_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Real part"),
//                SetFillColor(0),
//                SetLineColor(kBlue)),
//           SET3(
//                draw([&](double m)->double{return imag(f2_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Imag part"),
//                SetFillColor(0),
//                SetLineColor(kOrange) ),
//           SET3(
//                draw([&](double m)->double{return abs(f2_iso.T(m*m));}, 2*PI_MASS, 2.2),
//                SetTitle("Abs part"),
//                SetFillColor(0),
//                SetLineColor(kGreen) ) );
//   m->SetTitle("Amplitude f2(s)");
//   m->Draw("al");
//   c1.BuildLegend();
//   c1.SaveAs("/tmp/test_isobar_shale.pdf)");
//   c1.Clear();
// 
//   return 0;
// }
