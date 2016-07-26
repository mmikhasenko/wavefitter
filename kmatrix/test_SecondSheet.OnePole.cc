// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MTwoBodyChannel.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.h"

#include "MatrixInverse.h"

#include "TH2D.h"


int main(int argc, char *argv[]) {
  // stable
  MTwoBodyChannel rho_pi(0.77, 0.14);
  rho_pi.makeDisperseLookupTable(0.1, 10., 110);
  std::vector<MChannel*> channels = {&rho_pi};
  MmatrixK km(channels, 1);
  // masses
  MParKeeper::gI()->set("m0", 1.402);
  // MParKeeper::gI()->set("m1", 1.601);
  // couplings
  MParKeeper::gI()->set("g0", 4.5);
  // MParKeeper::gI()->set("h0", 4.);
  MParKeeper::gI()->printAll();

  // canva
  TCanvas c1("c1");
  // hist
  TH2D hssh("ssh", "Second sheet",
                    98, 0.9, 4.0,
                    98, -1.0, 1.0);
  hssh.SetStats(0);

  std::vector<cd> complex_s_values(hssh.GetNbinsX()*hssh.GetNbinsY());
  std::vector<cd> cd_values(complex_s_values.size());
  std::vector<double> re_values(complex_s_values.size());
  std::vector<double> im_values(complex_s_values.size());
  std::vector<double> abs_values(complex_s_values.size());

  for (int i=1; i <= hssh.GetNbinsX(); i++)
    for (int j=1; j <= hssh.GetNbinsY(); j++) {
      cd s(hssh.GetXaxis()->GetBinCenter(i),
           hssh.GetYaxis()->GetBinCenter(j));
      complex_s_values[(i-1)*hssh.GetNbinsY()+(j-1)] = s;
    }

  // function to fill histogram
  std::function<void(const std::vector<double>&)> fill_hist =
    [&](const std::vector<double> &arr)->void{
    for (int i=1; i <= hssh.GetNbinsX(); i++)
      for (int j=1; j <= hssh.GetNbinsY(); j++)
        hssh.SetBinContent(i, j, arr[(i-1)*hssh.GetNbinsY()+(j-1)]);
  };

  combine(
          SET1(
               draw([&](double s)->double{return real(rho_pi.rholtilde(cd(s, 1e-5)));}, 0.2, 4., 300),
               SetLineColor(kBlack) ),
          SET1(
               draw([&](double s)->double{return imag(rho_pi.rholtilde(cd(s, 1e-5)));}, 0.2, 4., 300),
               SetLineColor(kRed) ) )->Draw("al");
  const char* pdfout = "/tmp/plots.SecondSheet.OnePole.pdf";
  c1.Print(TString::Format("%s(",pdfout));
  /************ Rho Isobar **************/

  combine(
          SET1(
               draw([&](double s)->double{return real(rho_pi.rho(cd(s, 1e-6)));}, 0.2, 4., 300),
               SetLineColor(kBlack) ),
          SET1(
               draw([&](double s)->double{return imag(rho_pi.rho(cd(s, 1e-6)));}, 0.2, 4., 300),
               SetLineColor(kRed) ) )->Draw("al");
  c1.Print(pdfout);

  /************ Rho L tilde *************/
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s) { return rho_pi.rholtilde(s); });
  std::transform(cd_values.begin(), cd_values.end(),
                 re_values.begin(), [&](cd v) { return real(v); });
  std::transform(cd_values.begin(), cd_values.end(),
                 im_values.begin(), [&](cd v) { return imag(v); });
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return abs(v); });

  // Real part
  hssh.SetTitle("real part of #tilde{#rho_{l}}");
  fill_hist(re_values);
  hssh.Draw("lego2");
  c1.Print(pdfout);

  // Imaginary part
  hssh.SetTitle("imag part of #rholtilde");
  fill_hist(im_values);
  hssh.Draw("lego2");
  c1.Print(pdfout);

  /************ First sheet *************/
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{ return 1./det_fast(km.getFSdenominator(s));});
  std::transform(cd_values.begin(), cd_values.end(),
                 re_values.begin(), [&](cd v)->double{return real(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 im_values.begin(), [&](cd v)->double{return imag(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return log10(abs(v)); });

  // Real part
  hssh.SetTitle("Real part of D_{I}^{-1}");
  fill_hist(re_values);
  hssh.Draw("colz");
  c1.Print(pdfout);

  // Imaginary part
  hssh.SetTitle("Imag part of D_{I}^{-1}");
  fill_hist(im_values);
  hssh.Draw("colz");
  c1.Print(pdfout);

  // Abs part
  hssh.SetTitle("Ln@Abs part of log @ D_{I}^{-1}");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.Print(pdfout);

  /************ Second sheet *************/
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{ return 1./det_fast(km.getSSdenominator(s));});
  std::transform(cd_values.begin(), cd_values.end(),
                 re_values.begin(), [&](cd v)->double{return real(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 im_values.begin(), [&](cd v)->double{return imag(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return log10(abs(v)); });

  // Real part
  hssh.SetTitle("Real part of D_{II}^{-1}");
  fill_hist(re_values);
  hssh.Draw("colz");
  c1.Print(pdfout);

  // Imaginary part
  hssh.SetTitle("Imag part of D_{II}^{-1}");
  fill_hist(im_values);
  hssh.Draw("colz");
  c1.Print(pdfout);

  // Abs part
  hssh.SetTitle("Ln@Abs part of D_{II}^{-1}");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.Print(pdfout);

  /*****************************************************/
  // Combine two pictures;
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{
                   if (imag(s) > 0) { return 1./det_fast(km.getFSdenominator(s));
                   } else { return 1./det_fast(km.getSSdenominator(s));} } );
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return log10(abs(v)); });

  // Abs part
  hssh.SetTitle("Ln@Abs part of D_{I}^{-1}~D_{II}^{-1}");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.Print(pdfout);

  /*****************************************************/

  combine(
          SET1(
               draw([&](double s)->double{return real(km.getValue(s)(0, 0));}, 0.2, 4., 300),
               SetLineColor(kBlack) ),
          SET1(
               draw([&](double s)->double{return imag(km.getValue(s)(0, 0));}, 0.2, 4., 300),
               SetLineColor(kRed) ),
          SET1(
               draw([&](double s)->double{return  abs(km.getValue(s)(0, 0));}, 0.2, 4., 300),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.Print(TString::Format("%s)",pdfout));

  return 1;
}
