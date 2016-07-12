// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "MIsobar.h"
#include "MTwoBodyChannel.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.h"

#include "TH2D.h"


int main(int argc, char *argv[]) {
  // stable
  MTwoBodyChannel rho_pi(0.77, 0.14);
  rho_pi.makeDisperseLookupTable(0.1, 10., 110);
  std::vector<MChannel*> channels = {&rho_pi};
  MmatrixK km(channels, 1);
  // masses
  MParKeeper::gI()->set("m1", 1.402);
  // MParKeeper::gI()->set("m2", 1.601);
  // couplings
  MParKeeper::gI()->set("g1", 2.5);
  // MParKeeper::gI()->set("h1", 4.);
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
  c1.SaveAs("/tmp/11.pdf");

  /************ Rho Isobar **************/

  combine(
          SET1(
               draw([&](double s)->double{return real(rho_pi.rho(cd(s, 1e-1)));}, 0.2, 4., 300),
               SetLineColor(kBlack) ),
          SET1(
               draw([&](double s)->double{return imag(rho_pi.rho(cd(s, 1e-1)));}, 0.2, 4., 300),
               SetLineColor(kRed) ) )->Draw("al");
  c1.SaveAs("/tmp/rho.SS.pdf");

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
  c1.SaveAs("/tmp/real.rholtilde.test_Secondsheet.pdf");

  // Imaginary part
  hssh.SetTitle("imag part of #rholtilde");
  fill_hist(im_values);
  hssh.Draw("lego2");
  c1.SaveAs("/tmp/imag.rholtilde.test_Secondsheet.pdf");

  /************ km.getValue(0,0) *************/
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{ return km.getValue(s)(0, 0);});
  std::transform(cd_values.begin(), cd_values.end(),
                 re_values.begin(), [&](cd v)->double{return real(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 im_values.begin(), [&](cd v)->double{return imag(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return abs(1./v); });

  // Real part
  hssh.SetTitle("real part of T_{I}");
  fill_hist(re_values);
  hssh.Draw("colz");
  c1.SaveAs("/tmp/real.T.test_Secondsheet.pdf");

  // Imaginary part
  hssh.SetTitle("imag part of T_{I}");
  fill_hist(im_values);
  hssh.Draw("colz");
  c1.SaveAs("/tmp/imag.T.test_Secondsheet.pdf");

  // Abs part
  hssh.SetTitle("abs part of T_{I}^{-1}");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.SaveAs("/tmp/abs.T1_inv.test_Secondsheet.pdf");

  /************ km.getSSInverseValue(0,0) *************/
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{ return km.getSSInverseValue(s)(0, 0);});
  std::transform(cd_values.begin(), cd_values.end(),
                 re_values.begin(), [&](cd v)->double{return real(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 im_values.begin(), [&](cd v)->double{return imag(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return abs(v); });

  // Real part
  hssh.SetTitle("real part of T_{II}");
  fill_hist(re_values);
  hssh.Draw("colz");
  c1.SaveAs("/tmp/real.T2_inv.test_Secondsheet.pdf");

  // Imaginary part
  hssh.SetTitle("imag part of T_{II}");
  fill_hist(im_values);
  hssh.Draw("colz");
  c1.SaveAs("/tmp/imag.T2_inv.test_Secondsheet.pdf");

  // Abs part
  hssh.SetTitle("abs part of T_{II}^{-1}");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.SaveAs("/tmp/abs.T2_inv.test_Secondsheet.pdf");

  /*****************************************************/
  // Combine two pictures;
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{
                   if (imag(s) > 0) { return 1./km.getValue(s)(0, 0);
                   } else { return km.getSSInverseValue(s)(0, 0);} } );
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return abs(v); });

  // Abs part
  hssh.SetTitle("abs part of T_{II}^{-1}");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.SaveAs("/tmp/abs.T1-2_inv.test_Secondsheet.pdf");

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
  c1.SaveAs("/tmp/T.physical.test_Secondsheet.pdf");

  return 1;
}
