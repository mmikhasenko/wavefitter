// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "MIsobar.h"
#include "MTwoBodyChannel.h"
#include "MIsobarChannel.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.h"

#include "MatrixInverse.h"

#include "TH2D.h"

int main(int argc, char *argv[]) {
  // stable
  MIsobar rho(0.77, 0.15, 0.14, 0.14, 1, 5.);
  MIsobar  f2(1.23, 0.2,  0.14, 0.14, 2, 5.);
  // rho-pi
  // MTwoBodyChannel rho_pi(rho.GetM(), 0.14);
  MIsobarChannel rho_pi(rho, 0.14);
  rho_pi.makeLookupTable(rho_pi.sth(), 5., 100);
  rho_pi.makeDisperseLookupTable(0.1, 10., 110);
  // f2-pi
  // MTwoBodyChannel  f2_pi(f2.GetM(),  0.14);
  MIsobarChannel  f2_pi(f2, 0.14);
  f2_pi.makeLookupTable(f2_pi.sth(), 5., 100);
  f2_pi.makeDisperseLookupTable(0.1, 10., 110);
  // together
  std::vector<MChannel*> channels = {&rho_pi, &f2_pi};
  MmatrixK km(channels, 2);
  // masses
  MParKeeper::gI()->set("m0", 1.402);
  MParKeeper::gI()->set("m1", 1.702);
  // couplings
  MParKeeper::gI()->set("g0", 4.5);
  MParKeeper::gI()->set("g1", 0.5);
  MParKeeper::gI()->set("h0", 1.5);
  MParKeeper::gI()->set("h1", 3.5);
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

  for (uint i=1; i <= (uint)hssh.GetNbinsX(); i++)
    for (uint j=1; j <= (uint)hssh.GetNbinsY(); j++) {
      cd s(hssh.GetXaxis()->GetBinCenter(i),
           hssh.GetYaxis()->GetBinCenter(j));
      complex_s_values[(i-1)*hssh.GetNbinsY()+(j-1)] = s;
    }

  // function to fill histogram
  std::function<void(const std::vector<double>&)> fill_hist =
    [&](const std::vector<double> &arr)->void{
    for (uint i=1; i <= (uint)hssh.GetNbinsX(); i++)
      for (uint j=1; j <= (uint)hssh.GetNbinsY(); j++)
        hssh.SetBinContent(i, j, arr[(i-1)*hssh.GetNbinsY()+(j-1)]);
  };

  /*****************************************************/

  std::cout << "Second part\n";
  //  TCanvas c2("c2", "title");
  const uint Nch = channels.size();
  c1.Divide(Nch, Nch);
  for (uint i=0; i < Nch; i++)
    for (uint j=0; j < Nch; j++) {
      if (i > j) continue;
      c1.cd(i*Nch+j+1);
      TMultiGraph *m =
        combine(
                SET1(
                     draw([&](double s)->double{return real(km.getValue(s)(i, j));}, 0.2, 4., 300),
                     SetLineColor(kBlack) ),
                SET1(
                     draw([&](double s)->double{return imag(km.getValue(s)(i, j));}, 0.2, 4., 300),
                     SetLineColor(kRed) ),
                SET1(
                     draw([&](double s)->double{return  abs(km.getValue(s)(i, j));}, 0.2, 4., 300),
                     SetLineColor(kGreen) ) );
      m->SetTitle(TString::Format("Re, Im, Abs of T_{%d%d};s(GeV)", i, j));
      m->Draw("al");
    }
  c1.SaveAs("/tmp/T.physical.test_Secondsheet.pdf");

  c1.Divide(1);
  c1.Clear();
  /************ Rho L tilde *************/
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s) { return rho_pi.rho(s); });
  std::transform(cd_values.begin(), cd_values.end(),
                 re_values.begin(), [&](cd v) { return real(v); });
  std::transform(cd_values.begin(), cd_values.end(),
                 im_values.begin(), [&](cd v) { return imag(v); });
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return abs(v); });

  // Real part
  hssh.SetTitle("real part of #rho(s)");
  fill_hist(re_values);
  hssh.Draw("lego2");
  c1.SaveAs("/tmp/real.rholtilde.test_Secondsheet.pdf");

  // Imaginary part
  hssh.SetTitle("imag part of #rho(s)");
  fill_hist(im_values);
  hssh.Draw("lego2");
  c1.SaveAs("/tmp/imag.rholtilde.test_Secondsheet.pdf");

  /************ km.getValue(0,0) *************/
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{ return det_fast(km.getValue(s));});
  std::transform(cd_values.begin(), cd_values.end(),
                 re_values.begin(), [&](cd v)->double{return real(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 im_values.begin(), [&](cd v)->double{return imag(v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return abs(1./v); });

  // Real part
  hssh.SetTitle("real part of Det[T_{I}]");
  fill_hist(re_values);
  hssh.Draw("colz");
  c1.SaveAs("/tmp/real.T.test_Secondsheet.pdf");

  // Imaginary part
  hssh.SetTitle("imag part of Det[T_{I}]");
  fill_hist(im_values);
  hssh.Draw("colz");
  c1.SaveAs("/tmp/imag.T.test_Secondsheet.pdf");

  // Abs part
  hssh.SetTitle("abs part of Det[T_{I}^{-1}]");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.SaveAs("/tmp/abs.T1_inv.test_Secondsheet.pdf");

  /************ km.getSSInverseValue(0,0) *************/
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{ return det_fast(km.getSSInverseValue(s));});
  std::transform(cd_values.begin(), cd_values.end(),
                 re_values.begin(), [&](cd v)->double{return real(1./v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 im_values.begin(), [&](cd v)->double{return imag(1./v);});
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return abs(v); });

  // Real part
  hssh.SetTitle("real part of Det[T_{II}]");
  fill_hist(re_values);
  hssh.Draw("colz");
  c1.SaveAs("/tmp/real.T2.test_Secondsheet.pdf");

  // Imaginary part
  hssh.SetTitle("imag part of Det[T_{II}]");
  fill_hist(im_values);
  hssh.Draw("colz");
  c1.SaveAs("/tmp/imag.T2.test_Secondsheet.pdf");

  // Abs part
  hssh.SetTitle("abs part of Det[T_{II}^{-1}]");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.SaveAs("/tmp/abs.T2_inv.test_Secondsheet.pdf");

  /*****************************************************/
  // Combine two pictures;
  std::transform(complex_s_values.begin(), complex_s_values.end(),
                 cd_values.begin(), [&](cd s)->cd{
                   if (imag(s) > 0) { return 1./det_fast(km.getValue(s));
                   } else { return det_fast(km.getSSInverseValue(s));} } );
  std::transform(cd_values.begin(), cd_values.end(),
                 abs_values.begin(), [&](cd v) { return abs(v); });

  // Abs part
  hssh.SetTitle("merge I,II-sheets: abs part of Det[T^{-1}]");
  fill_hist(abs_values);
  hssh.Draw("colz");
  hssh.Draw("cont3 same");
  c1.SaveAs("/tmp/abs.T1-2_inv.test_Secondsheet.pdf");

  return 1;
}
