// Copyright [2016] Mikhail Mikhasenko

#include <functional>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TGraph.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH1D.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "mstructures.h"

TH1D *draw(const std::vector<double> &v, double ll = 1., double rl = -1);
TH1D *draw(const std::vector<double> &v, double ll, double rl) {
  double max = (ll < rl) ? rl : *std::max_element(v.begin(), v.end());  // , [](double v1, double v2)->bool{return (v1>v2);}
  double min = (ll < rl) ? ll : *std::min_element(v.begin(), v.end());
  TH1D *h = new TH1D("h1","title", 100, min, max);
  std::for_each(v.begin(), v.end(), [&](double x)->void{h->Fill(x);});
  return h;
}

int main(int argc, char *argv[]) {

  const uint Np = 5;
  double e[Np] = {2.1, 2.2, 2.3, 2.4, 2.5};

  /********************************* M O D E L **********************************/
  double m = 2.25;
  double G = 0.2;
  double c = 1.4;
  TF1 fphi("fphi",
           [&](double* x, double* par)->double{
             double eX = x[0];
             return atan2(m*G, m*m-eX*eX);
           }, e[0], e[Np-1], 0);
  TF1 famp("famp",
           [&](double* x, double* par)->double{
             double eX = x[0];
             return c/sqrt(m*G*m*G+(m*m-eX*eX)*(m*m-eX*eX));
           }, e[0], e[Np-1], 0);
  TF1 frev("frev",
           [&](double* x, double* par)->double{
             double eX = x[0];
             return c*(m*m-eX*eX)/(m*G*m*G+(m*m-eX*eX)*(m*m-eX*eX));
           }, e[0], e[Np-1], 0);
  TF1 fimv("fimv",
           [&](double* x, double* par)->double{
             double eX = x[0];
             return c*m*G/(m*G*m*G+(m*m-eX*eX)*(m*m-eX*eX));
           }, e[0], e[Np-1], 0);
  std::cout << "Test how well the functions are defined!" << std::endl;
  /******************************************************************************/
  
  /******************************* D A T A **************************************/
  // phase
  double phi[Np]; for (uint i=0; i < Np; i++) phi[i] = fphi.EvalPar(&e[i]);
  double dphi[Np] = {0.05, 0.2, 0.2, 0.2, 0.05};
  // amplitude
  double A[Np]; for (uint i=0; i < Np; i++) A[i] = famp.EvalPar(&e[i]);
  double dA[Np] = {0.1, 0.5, 0.5, 0.1, 0.1};

  /***************************************/
  gRandom->SetSeed(0);
  // phase
  TGraphErrors gphi(Np, e, phi, 0, dphi); gphi.SetTitle("Phase");
  for (int i = 0; i < gphi.GetN(); i++) { gphi.GetY()[i] = phi[i]; gphi.GetEY()[i] = dphi[i]; }
  // mod
  TGraphErrors gamp(Np, e, A, 0, dA);  gamp.SetTitle("absolute value");
  for (int i = 0; i < gphi.GetN(); i++) { gamp.GetY()[i] = A[i]; gamp.GetEY()[i] = dA[i]; }

  // real
  TGraphErrors grev(Np, e, A, 0, dA);  grev.SetTitle("real value");
  for (int i = 0; i < grev.GetN(); i++) {
    grev.GetY()[i] = gamp.GetY()[i]*cos(gphi.GetY()[i]);
    grev.GetEY()[i] = sqrt(pow(gamp.GetEY()[i]*cos(gphi.GetY()[i]), 2) +
                           pow(gamp.GetY()[i]*sin(gphi.GetY()[i])*gphi.GetEY()[i], 2));
  }
  // imag
  TGraphErrors gimv(Np, e, A, 0, dA);  gimv.SetTitle("imag value");
  for (int i = 0; i < gimv.GetN(); i++) {
    gimv.GetY()[i] = gamp.GetY()[i]*sin(gphi.GetY()[i]);
    gimv.GetEY()[i] = sqrt(pow(gamp.GetEY()[i]*sin(gphi.GetY()[i]), 2) +
                           pow(gamp.GetY()[i]*cos(gphi.GetY()[i])*gphi.GetEY()[i], 2));
  }
  // covariance
  TGraph gcov(Np, e, A);
  for (int i = 0; i < gimv.GetN(); i++) {
    double ph = gphi.GetY()[i];
    gcov.GetY()[i] = cos(ph)*sin(ph)*(pow(gamp.GetEY()[i], 2) +
                                      pow(gamp.GetY()[i]*gphi.GetEY()[i], 2));
  }
  
  /****************************************************************************************/
  const uint Nmethods = 3;
  std::function<double(const double*)> fchi2[Nmethods];

  fchi2[0] = [&](const double *pars)->double {
      double chi2 = 0;
      m = pars[0];
      G = pars[1];
      c = pars[2];
      for (int i = 0; i < gphi.GetN(); i++) chi2 += pow((grev.GetY()[i]-frev.EvalPar(grev.GetX()+i)), 2); // 
      for (int i = 0; i < gamp.GetN(); i++) chi2 += pow((gimv.GetY()[i]-fimv.EvalPar(gimv.GetX()+i)), 2); // 
      // std::cout << "chi2 = " << chi2 << "\n";
      return chi2;
  };
  fchi2[1] = [&](const double *pars)->double {
      double chi2 = 0;
      m = pars[0];
      G = pars[1];
      c = pars[2];
      for (int i = 0; i < gphi.GetN(); i++) chi2 += pow((grev.GetY()[i]-frev.EvalPar(grev.GetX()+i))/grev.GetEY()[i], 2); // 
      for (int i = 0; i < gamp.GetN(); i++) chi2 += pow((gimv.GetY()[i]-fimv.EvalPar(gimv.GetX()+i))/gimv.GetEY()[i], 2); // 
      // std::cout << "chi2 = " << chi2 << "\n";
      return chi2;
  };
  fchi2[2] = [&](const double *pars)->double {
      double chi2 = 0;
      m = pars[0];
      G = pars[1];
      c = pars[2];
      for (int i = 0; i < gphi.GetN(); i++) {
        double sigmaX = grev.GetEY()[i];
        double sigmaY = gimv.GetEY()[i];
        double covXY = 0.8*sigmaX*sigmaY; // gcov.GetY()[i];  // 
        double deltaX = grev.GetY()[i]-frev.EvalPar(grev.GetX()+i);
        double deltaY = gimv.GetY()[i]-fimv.EvalPar(gimv.GetX()+i);
        chi2 += (deltaX*deltaX*sigmaY*sigmaY+
                 deltaY*deltaY*sigmaX*sigmaX-
                 2*deltaX*deltaY*covXY) / (sigmaX*sigmaX*sigmaY*sigmaY-covXY*covXY);
        // chi2 += pow(deltaX/sigmaX,2);
        // chi2 += pow(deltaY/sigmaY,2);
      }
      return chi2;
  };
  /**************************************MINIMIZE******************************************/
  ROOT::Math::Minimizer* min[Nmethods];
  ROOT::Math::Functor *functor[Nmethods];
  for (uint met = 0; met < Nmethods; met++) {
    // Build minimizer
    min[met] =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    
    // set tolerance , etc...
    min[met]->SetMaxFunctionCalls(100000);
    min[met]->SetTolerance(0.001);
    min[met]->SetStrategy(1);
    min[met]->SetPrintLevel(3);
    min[met]->Options().Print();
  
    // Create funciton wrapper for minmizer a IMultiGenFunction type
    functor[met] = new ROOT::Math::Functor(fchi2[met], 3);
    min[met]->SetFunction(*functor[met]);
    min[met]->SetVariable(0, "m", m, 0.1);
    min[met]->SetVariable(1, "G", G, 0.1);
    min[met]->SetVariable(2, "c", c, 0.1);
  };
    
  TCanvas c1("can");
  c1.Divide(2, 2);
  const double Natt = 10000;
  std::vector<double> pw[Nmethods][3];
  for (uint e = 0; e < Natt; e++) {
    // redistribute data
    for (int i = 0; i < gphi.GetN(); i++) gphi.GetY()[i] = gRandom->Gaus(phi[i], dphi[i]);
    for (int i = 0; i < gphi.GetN(); i++) gamp.GetY()[i] = gRandom->Gaus(A[i], dA[i]);
    for (int i = 0; i < grev.GetN(); i++) grev.GetY()[i] = gamp.GetY()[i]*cos(gphi.GetY()[i]);
    for (int i = 0; i < gimv.GetN(); i++) gimv.GetY()[i] = gamp.GetY()[i]*sin(gphi.GetY()[i]);
    // minimize
    for (uint met = 0; met < Nmethods; met++) {
      m = 2.5; G = 0.03; c = 0.5;
      min[met]->Minimize();
      // plot
      if (m<0 || m>5) continue;
      if (e == 53) {
        c1.cd(1); gphi.Draw("ap"); fphi.Draw("same");
        c1.cd(2); gamp.Draw("ap"); famp.Draw("same");
        c1.cd(3); grev.Draw("ap"); frev.Draw("same");
        c1.cd(4); gimv.Draw("ap"); fimv.Draw("same");
        if  (met == 0) c1.Print("/tmp/test_Fit_ReIm.pdf(");
        else c1.Print("/tmp/test_Fit_ReIm.pdf");
      }
      pw[met][0].push_back(m);
      pw[met][1].push_back(G);
      pw[met][2].push_back(c);
    }
  }
  c1.Clear(); c1.Divide(3, 3);
  TH1D *hw[Nmethods];
  for (uint i=0; i < 3; i++) {
    c1.cd(i+1);
    hw[0]  = draw(pw[0][i]); hw[0]->SetFillColor(kGray);   hw[0]->SetTitle(TString::Format("par[%d] distribution", i)); hw[0] ->Draw();
    c1.cd(3+i+1); // ->SetLogy();
    hw[1]  = draw(pw[1][i], hw[0]->GetXaxis()->GetXmin(), hw[0]->GetXaxis()->GetXmax()); hw[1]->SetFillColorAlpha(kOrange, 0.5); hw[1]->Draw("same");
    c1.cd(6+i+1);  // ->SetLogy();
    hw[2]  = draw(pw[2][i], hw[0]->GetXaxis()->GetXmin(), hw[0]->GetXaxis()->GetXmax()); hw[2]->SetFillColorAlpha(kRed, 0.5);    hw[2]->Draw("same");
  }
  c1.Print("/tmp/test_Fit_ReIm.pdf)");
  // c1.SaveAs("/tmp/b.test_Errors.pdf");
  
  return 0;
}
