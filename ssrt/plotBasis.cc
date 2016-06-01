#include <iostream>
#include <complex>
#include <fstream>
#include <initializer_list>
#include <MIsobar.h>
#include <NoverD.h>
#include <constants.h>
#include <deflib.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>

using namespace std; 
typedef std::complex<double> cd;

TMultiGraph *mg(TGraph &g1,TGraph &g2);
TMultiGraph *mgm(std::initializer_list<TGraph*> a_args);

TGraph *grRe   (vector<double> &x, cd *arr);
TGraph *grIm   (vector<double> &x, cd *arr);
TGraph *grArgon(vector<double> &x, cd *arr);

cd rho(cd s);
cd rho3(cd s);
int main() {
  
  gROOT->ProcessLine(".x ~/Documents/root-scripts/cv_n.C");

  MIsobar RhoPi(F2_MASS,F2_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  const int N=100;
  double m_l=3*PI_MASS, m_r = 2;
  double step = (m_r-m_l)/(N-1);
  double x[N],y[N],yc[N];
  for(int i=0;i<N;i++) {
    x[i] = m_l+i*step;
    y [i] = real(RhoPi.rho3(x[i]*x[i]));
    yc[i] = (x[i]>F2_MASS+PI_MASS) ? real(rho(x[i]*x[i])) : 0;
  }
  
  TGraph gr (N,x,y);  gr .SetLineWidth(2); gr.SetTitle("Phase space;#sqrt(s)");    
  TGraph grc(N,x,yc); grc.SetLineWidth(2); grc.SetLineColor(kRed); 
  TCanvas c1("c1","can",0,0,500,500);
  gr.Draw("alp"); grc.Draw("same");
  c1.SaveAs("/tmp/rho.png");
  
  //NoverD
  const int Nexp = 3;
  NoverD nd1(Nexp,pow(3*PI_MASS,2)-pow(PI_MASS,2),rho, pow(F2_MASS+PI_MASS,2),1000,500);
  NoverD nd2(Nexp,pow(3*PI_MASS,2)-pow(PI_MASS,2),       rho3,pow(3*PI_MASS,2)       ,1000,500);
  nd1.npar[0]=1.; nd1.npar[1]=1.; nd1.npar[2]=1.;
  nd2.npar[0]=1.; nd2.npar[1]=1.; nd2.npar[2]=1.;

  cout << "*************************************************************************" << endl;
  cout << "************************Calculation of real axis*************************" << endl;
  const int Nm=1000;
  vector<double> m1(Nm), m2(Nm);
  cd vi1[Nexp][Nm], vi2[Nexp][Nm];
  for(int i=0;i<Nm;i++) {
    m1[i] = 1e-3+1./sqrt(nd1.GetThreshold())/(Nm-1)*i;
    m2[i] = 1e-3+1./sqrt(nd2.GetThreshold())/(Nm-1)*i;
    for(int k=0;k<Nexp;k++) {
      cd vk1 = nd1.getvi(1./(m1[i]*m1[i]),k); vi1[k][i] = vk1;
      cd vk2 = nd2.getvi(1./(m2[i]*m2[i]),k); vi2[k][i] = vk2;
    }
  }
  ///////////////////
  TCanvas c2("c2","can2",0,0,1500,Nexp*500);
  c2.Divide(3,Nexp);
  for(int k=0;k<Nexp;k++) {
    c2.cd(Nexp*k+1); grRe   (m1,vi1[k])->Draw("apl");
    c2.cd(Nexp*k+2); grIm   (m1,vi1[k])->Draw("apl");
    c2.cd(Nexp*k+3); grArgon(m1,vi1[k])->Draw("apl");
  }   
  c2.SaveAs("/tmp/integrals1.png");
  c2.SaveAs("/tmp/integrals1.pdf");
  ///////////////////
  TCanvas c3("c3","can3",0,0,1500,Nexp*500);
  c3.Divide(3,Nexp);
  for(int k=0;k<Nexp;k++) {
    c3.cd(Nexp*k+1); grRe   (m2,vi2[k])->Draw("apl");
    c3.cd(Nexp*k+2); grIm   (m2,vi2[k])->Draw("apl");
    c3.cd(Nexp*k+3); grArgon(m2,vi2[k])->Draw("apl");
  }   
  c3.SaveAs("/tmp/integrals2.png");
  c3.SaveAs("/tmp/integrals2.pdf");
  ///////////////////
  TCanvas c4("c3","can3",0,0,1500,Nexp*500);
  c4.Divide(3,Nexp);
  for(int k=0;k<Nexp;k++) {
    TGraph *gr1,*gr2;
    gr2 = grRe   (m2,vi2[k]); gr2->SetLineStyle(1); 
    gr1 = grRe   (m1,vi1[k]); gr1->SetLineStyle(2); 
    c4.cd(Nexp*k+1); mg(*gr1,*gr2)->Draw("al");
    gr2 = grIm   (m2,vi2[k]); gr2->SetLineStyle(1); 
    gr1 = grIm   (m1,vi1[k]); gr1->SetLineStyle(2); 
    c4.cd(Nexp*k+2); mg(*gr1,*gr2)->Draw("al");
    gr2 = grArgon(m2,vi2[k]); gr2->SetLineStyle(1); 
    gr1 = grArgon(m1,vi1[k]); gr1->SetLineStyle(2);
    c4.cd(Nexp*k+3); mg(*gr1,*gr2)->Draw("al");
  }   
  c4.SaveAs("/tmp/integrals12.png");
  c4.SaveAs("/tmp/integrals12.pdf");

  return 1;
}

cd rho(cd s) {
  if(real(s)>1000 && fabs(imag(s))<1e-5) return 1./(8*M_PI)*(1.0-(pow(F2_MASS,2)+pow(PI_MASS,2))/s);
  cd phspc2 = (s-pow(F2_MASS+PI_MASS,2))*(s-pow(F2_MASS-PI_MASS,2))/(s*s);
  return 1./(8*M_PI)*sqrtPi(phspc2);
}

cd rho3(cd s) {
  if(real(s)>1000 && fabs(imag(s))<1e-5) return 1./(8*M_PI)*(1.0-(pow(F2_MASS,2)+pow(PI_MASS,2))/s);
  MIsobar RhoPi(F2_MASS,F2_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  cd value = RhoPi.rho3(s);
  //cout << "s = " << s << ", -> " << value << endl;
  return value;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

TGraph *grRe(vector<double> &x, cd *arr) {
  double yi[x.size()];
  double xi[x.size()];
  for(int i=0;i<x.size();i++) {xi[i] = x[i]; yi[i] = real(arr[i]);}
  return new TGraph(x.size(),xi,yi); 
} 
TGraph *grIm(vector<double> &x, cd *arr) {
  double yi[x.size()];
  double xi[x.size()];
  for(int i=0;i<x.size();i++) {xi[i] = x[i]; yi[i] = imag(arr[i]);}
  return new TGraph(x.size(),xi,yi); 
} 
TGraph *grArgon(vector<double> &x, cd *arr) {
  double yi[x.size()];
  double xi[x.size()];
  for(int i=0;i<x.size();i++) {xi[i] = real(arr[i]); yi[i] = imag(arr[i]);}
  return new TGraph(x.size(),xi,yi); 
} 


TMultiGraph *mg(TGraph &g1,TGraph &g2) {
  TMultiGraph *_mg = new TMultiGraph();
  _mg->Add(&g1);
  _mg->Add(&g2);
  return _mg;
}

TMultiGraph *mgm(std::initializer_list<TGraph*> a_args) {
  TMultiGraph *_mg = new TMultiGraph();
  for (auto i: a_args) _mg->Add(i);
  return _mg;
}
