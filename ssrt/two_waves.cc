//Misha Mikasenko
#include <iostream>
#include <sstream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>

#include <TH1D.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TVirtualFitter.h>
#include <Math/MinimizerOptions.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <constants.h>

#define NPAR_FIT1  2
#define APAR_FIT1  1
#define NPAR_FIT2  2
#define APAR_FIT2  1

#define S1 0.5
#define S2 0.5

using namespace std; 
typedef std::complex<double> cd;

cd BreitWignerA(double x,double m, double G0);
cd BreitWignerAndDeckA(double x,double m, double G0,double a, double b);
double intensity1(double *x, double *par);
double intensity2(double *x, double *par);
double dphase(double *x, double *par);

cd expand(cd s, cd s1, vector<cd> par);
cd vtable(double s,vector<pair<double,cd> > &arr);
cd Afull(cd s, vector<cd> apar, vector<cd> npar, cd s1, cd s2);
double globChi2(const double *par);
void prepareNames(int npar, int apar, string *var_name,const char *pref="");

void fillvectors(const double *par, vector<cd> &apar, vector<cd> &npar);
double normPhase(double a) {
  if(a>M_PI)  return a-2*M_PI*int( a/M_PI);
  if(a<-M_PI) return a+2*M_PI*int(-a/M_PI);
  return a;
}


vector<vector<pair<double,cd> > > *garr;
TH1D *gh1, *gh2, *gh3;

int main() {

  //**************** Load integrals ****************//
  const int nTable = 3;
  string table_name[nTable] = {"f0.txt","f1.txt","f2.txt"};
  
  TGraph2D *table[nTable];
  for(int i=0;i<nTable;i++) {
    table[i] = new TGraph2D(table_name[i].c_str());
    cout << "found first file: " << table_name[i].c_str()
	 << ", N = " << table[i]->GetN()
	 << endl;
  }
  vector<vector<pair<double,cd> > > data;
  for(int i=0;i<nTable;i++) {
    vector<pair<double,cd> > vj;
    for(int j=0;j<table[i]->GetN();j++) {
      cd value(table[i]->GetY()[j],
	       table[i]->GetZ()[j]);
      vj.push_back(make_pair<double,cd>(table[i]->GetX()[j],
					value));
    }
    data.push_back(vj);
  }
  garr=&data;
  cout << "test access " << (*garr)[0][5].second << endl; 
  
  //**************** Load the data *****************//
  double dataLeft = 1;  double dataRight = 2; double dataNbin = 100;
  TH1D *hist_data_int1  = new TH1D("data1","Intensity;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  TH1D *hist_data_int2  = new TH1D("data2","Intensity;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  TH1D *hist_data_dphi = new TH1D("phase","Phase;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  //first channel
  for(int i=1;i<=hist_data_int1->GetNbinsX();i++) {
    hist_data_int1->SetBinContent(i,
				  norm(BreitWignerAndDeckA(hist_data_int1->GetBinCenter(i),
							  R1_MASS,R1_WIDTH,1,M_PI/3)));
    hist_data_int1->SetBinError(i,0.1*sqrt(hist_data_int1->GetBinContent(i)));
  }
  //second channel
  for(int i=1;i<=hist_data_int2->GetNbinsX();i++) {
    hist_data_int2->SetBinContent(i,
				  norm(BreitWignerAndDeckA(hist_data_int2->GetBinCenter(i),
							  R1_MASS,R1_WIDTH,2,-M_PI/2)));
    hist_data_int2->SetBinError(i,0.1*sqrt(hist_data_int2->GetBinContent(i)));
  }
  //phase difference
  for(int i=1;i<=hist_data_dphi->GetNbinsX();i++) {
    double s = hist_data_dphi->GetBinCenter(i);
    cd res1 = BreitWignerAndDeckA(s,R1_MASS,R1_WIDTH,1,M_PI/3);
    cd res2 = BreitWignerAndDeckA(s,R1_MASS,R1_WIDTH,2,-M_PI/2);
    hist_data_dphi->SetBinContent(i,arg(res1*conj(res2)));
    hist_data_dphi->SetBinError(i,gRandom->Rndm()*M_PI/15.);
  }

  gh1 = hist_data_int1;
  gh2 = hist_data_int2;
  gh3 = hist_data_dphi;

  TCanvas can("c1");
  can.Divide(2,2);
 
  can.cd(1); gh1->Draw();
  can.cd(4); gh2->Draw();
  can.cd(2); gh3->Draw();

  //***************** Fit the data *****************//
  //Build minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  //min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  
  const int NumPar = 
    2*NPAR_FIT1 + 2*APAR_FIT1 +
    2*NPAR_FIT2 + 2*APAR_FIT2;
  
  // create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor f(&globChi2,NumPar); 
  min->SetFunction(f);
  
  // starting point
  double start_par[NumPar];
  gRandom->SetSeed(0);
  for(int i=0;i<NumPar;i++) start_par[i] = gRandom->Rndm(); 
  
  double step = 0.01;
  //form names
  string var_name[NumPar];
  prepareNames(NPAR_FIT1,APAR_FIT1,&var_name[0],"ch1_");
  prepareNames(NPAR_FIT2,APAR_FIT2,&var_name[2*(NPAR_FIT1+APAR_FIT1)],"ch2_");
  // Set the free variables to be minimized!
  for(int i=0;i<NumPar;i++) {//var_name[i].c_str()
    cout << "name = " << var_name[i] << endl;
    min->SetVariable(i,var_name[i],start_par[i], step);
  }
  // do the minimization
  min->Minimize(); 
  const double *final_pars = min->X();

  TF1 *fInt1 = new TF1("int1",intensity1,1.2,1.8,2*(NPAR_FIT1+APAR_FIT1)); fInt1->SetParameters(final_pars);
  TF1 *fInt2 = new TF1("int2",intensity2,1.2,1.8,2*(NPAR_FIT2+APAR_FIT2)); fInt2->SetParameters(&final_pars[2*(NPAR_FIT1+APAR_FIT1)]);
  TF1 *fPhi  = new TF1("dphi",dphase,    1.2,1.8,NumPar);              fPhi->SetParameters(final_pars);

  can.cd(1); fInt1->Draw("same");
  can.cd(2); fPhi ->Draw("same");
  can.cd(4); fInt2->Draw("same");
  
  can.SaveAs("c1.png");
  
  cout << "finished" << endl;
  for(int i=0;i<nTable;i++) delete table[i];
  return 0;
}


////////////////////////////////////////////////////////////////////////
//                     /////////////////////////////////////////////////
////////////////////                   ///////////////            //////
/////////////////////////////////////                    ///////////////
////////////////////////////////////////////////////////////////////////

double globChi2(const double *par) {

  vector<cd> npar1(NPAR_FIT1), apar1(APAR_FIT1); fillvectors(&par[0],apar1,npar1);
  vector<cd> npar2(NPAR_FIT2), apar2(APAR_FIT2); fillvectors(&par[2*(NPAR_FIT1+APAR_FIT1)],apar2,npar2);

  // *********************** Calculate chi2 ********************** //
  double chi2 = 0;
  //Intensity1
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double s = gh1->GetBinCenter(i);
    cd amp = Afull(s,apar1,npar1,S1,S2);
    double intens = norm(amp);
    double err = gh1->GetBinError(i);
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //Intensity2
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double s = gh2->GetBinCenter(i);
    cd amp = Afull(s,apar2,npar2,S1,S2);
    double intens = norm(amp);
    double err = gh2->GetBinError(i);
    chi2 += pow(intens - gh2->GetBinContent(i),2)/(err*err);
  }
  //dPhase
  for(int i=1;i<=gh3->GetNbinsX();i++) {
    double s = gh3->GetBinCenter(i);
    cd amp1 = Afull(s,apar1,npar1,S1,S2);
    cd amp2 = Afull(s,apar2,npar2,S1,S2);
    double phi = arg(amp1*conj(amp2));
    double err = gh3->GetBinError(i);
    double dPhi = normPhase(phi-gh3->GetBinContent(i));
    chi2 += pow(dPhi,2)/(err*err);
  }
  return chi2;
}

void fillvectors(const double *par, vector<cd> &apar, vector<cd> &npar) {

  //**************** Convert parameters to format **************** //
  //Constract parameters N
  for(int i=0;i<npar.size();i++) {
    cd ci(par[2*i],par[2*i+1]); 
    npar[i] = ci;
  }
  //Constract Alpha
  for(int i=0;i<apar.size();i++) {
    cd normi(par[2*npar.size()+2*i],par[2*npar.size()+2*i+1]);  
    apar[i] = normi;
  }

}


void prepareNames(int npar, int apar, string *var_name,const char *pref) {
  for(int i=0;i<npar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rN"<<i; iIO<<pref<<"iN"<<i; var_name[2*i] = rIO.str(); var_name[2*i+1] = iIO.str();}
  for(int i=0;i<apar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rA"<<i; iIO<<pref<<"iA"<<i; var_name[2*npar+2*i] = rIO.str(); var_name[2*npar+2*i+1] = iIO.str();}
}

double intensity1(double *x, double *par) {
  vector<cd> npar(NPAR_FIT1), apar(APAR_FIT1);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0],apar,npar,S1,S2);
  return norm(amp);
}
double intensity2(double *x, double *par) {
  vector<cd> npar(NPAR_FIT2), apar(APAR_FIT2);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0],apar,npar,S1,S2);
  return norm(amp);
}

double dphase(double *x, double *par) {
  vector<cd> npar1(NPAR_FIT1), apar1(APAR_FIT1);
  vector<cd> npar2(NPAR_FIT2), apar2(APAR_FIT2);
  fillvectors(&par[0],apar1,npar1);
  fillvectors(&par[2*(NPAR_FIT1+APAR_FIT1)],apar2,npar2);
  cd amp1 = Afull(x[0],apar1,npar1,S1,S2);
  cd amp2 = Afull(x[0],apar2,npar2,S1,S2);
  return arg(amp1*conj(amp2));
}

cd BreitWignerA(double x,double m, double G0) {
  cd unit(0,1); 
  double p  = sqrt(lambda(x,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*sqrt(x));
  double p0 = sqrt(lambda(m*m,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*m);
  double R = 5; double BW = (1+pow(R*p0,2))/(1+pow(R*p,2));
  double G = G0*p/p0* BW;
  cd A = sqrt(m*G)/(m*m-x-unit*m*G);
  return A;
} 

cd BreitWignerAndDeckA(double x,double m, double G0,double a, double b) {
  cd unit(0,1); 
  double p  = sqrt(lambda(x,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*sqrt(x));
  double p0 = sqrt(lambda(m*m,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*m);
  double R = 5; double BW = (1+pow(R*p0,2))/(1+pow(R*p,2));
  double G = G0*p/p0* BW;
  cd uni(0,1);
  cd A = sqrt(m*G)/(m*m-x-unit*m*G) + a*exp(-p/2)*exp(unit*b);
  return A;
} 

cd expand(cd s, cd s1, vector<cd> par) {
  cd res=0;
  for(int i=0;i<par.size();i++) res+=pow(w(s,s1),i)*par[i];
  return res;
}

cd vtable(double s,vector<pair<double,cd> > &arr) {
  int N = arr.size();
  double sf = arr[0].first;
  double sl = arr[N-1].first;
  double step = arr[1].first-arr[0].first;
  int Nstep = (s-sf)/step;
  return arr[Nstep].second+
    (arr[Nstep+1].second-arr[Nstep].second) / 
    (s                  -arr[Nstep].second);
}

cd Afull(cd s, vector<cd> apar, vector<cd> npar, 
	 cd s1, cd s2) {

  //Constract D
  cd dsum =0;
  for(int i=0;i<npar.size();i++) dsum += npar[i]*vtable(real(s),(*garr)[i]);
  cd Denomin = 1. - s / (2*M_PI) * dsum;
  //Constract N, Alpha and A.
  cd Numir = expand(s,s1,npar);
  cd Aprod = expand(s,s2,apar);
  cd Afull = Aprod*Numir/Denomin;

  return Afull;
}
