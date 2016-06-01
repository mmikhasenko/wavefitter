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

#define NPAR_FIT  4
#define APAR_FIT  1

#define S1 0.5
#define S2 0.5

#define LEFT_INT   RHO_MASS+PI_MASS
#define RIGHT_INT  1.4

using namespace std; 
typedef std::complex<double> cd;

cd BreitWignerA(double x,double m, double G0);
double intensity(double *x, double *par);

cd expand(cd s, cd s1, vector<cd> par);
cd vtable(double s,vector<pair<double,cd> > &arr);
cd Afull(cd s, vector<cd> apar, vector<cd> npar, cd s1, cd s2);
double globChi2(const double *par);

void fillvectors(const double *par, vector<cd> &apar, vector<cd> &npar);
double intensity(double *x, double *par);
double phase(double *x, double *par);
double normPhase(double a) {
  if(a>M_PI)  return a-2*M_PI*int( a/M_PI);
  if(a<-M_PI) return a+2*M_PI*int(-a/M_PI);
  return a;
}
double rho(double e, double m1, double m2) {
  return sqrt(lambda(e*e,m1*m1,m2*m2))/(e*e);
}


vector<vector<pair<double,cd> > > *garr;
TH1D *gh1, *gh2;

int main() {

  //**************** Load integrals ****************//
  string table_name[] = {"f0.txt","f1.txt","f2.txt","f3.txt","f4.txt"};
  
  TGraph2D *table[NPAR_FIT];
  for(int i=0;i<NPAR_FIT;i++) {
    table[i] = new TGraph2D(table_name[i].c_str());
    cout << "found first file: " << table_name[i].c_str()
	 << ", N = " << table[i]->GetN()
	 << endl;
  }
  vector<vector<pair<double,cd> > > data;
  for(int i=0;i<NPAR_FIT;i++) {
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
  double dataLeft = RHO_MASS+PI_MASS;  double dataRight = 1.5; double dataNbin = 100;
  TH1D *hist_data = new TH1D("data","Intensity;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  TH1D *hist_data_phase = new TH1D("phase","Phase;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  for(int i=1;i<=hist_data->GetNbinsX();i++) {
    hist_data->SetBinContent(i,
			     gRandom->Gaus(norm(BreitWignerA(hist_data->GetBinCenter(i),
							     A1_MASS,A1_WIDTH)),
					   1)
			     );
  }
  for(int i=1;i<=hist_data->GetNbinsX();i++) {
    hist_data_phase->SetBinContent(i,
				   arg(BreitWignerA(hist_data->GetBinCenter(i),
						    A1_MASS,A1_WIDTH)));

  }

  
  gh1 = hist_data;
  gh2 = hist_data_phase;

  TCanvas can("c1");
  can.Divide(2,1);
 
  can.cd(1); gh1->Draw();
  can.cd(2); gh2->Draw();

  //***************** Fit the data *****************//
  //Build minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  //min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  
  const int NumPar = 2*NPAR_FIT + 2*APAR_FIT;
  
  // create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor f(&globChi2,NumPar); 
  min->SetFunction(f);
  
  // starting point
  double start_par[NumPar];
  gRandom->SetSeed(0);
  for(int i=0;i<NumPar;i++) start_par[i] = gRandom->Rndm(); 
  
  double step = 0.01;
  //form names
  string var_name[NumPar]; // = {"rN0","iN0","rN1","iN1","rL0","rL0"};
  for(int i=0;i<NPAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rN"<<i; iIO<<"iN"<<i; var_name[2*i] = rIO.str(); var_name[2*i+1] = iIO.str();}
  for(int i=0;i<APAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rA"<<i; iIO<<"iA"<<i; var_name[2*NPAR_FIT+2*i] = rIO.str(); var_name[2*NPAR_FIT+2*i+1] = iIO.str();}
  // Set the free variables to be minimized!
  for(int i=0;i<NumPar;i++) min->SetVariable(i,var_name[i],start_par[i], step);
 
  // do the minimization
  min->Minimize(); 
  const double *final_pars = min->X();

  TF1 *fInt = new TF1("int",intensity,LEFT_INT,RIGHT_INT,NumPar); fInt->SetParameters(final_pars);
  TF1 *fPhi = new TF1("phi",phase,    LEFT_INT,RIGHT_INT,NumPar); fPhi->SetParameters(final_pars);

  can.cd(1); fInt->Draw("same");
  can.cd(2); fPhi->Draw("same");
  
  can.SaveAs("c1.png");
  
  cout << "finished" << endl;
  for(int i=0;i<NPAR_FIT;i++) delete table[i];
  return 0;
}


////////////////////////////////////////////////////////////////////////
//                     /////////////////////////////////////////////////
////////////////////                   ///////////////            //////
/////////////////////////////////////                    ///////////////
////////////////////////////////////////////////////////////////////////

double globChi2(const double *par) {

  vector<cd> npar, apar;
  fillvectors(&par[0],apar,npar);

  // *********************** Calculate chi2 ********************** //
  double chi2 = 0;
  //Intensity
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT || mass>RIGHT_INT) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar,npar,S1,S2);
    double intens = norm(amp)*rho(mass,RHO_MASS,PI_MASS);;
    double err = 1;
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //Phase
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT || mass>RIGHT_INT) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar,npar,S1,S2);
    double phase = arg(amp);
    double err = 1;//gh2->GetBinError(i);
    double dPhi = normPhase(phase-gh2->GetBinContent(i));
    chi2 += pow(dPhi,2)/(err*err);
  }
  return chi2;
}

void fillvectors(const double *par, vector<cd> &apar, vector<cd> &npar) {

  //**************** Convert parameters to format **************** //
  //Constract parameters N
  for(int i=0;i<NPAR_FIT;i++) {
    cd ci(par[2*i],par[2*i+1]); 
    npar.push_back(ci);
  }
  //Constract Alpha
  for(int i=0;i<APAR_FIT;i++) {
    cd normi(par[2*NPAR_FIT+2*i],par[2*NPAR_FIT+2*i+1]);  
    apar.push_back(normi);
  }

}

double intensity(double *x, double *par) {

  vector<cd> npar, apar;
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,S1,S2);
  return norm(amp)*rho(x[0],RHO_MASS,PI_MASS);
}

double phase(double *x, double *par) {

  vector<cd> npar, apar;
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,S1,S2);
  return arg(amp);
}

cd BreitWignerA(double x,double m, double G0) {
  cd unit(0,1); 
  double p  = sqrt(lambda(x*x,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*x);
  double p0 = sqrt(lambda(m*m,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*m);
  double R = 5; double BW = (1+pow(R*p0,2))/(1+pow(R*p,2));
  double G = G0*pow(p/p0,2)*BW*rho(x,RHO_MASS,PI_MASS)/rho(m,RHO_MASS,PI_MASS);
  //cd back(2,2);
  cd A = sqrt(2*m*G)/(m*m-x*x-unit*m*G);// + back*p*exp(-p/2);
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
