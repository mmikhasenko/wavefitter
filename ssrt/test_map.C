//Misha Mikasenko
#include <iostream>
#include <complex>

#include <TH1D.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <Math/MinimizerOptions.h>
#include <TVirtualFitter.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include <constants.h>

#define NPAR_FIT  2
#define APAR_FIT  1

#define S1 0.5
#define S2 0.5

using namespace std; 
typedef std::complex<double> cd;

double BreitWignerI(double x,double m, double G0);
double intensity(double *x, double *par);

cd expand(cd s, cd s1, vector<cd> par);
cd vtable(double s,vector<pair<double,cd> > &arr);
cd Afull(cd s, vector<cd> apar, vector<cd> npar, cd s1, cd s2);

vector<vector<pair<double,cd> > > *garr;

int main() {

  //**************** Load integrals ****************//
  string table_name[] = {"f0.txt","f1.txt","f2.txt"};
  
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
  TH1D *hist_data = new TH1D("data","M_{3#pi}^{2}",100,1,2);
  for(int i=1;i<=hist_data->GetNbinsX();i++) {
    hist_data->SetBinContent(i,
			     BreitWignerI(hist_data->GetBinCenter(i),
					  A1_MASS,A1_WIDTH));

  }
  
  TCanvas can("c1");
  hist_data->Draw();
  
  //***************** Fit the data *****************//
  const int NumPar = 2*NPAR_FIT + 2*APAR_FIT - 1;
  TF1 *f = new TF1("fit",intensity,1.2,1.8,NumPar);
  double start_par[NumPar]; 
  gRandom->SetSeed(0);
  for(int i=0;i<NumPar;i++) start_par[i] = gRandom->Rndm();
  f->SetParameters(start_par);
  f->SetNpx(500);
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
  TVirtualFitter::SetMaxIterations(5000);
  hist_data->Fit(f,"R0");
  hist_data->Fit(f,"R0");
  //f->SetParameters(f->GetParameters());
  //f->Draw("same");

  TF1 *f2 = new TF1("fit_draw",intensity,1,2,NumPar);
  f2->SetParameters(f->GetParameters());
  f2->SetLineColor(46);
  f2->Draw("same");
  
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

double intensity(double *x, double *par) {

  //Constract N
  vector<cd> npar;
  for(int i=0;i<NPAR_FIT;i++) {
    cd ci(par[2*i],par[2*i+1]); 
    npar.push_back(ci);
  }
  //Constract Alpha
  vector<cd> apar;
  cd norm1(par[2*NPAR_FIT],0); apar.push_back(norm1);
  for(int i=1;i<APAR_FIT;i++) {
    cd normi(par[2*NPAR_FIT+2*i+1],par[2*NPAR_FIT+2*i+2]);  
    apar.push_back(normi);
  }
  cd amp = Afull(x[0],apar,npar,S1,S2);
  return norm(amp);
}

double BreitWignerI(double x,double m, double G0) {
  cd unit(0,1); 
  double p  = sqrt(lambda(x,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*sqrt(x));
  double p0 = sqrt(lambda(m*m,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*m);
  double G = G0*p/p0;
  cd back(2,2);
  cd A = sqrt(m*G)/(m*m-x-unit*m*G) + back*p*exp(-p/2);
  return norm(A);
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
