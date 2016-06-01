//Misha Mikasenko
#include <iostream>
#include <sstream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>

#include <TFile.h>
#include <TTree.h>
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
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include <constants.h>

#define NPAR_FIT  2
#define APAR_FIT  1

//#define S1 -10.0
#define S2 0.5

#define LEFT_INT   RHO_MASS+PI_MASS
#define RIGHT_INT  1.4

#define VFUNC_NPOINT 1000
#define VFUNC_RB pow(2.5,2)
#define VFUNC_LB pow(RHO_MASS+PI_MASS,2)


using namespace std; 
typedef std::complex<double> cd;

cd BreitWignerA(double x,double m, double G0);
double phase(double *x, double *par);
double intensity(double *x, double *par);

double intensity(double *x, double *par);
cd expand(cd s, cd s1, vector<cd> par);
cd vtable(double s,vector<pair<double,cd> > &arr);
cd Afull(cd s, vector<cd> apar, vector<cd> npar, double s1, double s2);
double globChi2(const double *par);

void build_integrals_table(vector<vector<pair<double,cd> > > &data);
void fillvectors(const double *par, vector<cd> &apar, vector<cd> &npar);
double normPhase(double a) {
  if(a>M_PI)  return a-2*M_PI*int( a/M_PI);
  if(a<-M_PI) return a+2*M_PI*int(-a/M_PI);
  return a;
}
double rho(double e, double m1, double m2) {
  return sqrt(lambda(e*e,m1*m1,m2*m2))/(e*e);
}
////////////////////////////////////////////
cd vint_sub(cd sp, cd s, double k, double s0, double sl, double sr);
double re_vint_sub(double *x, double *par);
double im_vint_sub(double *x, double *par);
cd integral(cd s, double k, double s0, double sl, double sr);
cd vf(cd s, double k, double s1, double sl, double sr);
////////////////////////////////////////////

vector<vector<pair<double,cd> > > *garr;
TH1D *gh1, *gh2;
double gs1;


int main(int ac, char **av) {

  const int nAttempt = (ac > 1) ? atoi(av[1]) : 1;
  const char *fout_name = (ac>2) ? av[2] : "/tmp/test.root";

  //**************** Calculate integrals ****************//
  vector<vector<pair<double,cd> > > data;

  gs1 = -9.;
  build_integrals_table(data);

  garr=&data;
  cout << "test access " << (*garr)[0][5].second << endl; 
  
  //**************** Load the data *****************//
  double dataLeft = RHO_MASS+PI_MASS;  double dataRight = 1.5; double dataNbin = 100;
  TH1D *hist_data = new TH1D("data","Intensity;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  TH1D *hist_data_phase = new TH1D("phase","Phase;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  for(int i=1;i<=hist_data->GetNbinsX();i++) {
    hist_data->SetBinContent(i,
			     norm(BreitWignerA(hist_data->GetBinCenter(i),
					       A1_MASS,A1_WIDTH)));
  }
  for(int i=1;i<=hist_data->GetNbinsX();i++) {
    hist_data_phase->SetBinContent(i,
				   arg(BreitWignerA(hist_data->GetBinCenter(i),
						    A1_MASS,A1_WIDTH)));

  }

  
  gh1 = hist_data;
  gh2 = hist_data_phase;

  TCanvas *can = new TCanvas("c1","plots",0,0,1000,600);
  can->Divide(2,1);
 
  //***************** Fit the data *****************//
  //Build minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  //min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  
  const int NumPar = 2*NPAR_FIT + 2*APAR_FIT + 1;
  
  // create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor f(&globChi2,NumPar); 
  min->SetFunction(f);
    
  double step = 0.01;
  //form names
  string var_name[NumPar]; // = {"rN0","iN0","rN1","iN1","rL0","rL0"};
  for(int i=0;i<NPAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rN"<<i; iIO<<"iN"<<i; var_name[2*i] = rIO.str(); var_name[2*i+1] = iIO.str();}
  for(int i=0;i<APAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rA"<<i; iIO<<"iA"<<i; var_name[2*NPAR_FIT+2*i] = rIO.str(); var_name[2*NPAR_FIT+2*i+1] = iIO.str();}
  var_name[NumPar-1] = "sn1"; 
  // Set the free variables to be minimized!
  for(int i=0;i<NumPar;i++) min->SetVariable(i,var_name[i],0, step);
  //right singularities
  min->SetVariableUpperLimit(NumPar-1,pow(RHO_MASS+PI_MASS,2));

  //Do minimisation several times
  TFile fout(fout_name,"RECREATE");
  TTree t("mins","several_mimimas");
  double final_pars[NumPar];
  for(int i=0;i<NumPar;i++) t.Branch(var_name[i].c_str(),&final_pars[i]);
  t.Branch("canva","TCanvas",&can);
  double chi2; t.Branch("chi2",   &chi2);
  int status; t.Branch("status",&status);
  for(int e=0;e<nAttempt;e++) {

    //reload step size
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,step);
    min->SetVariableStepSize(NumPar-1,0.1);
    // starting point
    double start_par[NumPar]; gRandom->SetSeed(0);
    for(int i=0;i<NumPar;i++) start_par[i] = 2*gRandom->Rndm();
    start_par[NumPar-1] = pow(RHO_MASS+PI_MASS,2)-1./gRandom->Rndm();
    min->SetVariableValues(start_par);

    min->FixVariable(NumPar-1); min->Minimize(); 
    min->ReleaseVariable(NumPar-1); for(int i=0;i<NumPar-1;i++) min->FixVariable(i);     min->Minimize(); 
    min->FixVariable(NumPar-1);     for(int i=0;i<NumPar-1;i++) min->ReleaseVariable(i); min->Minimize(); 
    min->ReleaseVariable(NumPar-1);

    // do the minimization
    bool fit_stat = min->Minimize(); 
    status = min->Status(); if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);

    TF1 *fInt = new TF1("int",intensity,LEFT_INT,RIGHT_INT,NumPar); fInt->SetParameters(final_pars);
    TF1 *fPhi = new TF1("phi",phase,    LEFT_INT,RIGHT_INT,NumPar); fPhi->SetParameters(final_pars);

    can->cd(1); gh1->Draw(); fInt->Draw("same");
    can->cd(2); gh2->Draw(); fPhi->Draw("same");

    if(nAttempt==1) can->SaveAs("c1.png");
    t.Fill();
  }

  t.Write();
  fout.Close();
  
  cout << "finished" << endl;
  return 0;
}
  

////////////////////////////////////////////////////////////////////////
//                     /////////////////////////////////////////////////
////////////////////                   ///////////////            //////
/////////////////////////////////////                    ///////////////
////////////////////////////////////////////////////////////////////////

double globChi2(const double *par) {

  vector<cd> npar(NPAR_FIT), apar(APAR_FIT);
  fillvectors(&par[0],apar,npar);

  // *********************** Calculate chi2 ********************** //
  double chi2 = 0;
  //Intensity
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT || mass>RIGHT_INT) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar,npar,par[2*NPAR_FIT+2*APAR_FIT],S2);
    double intens = norm(amp)*rho(mass,RHO_MASS,PI_MASS);;
    double err = 1;
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //Phase
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT || mass>RIGHT_INT) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar,npar,par[2*NPAR_FIT+2*APAR_FIT],S2);
    double phase = arg(amp);
    double err = 1;//gh2->GetBinError(i);
    double dPhi = normPhase(phase-gh2->GetBinContent(i));
    chi2 += pow(dPhi,2)/(err*err);
  }
  //cout << "chi2 = " << chi2<<endl;
  return chi2;
}


void build_integrals_table(vector<vector<pair<double,cd> > > &data) {
  data.clear();
  for(int i=0;i<NPAR_FIT;i++) {
    vector<pair<double,cd> > vj;
    double step = 1.0*(VFUNC_RB - VFUNC_LB)/(VFUNC_NPOINT-1);
    for(int j=0;j<VFUNC_NPOINT;j++) {
      double s = VFUNC_LB + step*j;
      cd value = vf(s,i,gs1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));
      vj.push_back(make_pair<double,cd>(s,value));
    }
    data.push_back(vj);
  }
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

double intensity(double *x, double *par) {

  vector<cd> npar(NPAR_FIT), apar(APAR_FIT);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,par[2*NPAR_FIT+2*APAR_FIT],S2);
  return norm(amp)*rho(x[0],RHO_MASS,PI_MASS);
}

double phase(double *x, double *par) {

  vector<cd> npar(NPAR_FIT), apar(APAR_FIT);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,par[2*NPAR_FIT+2*APAR_FIT],S2);
  return arg(amp);
}

cd BreitWignerA(double x,double m, double G0) {
  cd unit(0,1); 
  double p  = sqrt(lambda(x*x,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*x);
  double p0 = sqrt(lambda(m*m,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*m);
  double R = 5; double BW = (1+pow(R*p0,2))/(1+pow(R*p,2));
  double G = G0*pow(p/p0,2)*BW*rho(x,RHO_MASS,PI_MASS)/rho(m,RHO_MASS,PI_MASS);

  cd A = sqrt(2*m*G)/(m*m-x*x-unit*m*G);
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
    (arr[Nstep+1].first -arr[Nstep].first ) * (s-arr[Nstep].first);
}

cd Afull(cd s, vector<cd> apar, vector<cd> npar, 
	 double s1, double s2) {

  //Constract D
  cd dsum =0;
  if(s1!=gs1) {/*cout<<"------------->--rebuild all integrals, s1="<<s1<<"--<---------------"<<endl;*/ gs1=s1; build_integrals_table(*garr); }
  for(int i=0;i<npar.size();i++) dsum += npar[i]*vtable(real(s),(*garr)[i]);//vf(s,i,s1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));//
  cd Denomin = 1. - s / (2*M_PI) * dsum;
  //Constract N, Alpha and A.
  cd Numir = expand(s,s1,npar);
  cd Aprod = expand(s,s2,apar);
  cd Afull = Aprod*Numir/Denomin;

  return Afull;
}

//////////////////////////////////////////////////////////////////////////////                  
///  //   /  ///      //      ///     ///      ///      ///  ///////     /////
///  //      /////  ////  ///////  //////  //  ///  //  ///  ///////   ///////
///  //  /   /////  ////    /////  /  ///    /////      ///  /////////   /////
///  //  //  /////  ////      ///     ///  /   ///  //  ///      ///     /////
//////////////////////////////////////////////////////////////////////////////


cd vint_sub(cd sp, cd s, double k, double s0, double sl, double sr) {
  cd rho_sp2 = (sp-sl)*(sp-sr)/(sp*sp);
  cd rho_s2  = (s -sl)*(s -sr)/(s*s);
  return (pow(w(sp,s0),k)*sqrt(rho_sp2)-pow(w(s,s0),k)*sqrt(rho_s2))/(sp*(sp-s));
}

double re_vint_sub(double *x, double *par) {
  cd spr(1/x[0],0.), sr(par[0],par[1]);
  return real(vint_sub(spr,sr,par[2],par[3],par[4],par[5]))/(x[0]*x[0]);
}
double im_vint_sub(double *x, double *par) {
  cd spr(1/x[0],0.), sr(par[0],par[1]);
  return imag(vint_sub(spr,sr,par[2],par[3],par[4],par[5]))/(x[0]*x[0]);
}

cd integral(cd s, double k, double s0, double sl, double sr) {

  TF1 fre("fre", re_vint_sub, 0, 1./sr,6);
  TF1 fim("fim", im_vint_sub, 0, 1./sr,6);
  double pars[] = {real(s),imag(s),k,s0,sl,sr};
  fre.SetParameters(pars);
  fim.SetParameters(pars);

  ROOT::Math::WrappedTF1 wre(fre);  ROOT::Math::GaussIntegrator igre;
  ROOT::Math::WrappedTF1 wim(fim);  ROOT::Math::GaussIntegrator igim;
  
  igre.SetFunction(wre); igre.SetRelTolerance(0.01);
  igim.SetFunction(wim); igim.SetRelTolerance(0.01);

  cd res(igre.Integral(0, 1./sr),igim.Integral(0, 1./sr));
  return res;
}

cd vf(cd s, double k, double s0, double sl, double sr) {
  //cout << "integral" << endl;
  cd rho_s2  = (s -sl)*(s -sr)/(s*s);
  cd mult = pow(w(s,s0),k)*sqrt(rho_s2)/s;
  cd unit(0,1);
  cd alog = 1.-(s+unit*EPSILON)/sr;
  cd first_int = integral(s,k,s0,sl,sr);
  return first_int-mult*log(alog);
}
