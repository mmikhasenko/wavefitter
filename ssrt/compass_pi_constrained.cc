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
#include <TFile.h>
#include <TTree.h>
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

//region of fit
#define LEFT_INT1   1.0
#define RIGHT_INT1  2.4
#define LEFT_INT2   1.1
#define RIGHT_INT2  2.4
#define LEFT_PHI    1.5
#define RIGHT_PHI   2.4

using namespace std; 
typedef std::complex<double> cd;

cd BreitWignerA(double x,double m, double G0);
cd BreitWignerAndDeckA(double x,double m, double G0,double a, double b);
double intensity1(double *x, double *par);
double intensity2(double *x, double *par);
double dphase(double *x, double *par);

cd expand(cd s, cd s1, vector<cd> par);
cd vtable(double s,vector<pair<double,cd> > &arr);
cd Afull(cd s, vector<cd> apar, vector<cd> npar, cd s1, cd s2, double sth);
double globChi2(const double *par);
void prepareNames(int npar, int apar, string *var_name,const char *pref="");

void fillvectors(const double *par, vector<cd> &apar, vector<cd> &npar);
double normPhase(double a) {
  if(a > M_PI)  return normPhase(a-2*M_PI);
  if(a < -M_PI) return normPhase(a+2*M_PI);
  return a;
}
int find_zeros();


vector<vector<pair<double,cd> > > *garr;
TH1D *gh1, *gh2, *gh3;

int main(int ac, char** av) {  

  const int nAttempt = (ac > 1) ? atoi(av[1]) : 1;
  const char *fout_name = (ac>2) ? av[2] : "/tmp/test.root";

  //**************** Load integrals ****************//
  const int nTable = 7;
  string table_name[nTable] = {"f0.txt","f1.txt","f2.txt","f3.txt","f4.txt","f5.txt","f6.txt"};
  
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
  TGraph2D f2piS("data/f2piS.txt");  
  TGraph2D f2piD("data/f2piD.txt");  
  TGraph2D phi_f2piD_f2piS("data/phi_f2piD_f2piS.txt");
  const int nBins = 100;
  TH1D *f2piS_int = new TH1D("f2piS","f_{2}#pi S;M_{3#pi}",nBins,0.5,2.5);
  TH1D *f2piD_int = new TH1D("f2piD","f_{2}#pi D;M_{3#pi}",nBins,0.5,2.5);
  TH1D *f2pi_dphi = new TH1D("phi_f2piD_f2piS","#Delta(f_{2}#pi D and f_{2}#pi S);M_{3#pi}",nBins,0.5,2.5);
  if(f2piS.GetN() != f2piS_int->GetNbinsX()) {cout << "Very strange!" << endl; return 0;} 
  for(int i=1;i<=nBins;i++) {
    f2piS_int->SetBinContent(i,f2piS.GetY()[i]);
    f2piS_int->SetBinError  (i,f2piS.GetZ()[i]);
  }
  for(int i=1;i<=nBins;i++) {
    f2piD_int->SetBinContent(i,f2piD.GetY()[i]);
    f2piD_int->SetBinError  (i,f2piD.GetZ()[i]);
  }
  for(int i=1;i<=nBins;i++) {
    f2pi_dphi->SetBinContent(i,-normPhase(phi_f2piD_f2piS.GetY()[i]/180.*M_PI));
    f2pi_dphi->SetBinError  (i,phi_f2piD_f2piS.GetZ()[i]/180.*M_PI);
  }

  gh1 = f2piS_int;
  gh2 = f2piD_int;
  gh3 = f2pi_dphi;

  TCanvas *can = new TCanvas("c1","plots",0,0,1000,600);
  can->Divide(2,2);
 
  can->cd(1); gh1->Draw();
  can->cd(4); gh2->Draw();
  can->cd(2); gh3->Draw();


  //***************** Fit the data *****************//
  //Build minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  //min->SetMaxIterations(10000);  // for GSL 
  //min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  
  const int NumPar = 
    2*NPAR_FIT1 + 2*APAR_FIT1 +
    2*NPAR_FIT2 + 2*APAR_FIT2;
  
  // create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor f(&globChi2,NumPar); 
  min->SetFunction(f);
  

  double step = 0.01;
  //form names
  string var_name[NumPar];
  prepareNames(NPAR_FIT1,APAR_FIT1,&var_name[0],"ch1_");
  prepareNames(NPAR_FIT2,APAR_FIT2,&var_name[2*(NPAR_FIT1+APAR_FIT1)],"ch2_");
  // Set the free variables to be minimized!
  for(int i=0;i<NumPar;i++) min->SetVariable(i,var_name[i],0, step);
  

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
    // starting point
    double start_par[NumPar];
    gRandom->SetSeed(0);
    for(int i=0;i<NumPar;i++) start_par[i] = 2*gRandom->Rndm();
    /* remove umbig phase */  start_par[2*NPAR_FIT1+1]=0;    
    min->SetVariableValues(start_par);
    
    //Scheam of fixation and releasation parameters
    for(int i=0;i<NumPar;i++) min->FixVariable(i); 
    min->ReleaseVariable(2*NPAR_FIT1);       min->ReleaseVariable(2*NPAR_FIT1+2*APAR_FIT1+2*NPAR_FIT2);
    /*min->ReleaseVariable(2*NPAR_FIT1+1);*/ min->ReleaseVariable(2*NPAR_FIT1+2*APAR_FIT1+2*NPAR_FIT2+1);
    min->Minimize();
    for(int i=0;i<2*NPAR_FIT1&&i<4;i++) min->ReleaseVariable(i);
    for(int i=0;i<2*NPAR_FIT2&&i<4;i++) min->ReleaseVariable(2*NPAR_FIT1+2*APAR_FIT1+i);
    min->Minimize();
    for(int i=0;i<2*NPAR_FIT1;i++) min->ReleaseVariable(i);
    for(int i=0;i<2*NPAR_FIT2;i++) min->ReleaseVariable(2*NPAR_FIT1+2*APAR_FIT1+i);
    min->Minimize();
    for(int i=0;i<2*NPAR_FIT1;i++) min->FixVariable(i);
    for(int i=0;i<2*NPAR_FIT2;i++) min->FixVariable(2*NPAR_FIT1+2*APAR_FIT1+i);
    for(int i=0;i<2*APAR_FIT1;i++) if(i!=1) min->ReleaseVariable(2*NPAR_FIT1+i);
    for(int i=0;i<2*APAR_FIT2;i++) min->ReleaseVariable(2*NPAR_FIT1+2*APAR_FIT1+2*NPAR_FIT2+i);
    min->Minimize();
    //release all
    for(int i=0;i<NumPar;i++) min->ReleaseVariable(i); 
    min->FixVariable(2*NPAR_FIT1+1);
    bool fit_stat = min->Minimize();
    //done
    
    status = min->Status();
    if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);
    for(int i=0;i<NumPar;i++) cout << /*min->VariableName(i) << " = " << */final_pars[i] << endl;
    
    TF1 fInt1("int1",intensity1,LEFT_INT1,RIGHT_INT1,2*(NPAR_FIT1+APAR_FIT1)); fInt1.SetParameters(final_pars);
    TF1 fInt2("int2",intensity2,LEFT_INT2,RIGHT_INT2,2*(NPAR_FIT2+APAR_FIT2)); fInt2.SetParameters(&final_pars[2*(NPAR_FIT1+APAR_FIT1)]);
    TF1 fPhi("dphi",dphase,     LEFT_PHI, RIGHT_PHI, NumPar);                   fPhi.SetParameters(final_pars);
    
    TF1 fInt1_fr(fInt1); fInt1_fr.SetRange(RHO_MASS+PI_MASS,gh1->GetBinLowEdge(gh1->GetNbinsX())+gh1->GetBinWidth(1)); fInt1_fr.SetLineColor(kOrange);
    TF1 fInt2_fr(fInt2); fInt2_fr.SetRange(RHO_MASS+PI_MASS,gh1->GetBinLowEdge(gh1->GetNbinsX())+gh1->GetBinWidth(1)); fInt2_fr.SetLineColor(kOrange);
    TF1  fPhi_fr(fPhi);   fPhi_fr.SetRange(RHO_MASS+PI_MASS,gh1->GetBinLowEdge(gh1->GetNbinsX())+gh1->GetBinWidth(1));  fPhi_fr.SetLineColor(kOrange);
    can->cd(1); fInt1_fr.Draw("same"); fInt1.Draw("same");
    can->cd(2);  fPhi_fr.Draw("same");  fPhi.Draw("same");
    can->cd(4); fInt2_fr.Draw("same"); fInt2.Draw("same");
    
    if(nAttempt==1) can->SaveAs("c1.png");
    t.Fill();
    
  }

  t.Write();
  fout.Close();
  
  cout << "finished" << endl;
  for(int i=0;i<nTable;i++) delete table[i];
  return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//                     //////////////////////////////////////////            /////////////////////
////////////////////                   ///////////////            ////////////////////////////////
/////////////////////////////////////                    /////////////                      //////
//////////////////////////////////////////////////////////////////////////////////////////////////


int find_zeros() {

  //const char *fin_name = (ac >2) ? av[2] : "/tmp/test.root";
  
  //TFile *fin = TFile::Open(fin_name);
  
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//                     //////////////////////////////////////////            /////////////////////
////////////////////                   ///////////////            ////////////////////////////////
/////////////////////////////////////                    /////////////                      //////
//////////////////////////////////////////////////////////////////////////////////////////////////

double globChi2(const double *par) {

  vector<cd> npar1(NPAR_FIT1), apar1(APAR_FIT1); fillvectors(&par[0],apar1,npar1);
  vector<cd> npar2(NPAR_FIT2), apar2(APAR_FIT2); fillvectors(&par[2*(NPAR_FIT1+APAR_FIT1)],apar2,npar2);

  // *********************** Calculate chi2 ********************** //
  double chi2 = 0;
  const int startBin = 40;
  const int lastBin = 95;

  //Intensity1
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT1 || mass>RIGHT_INT1) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar1,npar1,S1,S2);
    double intens = norm(amp);
    double err = gh1->GetBinError(i);
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //Intensity2
  for(int i=startBin;i<=gh2->GetNbinsX();i++) {
    double mass = gh2->GetBinCenter(i);
    if(mass<LEFT_INT2 || mass>RIGHT_INT2) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar2,npar2,S1,S2);
    double intens = norm(amp);
    double err = gh2->GetBinError(i);
    chi2 += pow(intens - gh2->GetBinContent(i),2)/(err*err);
  }
  //dPhase
  for(int i=startBin;i<=gh3->GetNbinsX();i++) {
    double mass = gh3->GetBinCenter(i);
    if(mass<LEFT_PHI || mass>RIGHT_PHI) continue;
    double s = pow(mass,2);
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
  cd amp = Afull(x[0]*x[0],apar,npar,S1,S2);
  return norm(amp);
}
double intensity2(double *x, double *par) {
  vector<cd> npar(NPAR_FIT2), apar(APAR_FIT2);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,S1,S2);
  return norm(amp);
}

double dphase(double *x, double *par) {
  vector<cd> npar1(NPAR_FIT1), apar1(APAR_FIT1);
  vector<cd> npar2(NPAR_FIT2), apar2(APAR_FIT2);
  fillvectors(&par[0],apar1,npar1);
  fillvectors(&par[2*(NPAR_FIT1+APAR_FIT1)],apar2,npar2);
  cd amp1 = Afull(x[0]*x[0],apar1,npar1,S1,S2);
  cd amp2 = Afull(x[0]*x[0],apar2,npar2,S1,S2);
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

cd Afull(cd s, vector<cd> &apar, vector<cd> &npar, 
	 cd s1, cd s2, double sth) {

  //Constract D
  cd npar0 = 0; for(int i=0;i<npar.size();i++) npar0 -= expand(sth,s1,npar);
  cd dsum = npar0*vtable(real(s),(*garr)[0]);
  for(int i=0;i<npar.size();i++) dsum += npar[i]*vtable(real(s),(*garr)[i+1]);
  cd Denomin = 1. - s / (2*M_PI) * dsum;
  //Constract N, Alpha and A.
  cd Numir = npar0+expand(s,s1,npar);
  cd Aprod = expand(s,s2,apar);
  cd Afull = Aprod*Numir/Denomin;

  return Afull;
}
