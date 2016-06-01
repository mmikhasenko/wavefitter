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
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include <constants.h>

#define NPAR_FIT1  3
#define APAR_FIT1  2
#define NPAR_FIT2  3
#define APAR_FIT2  2

#define S1 -9.
#define S2 -0.1

//region of fit
#define LEFT_INT1   1.42
#define RIGHT_INT1  2.4
#define LEFT_INT2   1.42
#define RIGHT_INT2  2.4
#define LEFT_PHI    1.5
#define RIGHT_PHI   2.4

#define VFUNC_NPOINT 1000
#define VFUNC_RB pow(2.5,2)
#define VFUNC_LB pow(F2_MASS+PI_MASS,2)

using namespace std; 
typedef std::complex<double> cd;

double intensity1(double *x, double *par);
double intensity2(double *x, double *par);
double dphase(double *x, double *par);

cd expand(cd s, cd s1, vector<cd> par);
cd vtable(double s,vector<pair<double,cd> > &arr);
cd Afull(cd s, vector<cd> apar, vector<cd> npar, double s1, double s2);
double globChi2(const double *par);
void prepareNames(int npar, int apar, string *var_name,const char *pref="");

void fillvectors(const double *par, vector<cd> &apar, vector<cd> &npar);
double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);
  if(a<-M_PI) return normPhase(a+2*M_PI);
  return a;
}
double rho(double e, double m1, double m2) {
  return sqrt(LAMBDA(e*e,m1*m1,m2*m2))/(e*e);
}
////////////////////////////////////////////
void build_integrals_table(vector<vector<pair<double,cd> > > &data, double s1);
cd vint_sub(cd sp, cd s, double k, double s0, double sl, double sr);
double re_vint_sub(double *x, double *par);
double im_vint_sub(double *x, double *par);
cd integral(cd s, double k, double s0, double sl, double sr);
cd vf(cd s, double k, double s1, double sl, double sr);
////////////////////////////////////////////


vector<vector<pair<double,cd> > > *garr;
TH1D *gh1, *gh2, *gh3;
double gs1;

int main(int ac, char** av) {  

  const int nAttempt = (ac > 1) ? atoi(av[1]) : 1;
  const char *fout_name = (ac>2) ? av[2] : "/tmp/test.root";

  //**************** Calculate integrals ****************//
  vector<vector<pair<double,cd> > > data;

  gs1 = S1;
  build_integrals_table(data,gs1);

  garr=&data;
  cout << "test access " << (*garr)[0][5].second << endl; 
  
  //**************** Load the data *****************//
  TGraph2D f2piS("data/f2piS.txt"); f2piS.SetName("gf2piS");  
  TGraph2D f2piD("data/f2piD.txt"); f2piD.SetName("gf2piD"); 
  TGraph2D phi_f2piD_f2piS("data/phi_f2piD_f2piS.txt"); phi_f2piD_f2piS.SetName("gphi_f2piD_f2piS"); 
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
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  //min->SetMaxIterations(10000);  // for GSL 
  //min->SetTolerance(0.001);
  min->SetPrintLevel(2);
  
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
  // Calculate limits
  double limit[NumPar];
  for(int i=0;i<2*NPAR_FIT1;i++) limit[i]=1000.;//abs(vf(2.5,i/2,0,pow(F2_MASS-PI_MASS,2),pow(F2_MASS+PI_MASS,2)));
  for(int i=0;i<2*APAR_FIT1;i++) limit[2*NPAR_FIT1+i]=5 * sqrt(gh1->GetBinContent(gh1->GetMaximumBin()));
  for(int i=0;i<2*NPAR_FIT2;i++) limit[2*(NPAR_FIT1+APAR_FIT1)+i]=1000.;//abs(vf(2.5,i/2,0,pow(F2_MASS-PI_MASS,2),pow(F2_MASS+PI_MASS,2))); 
  for(int i=0;i<2*APAR_FIT2;i++) limit[2*(NPAR_FIT1+APAR_FIT1+NPAR_FIT2)+i]=5 * sqrt(gh2->GetBinContent(gh2->GetMaximumBin()));
  cout << "Limits: ";
  for(int i=0;i<NumPar;i++) cout << limit[i] << " ";
  cout << endl;		      
  // Set the free variables to be minimized
  for(int i=0;i<NumPar;i++) min->SetLimitedVariable(i,var_name[i],0, step,-limit[i], limit[i]);  

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
    gRandom->SetSeed(0); cout << "------------- Seed is set to "<<gRandom->GetSeed()<<" ----------------"<<endl;
    for(int i=0;i<NumPar;i++) if(i%2==0) start_par[i] = limit[i]*(2*gRandom->Rndm()-1); else start_par[i] = 0;
    /* remove umbig phase */  start_par[2*NPAR_FIT1+1]=0;    
    min->SetVariableValues(start_par);
    //fix imag parts
    for(int i=0;i<NPAR_FIT1;i++) min->FixVariable(2*i+1);
    for(int i=0;i<NPAR_FIT2;i++) min->FixVariable(2*NPAR_FIT1+2*APAR_FIT1+2*i+1);
    min->FixVariable(2*NPAR_FIT1+1);


    //Scheam of fixation and releasation parameters
    for(int i=NPAR_FIT1+APAR_FIT1;i<NumPar/2;i++) min->FixVariable(2*i);
    min->Minimize();
    for(int i=NPAR_FIT1+APAR_FIT1;i<NumPar/2;i++) min->ReleaseVariable(2*i);
    for(int i=0;i<NPAR_FIT1+APAR_FIT1;i++) min->FixVariable(2*i);
    min->Minimize();
    for(int i=0;i<NPAR_FIT1+APAR_FIT1;i++) min->ReleaseVariable(2*i);
   
    bool fit_stat = min->Minimize();
    //done
    
    status = min->Status();
    //if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);
    for(int i=0;i<NumPar;i++) cout << min->VariableName(i) << " = " << final_pars[i] << endl;
    
    TF1 fInt1("int1",intensity1,LEFT_INT1,RIGHT_INT1,2*(NPAR_FIT1+APAR_FIT1)); fInt1.SetParameters(final_pars);
    TF1 fInt2("int2",intensity2,LEFT_INT2,RIGHT_INT2,2*(NPAR_FIT2+APAR_FIT2)); fInt2.SetParameters(&final_pars[2*(NPAR_FIT1+APAR_FIT1)]);
    TF1 fPhi("dphi",dphase,     LEFT_PHI, RIGHT_PHI, NumPar);                   fPhi.SetParameters(final_pars);
    
    double mass_r = gh1->GetBinLowEdge(gh1->GetNbinsX())+gh1->GetBinWidth(1);
    //TF1 fInt1_fr(fInt1); fInt1_fr.SetRange(F2_MASS+PI_MASS,mass_r); fInt1_fr.SetLineColor(kOrange); fInt1_fr.SetParameters(final_pars);
    //TF1 fInt2_fr(fInt2); fInt2_fr.SetRange(F2_MASS+PI_MASS,mass_r); fInt2_fr.SetLineColor(kOrange); fInt2_fr.SetParameters(&final_pars[2*(NPAR_FIT1+APAR_FIT1)]);
    //TF1  fPhi_fr(fPhi);   fPhi_fr.SetRange(F2_MASS+PI_MASS,mass_r);  fPhi_fr.SetLineColor(kOrange);  fPhi_fr.SetParameters(final_pars);
    can->cd(1);/* fInt1_fr.Draw("same");*/ fInt1.Draw("same");
    can->cd(2);/*  fPhi_fr.Draw("same");*/  fPhi.Draw("same");
    can->cd(4);/* fInt2_fr.Draw("same");*/ fInt2.Draw("same");
    
    if(nAttempt==1) can->SaveAs("c1.png");
    t.Fill();
    
  }

  t.Write();
  fout.Close();
  
  cout << "finished" << endl;
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
    double intens = norm(amp)*rho(mass,F2_MASS,PI_MASS);
    double err = gh1->GetBinError(i);
    //cout << chi2<<", "<<intens<<", "<<err << endl;
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //cout << "chi2 = " << chi2 << endl;
  //Intensity2
  for(int i=startBin;i<=gh2->GetNbinsX();i++) {
    double mass = gh2->GetBinCenter(i);
    if(mass<LEFT_INT2 || mass>RIGHT_INT2) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar2,npar2,S1,S2);
    double intens = norm(amp)*rho(mass,F2_MASS,PI_MASS);
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

void build_integrals_table(vector<vector<pair<double,cd> > > &data, double s1) {
  data.clear();
  for(int i=0;i<NPAR_FIT1;i++) {
    vector<pair<double,cd> > vj;
    double step = 1.0*(VFUNC_RB - VFUNC_LB)/(VFUNC_NPOINT-1);
    for(int j=0;j<VFUNC_NPOINT;j++) {
      double s = VFUNC_LB + step*j;
      cd value = vf(s,i,s1,pow(F2_MASS-PI_MASS,2),pow(F2_MASS+PI_MASS,2));
      vj.push_back(make_pair<double,cd>(double(s),cd(value)));
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

void prepareNames(int npar, int apar, string *var_name,const char *pref) {
  for(int i=0;i<npar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rN"<<i; iIO<<pref<<"iN"<<i; var_name[2*i] = rIO.str(); var_name[2*i+1] = iIO.str();}
  for(int i=0;i<apar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rA"<<i; iIO<<pref<<"iA"<<i; var_name[2*npar+2*i] = rIO.str(); var_name[2*npar+2*i+1] = iIO.str();}
}

double intensity1(double *x, double *par) {
  vector<cd> npar(NPAR_FIT1), apar(APAR_FIT1);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,S1,S2);
  return norm(amp)*rho(x[0],F2_MASS,PI_MASS);
}
double intensity2(double *x, double *par) {
  vector<cd> npar(NPAR_FIT2), apar(APAR_FIT2);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,S1,S2);
  return norm(amp)*rho(x[0],F2_MASS,PI_MASS);
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
	 double s1, double s2) {

  //Constract D
  cd dsum =0;
  if(s1!=gs1) { cout << "rebuils integrals, s1 = "<<s1 << endl; gs1=s1; build_integrals_table(*garr, s1); }
  for(int i=0;i<npar.size();i++) dsum += npar[i]*vtable(real(s),(*garr)[i]);
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
  cd rho_s2  = (s -sl)*(s -sr)/(s*s);
  cd mult = pow(w(s,s0),k)*sqrt(rho_s2)/s;
  cd unit(0,1);
  cd alog = 1.-(s+unit*EPSILON)/sr;
  cd first_int = integral(s,k,s0,sl,sr);
  return first_int-mult*log(alog);
}
