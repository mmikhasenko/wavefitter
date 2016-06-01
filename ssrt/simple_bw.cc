#include <iostream>
#include <sstream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
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


#include <MIsobar.h>
#include <NoverD.h>

#include <constants.h>
#include <deflib.h>

//******Parameters of the model********//
#define DEFAULT_NPAR_FIT  4
#define APAR_FIT  1

#define DEFAULT_S1 0.5
#define S2 0.

#define MRANGE_NPAR 100.
#define MRANGE_APAR 1.

//******Parameters of the data*********//
#define R1_MASS  1.26
#define R1_WIDTH 0.2

#define LEFT_INT   RHO_MASS+PI_MASS
#define RIGHT_INT  2.

#define REL_SIGMA_MASS 0.02
#define REL_SIGMA_PHASE 0.02

using namespace std; 


//functions data
cd bw_swave_A(double x,double m, double G0);
double phase(double *x, double *par);
double intensity(double *x, double *par);

//functions fit
cd Afull(double s, const vector<cd> &apar, const vector<double> &npar, double s2);
cd rho(cd s);
double globChi2(const double *par);
void fillvectors(const double *par, vector<cd> &apar, vector<double> &npar);
double normPhase(double a);

//functions main
int fit_data(int nAttempt, const char *fout_name, double _s1, int _NPAR_FIT);
int find_poles(int nAttempt, const char *fin_name, const char *fout_name) {return 1;};
int plot_best_result(const char *fin_name) {return 1;};

int main(int ac, char **av) {

  if(ac<2) { cerr << "first arg is FIT or FIND or PLOT" << endl; return 1; }

  const int nAttempt = (ac > 2) ? atoi(av[2]) : 1;
  const char *file_name = (ac > 3) ? av[3] : "/tmp/test.root";
  const double s1 = (ac>4) ? atof(av[4]) : DEFAULT_S1;
  const int npar_fit = (ac>5) ? atoi(av[5]) : DEFAULT_NPAR_FIT;

  if(strcmp(av[1],"FIT")==0) {
    return fit_data(nAttempt,file_name,s1,npar_fit);
  } else if(strcmp(av[1],"FIND")==0) {
    return find_poles(nAttempt,file_name,"/tmp/find_poles.root");
  } else if(strcmp(av[1],"PLOT")==0) {
    return plot_best_result(file_name);
  } else { cerr << "first arg is FIT or FIND or PLOT" << endl; return 1; }
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

NoverD *gT;
TH1D *gh1, *gh2;
double gS1;
int gNPAR_FIT;

int fit_data(int nAttempt, const char *fout_name, double _s1, int _NPAR_FIT) {

  gS1 = _s1;
  gNPAR_FIT = _NPAR_FIT;

  //*****************Create the model***************//
  const double eTH = RHO_MASS+PI_MASS, sTH = pow(eTH,2);
  gT = new NoverD(gNPAR_FIT,gS1,rho,sTH,1000,300);

  //**************** Load the data *****************//
  double dataLeft = eTH;  double dataRight = 2.2; double dataNbin = 150;
  TH1D *hist_data_inten = new TH1D("intns","Intensity;M_{#rho#pi}^{2}",dataNbin,dataLeft,dataRight);
  TH1D *hist_data_phase = new TH1D("phase","Phase;M_{#rho#pi}^{2}",dataNbin,dataLeft,dataRight);
  double mass_inten = norm(bw_swave_A(R1_MASS,R1_MASS,R1_WIDTH))*gT->getrho(R1_MASS*R1_MASS);
  double sigma_int = REL_SIGMA_MASS  *mass_inten;
  double sigma_phi = REL_SIGMA_PHASE *2*M_PI;
  for(int i=1;i<=hist_data_inten->GetNbinsX();i++) {
    double e = hist_data_inten->GetBinCenter(i);  double s = pow(e,2);
    cd ampl =  bw_swave_A(e,R1_MASS,R1_WIDTH);
    double inten =  norm(ampl)*gT->getrho(s);
    hist_data_inten->SetBinContent(i, gRandom->Gaus(inten               ,sigma_int) );
    hist_data_phase->SetBinContent(i, gRandom->Gaus(normPhase(arg(ampl)),sigma_phi) );
    hist_data_inten->SetBinError  (i, sigma_int);
    hist_data_phase->SetBinError  (i, sigma_phi);
  }
  
  gh1 = hist_data_inten;
  gh2 = hist_data_phase; gh2->GetYaxis()->SetRangeUser(-M_PI,M_PI); 
  
  TCanvas *can = new TCanvas("c1","plots",0,0,1000,600);
  can->Divide(2,1);
  gh1->SetStats(kFALSE);
  gh2->SetStats(kFALSE);
  
  //***************** Fit the data *****************//
  //Build minimizer
  ROOT::Math::Minimizer* min = 
     ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");

  // set tolerance , etc...
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(0.001);
  min->SetStrategy(2);
  min->SetPrintLevel(1);
  
  const int NumPar = gNPAR_FIT + 2*APAR_FIT;
  
  // Create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor functor(&globChi2,NumPar);
  min->SetFunction(functor);
    
  const double step = 0.01;
  
  // Sorm names
  double limit[NumPar];
  string var_name[NumPar];
  for(int i=0;i<gNPAR_FIT;i++) { ostringstream rIO; rIO<<"rN"<<i; var_name[i] = rIO.str();}
  for(int i=0;i<APAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rA"<<i; iIO<<"iA"<<i; var_name[gNPAR_FIT+2*i] = rIO.str(); var_name[gNPAR_FIT+2*i+1] = iIO.str();}
  // Set the free variables to be minimized!

  for(int i=0;i<gNPAR_FIT;i++) limit[i] = MRANGE_NPAR/abs(gT->getmax_value_of_integral(i));
  for(int i=gNPAR_FIT;i<NumPar;i++) limit[i] = MRANGE_APAR*mass_inten;
  for(int i=0;i<NumPar;i++) min->SetLimitedVariable(i,var_name[i],0.0,step,-limit[i],limit[i]);

  //Do minimisation several times
  TFile fout(fout_name,"RECREATE");
  TTree t("mins","several_mimimas");
  double final_pars[NumPar];
  for(int i=0;i<NumPar;i++) t.Branch(var_name[i].c_str(),&final_pars[i]);
  t.Branch("canva","TCanvas",&can);
  double chi2; t.Branch("chi2",   &chi2);
  int status; t.Branch("status",&status);
  for(int e=0;e<nAttempt;e++) {
    cout << "---------- Attempt " << e << " -----------"<< endl;
    cout << "------------------------------------------"<< endl;
    // reload step size
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,step);
    // starting point
    double start_par[NumPar];
    gRandom->SetSeed(0);
    for(int i=0;i<gNPAR_FIT;i++) start_par[i] = limit[i]*(2*gRandom->Rndm()-1);
    for(int i=0;i<2*APAR_FIT;i++) start_par[gNPAR_FIT+i] = limit[i]*(2*gRandom->Rndm()-1);
    min->SetVariableValues(start_par);

    //do the minimization
    min->Minimize(); //min->Hesse();
    bool fit_stat = 1;
    status = min->Status(); //if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);

    TF1 *fInt = new TF1("int",intensity,LEFT_INT,RIGHT_INT,NumPar); fInt->SetParameters(final_pars);
    TF1 *fPhi = new TF1("phi",phase,    LEFT_INT,RIGHT_INT,NumPar); fPhi->SetParameters(final_pars);

    can->cd(1); gh1->Draw("e"); fInt->Draw("same");
    can->cd(2); gh2->Draw("e"); fPhi->Draw("same");
    for(int i=0;i<NumPar;i++) cout << var_name[i] << " = " << final_pars[i] << endl;

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

  vector<double> npar(gNPAR_FIT);
  vector<cd> apar(APAR_FIT);
  fillvectors(&par[0],apar,npar);

  // *********************** Calculate chi2 ********************** //
  double chi2 = 0;
  //Intensity
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT || mass>RIGHT_INT) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar,npar,S2);
    double intens = norm(amp)*gT->getrho(s);
    double err = gh1->GetBinError(i);
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //Phase
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT || mass>RIGHT_INT) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar,npar,S2);
    double phase = arg(amp);
    double err = gh2->GetBinError(i);
    double dPhi = normPhase(phase-gh2->GetBinContent(i));
    chi2 += pow(dPhi,2)/(err*err);
  }
  //  cout << "par[2] = "<< par[2] <<", chi2 = " << chi2<<endl;
  return chi2;
}

cd Afull(double s, const vector<cd> &apar, const vector<double> &npar, 
	 double s2) {
  for(int i=0;i<npar.size();i++) gT->npar[i] = npar[i];
  cd T = gT->A(s);
  cd Aprod = NoverD::cexpand(s,s2,apar);
  return Aprod*T;
}

cd rho(cd s) {
  if(real(s)>100 && fabs(imag(s))<1e-5) return 1./(8*M_PI)*(1.0-(pow(RHO_MASS,2)+pow(PI_MASS,2))/s);
  cd phspc2 = (s-pow(RHO_MASS+PI_MASS,2))*(s-pow(RHO_MASS-PI_MASS,2))/(s*s);
  return 1./(8*M_PI)*sqrtPi(phspc2);
}


////////////////////////////////////////////////////////////////////////
//                     /////////////////////////////////////////////////
////////////////////                   ///////////////            //////
/////////////////////////////////////                    ///////////////
////////////////////////////////////////////////////////////////////////

double intensity(double *x, double *par) {
  vector<double> npar(gNPAR_FIT);
  vector<cd> apar(APAR_FIT);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,S2);
  return norm(amp)*gT->getrho(x[0]*x[0]);
}

double phase(double *x, double *par) {
  vector<double> npar(gNPAR_FIT);
  vector<cd> apar(APAR_FIT);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(x[0]*x[0],apar,npar,S2);
  return arg(amp);
}

void fillvectors(const double *par, vector<cd> &apar, vector<double> &npar) {
  //Constract Numerator
  int np = npar.size();
  for(int i=0;i<np;i++) npar[i] = par[i];
  //Constract Alpha
  for(int i=0;i<apar.size();i++) apar[i] = cd(par[np+2*i],par[np+2*i+1]);
}

double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);
  if(a<-M_PI) return normPhase(a+2*M_PI);
  return a;
}

////////////////////////////////////////////////////////////////////////
//                     /////////////////////////////////////////////////
////////////////////                   ///////////////            //////
/////////////////////////////////////                    ///////////////
////////////////////////////////////////////////////////////////////////

cd bw_swave_A(double x,double m, double G0) { 
  double rho  = gT->getrho(x*x);
  double rho0 = gT->getrho(m*m);
  double G = G0*rho/rho0;
  cd A = sqrt(2*m*G0/rho0)/(m*m-x*x-cd(0,1)*m*G);
  return A;
} 
