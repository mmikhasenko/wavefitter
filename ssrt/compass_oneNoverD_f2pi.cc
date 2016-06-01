#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>

#include <TGraph2D.h>
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
#define S1 0.157
#define S2 0.0

#define MRANGE_NPAR 10.
#define MRANGE_APAR 10.

#define CUT_CHI2 1e5
#define MAX_INIT_ITERATION 1e3

//******Parameters of the data*********// 3*PI_MASS
#define LEFT_INT1   0.8
#define LEFT_INT2   0.8
#define LEFT_PHI    1.6
#define RIGHT_INT1  2.5
#define RIGHT_INT2  2.3
#define RIGHT_PHI   2.5

#define REL_SIGMA_MASS 0.02
#define REL_SIGMA_PHASE 0.02

using namespace std; 

//functions data
double dphase(double *x, double *par);
double intensity1(double *x, double *par);
double intensity2(double *x, double *par);

//functions fit
cd Afull(NoverD *lT,double s, const vector<cd> &apar, const vector<double> &npar, double s2);
cd rho(cd s);
double globChi2(const double *par);
void fillvectors(const double *par, vector<cd> &apar, vector<double> &npar);
void fillvectors(const double *par, vector<cd> &apar);
double normPhase(double a);
void prepareNames(int npar, int apar, string *var_name,const char *pref);
void printToFile(cd *arr, TH2D &t2d, ofstream &f);

//functions main
int fit_data(int nAttempt, const char *fout_name, int pars_all);
int plot_best_result(const char *fin_name) {;};
int plot_best_sheets(const char *fin_name);

int main(int ac, char **av) {

  if(ac<2) { cerr << "first arg is FIT or PLOT or SHEET" << endl; return 1; }

  const int nAttempt    = (ac > 2) ? atoi(av[2]) : 1;
  const char *file_name = (ac > 3) ? av[3] : "/tmp/test.root";
  const int pars_all    = (ac > 4) ? atoi(av[4]) : 4343;  

  if(strcmp(av[1],"FIT")==0) {
    return fit_data(nAttempt,file_name,pars_all);
  } else if(strcmp(av[1],"PLOT")==0) {
    return plot_best_result(file_name);
  } else if(strcmp(av[1],"SHEET")==0) {
    return plot_best_sheets(file_name);
  } else { cerr << "first arg is FIT or SHEET or PLOT" << endl; return 1; }
  
  return 0;
}

NoverD *gT;
TH1D *gh1, *gh2, *gh3;
int gNPAR_FIT, gAPAR_FIT1, gAPAR_FIT2;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int plot_best_sheets(const char *fin_name) {
  
  TFile *fin = TFile::Open(fin_name);           if(!fin) {cerr<<"No file!"<<endl;return 1;}
  TTree *t = (TTree*)gDirectory->Get("mins");   if(!t)   {cerr<<"No tree!"<<endl;return 1;}

  //Get parameters first
  TObjArray *br = t->GetListOfBranches(); br->ls();
  std::vector<TString> vName1;
  std::vector<TString> vName2;
  //real N pars
  int ipar=0;
  while(1) {
    TString svName = TString::Format("ch1_rN%d",ipar);
    if( br->FindObject(svName.Data()) ) vName1.push_back(svName); else break;
    std::cout << ipar << ": found ch1 N" << std::endl;
    ipar++;
  }
  //set result
  gNPAR_FIT = vName1.size();
  
  //Find entries number with minimal chi2
  const int N = t->GetEntries();
  double chi2; t->SetBranchAddress("chi2",&chi2); t->GetEntry(0);
  int emin=0;
  double chi2min = chi2;
  for(int i=1;i<N;i++) {
    if(i%100==0) std::cout << i << ", ";
    t->GetEntry(i);
    if(chi2<chi2min) { chi2min=chi2; emin=i;}
  } 
  cout << "\n-------------> Best chi2 = " << chi2min << endl;

  double vPar1[gNPAR_FIT];
  for(int i=0;i<vName1.size();i++) t->SetBranchAddress(vName1[i].Data(),&vPar1[i]);

  cout << "-------------> Creating models: " << gNPAR_FIT << endl;
  const double eTH = 3*PI_MASS, sTH = pow(eTH,2);
  gT = new NoverD(gNPAR_FIT,S1,rho,sTH,5000,500);

  //fill the model
  t->GetEntry(emin);
  for(int i=0;i<vName1.size();i++) {gT->npar[i] = vPar1[i]; cout << vPar1[i] << endl;}

  cout << "------------------> Calculation of sheets" << endl;
  const int Nbx=100, Nby=100;
  const double lx=0.5, rx=6.25, ly=-1.0, ry=1.0;
  TH2D t2d("t2d","t2d",Nbx,lx,rx,Nby,ly,ry);
  cd fs1[Nbx][Nby],ss1[Nbx][Nby],
     fs2[Nbx][Nby],ss2[Nbx][Nby];
  for(int i=1;i<=Nbx;i++) {
    cout << "Calculating " << i<<"/"<<Nbx << "..." << endl;
    for(int j=1;j<=Nby;j++) {
      cd s(t2d.GetXaxis()->GetBinCenter(i),t2d.GetYaxis()->GetBinCenter(j));
      cd DI1   = gT->DI(s);          fs1[i-1][j-1] = DI1;
      cd DII1  = DI1 - gT->Disc(s);  ss1[i-1][j-1] = DII1;
    }
  }

  ofstream tfout;
  tfout.open ("/tmp/sheets.out");
  tfout << "-----------> FIRST  SHEET OF 1" << endl;
  printToFile(&fs1[0][0],t2d,tfout);
  tfout << "-----------> SECOND SHEET OF 1" << endl;
  printToFile(&ss1[0][0],t2d,tfout);
  tfout.close();  

}


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int fit_data(int nAttempt, const char *fout_name, int pars_all) {

  gNPAR_FIT  = (pars_all % 1000) / 100;
  gAPAR_FIT1 = (pars_all % 100 ) / 10;
  gAPAR_FIT2 = (pars_all % 10  ) / 1;
  cout<<"Fit with model "<<gNPAR_FIT<<","<<gAPAR_FIT1<<gAPAR_FIT2<<endl;
  
  //**********************Seed**********************//
  gRandom->SetSeed(0); cout << "The seed " << gRandom->GetSeed() << " is used." << endl;

  //*****************Create the model***************//
  const double eTH = 3*PI_MASS, sTH = pow(eTH,2);
  gT = new NoverD(gNPAR_FIT,S1,rho,sTH,10000,500);

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
    f2piS_int->SetBinContent(i,f2piS.GetY()[i-1]);
    f2piS_int->SetBinError  (i,f2piS.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    f2piD_int->SetBinContent(i,f2piD.GetY()[i-1]);
    f2piD_int->SetBinError  (i,f2piD.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    f2pi_dphi->SetBinContent(i,-normPhase(phi_f2piD_f2piS.GetY()[i-1]/180.*M_PI));
    f2pi_dphi->SetBinError  (i,phi_f2piD_f2piS.GetZ()[i-1]/180.*M_PI);
  }

  gh1 = f2piS_int;
  gh2 = f2piD_int;
  gh3 = f2pi_dphi; gh3->GetYaxis()->SetRangeUser(-M_PI,M_PI); 
  //Normalise intensities
  double max = 1.0;
  max = gh1->GetBinContent(gh1->GetMaximumBin()); 
  for(int i=1;i<=gh1->GetNbinsX();i++) { gh1->SetBinContent(i,gh1->GetBinContent(i)/max); gh1->SetBinError(i,gh1->GetBinError(i)/max); }
  max = gh2->GetBinContent(gh2->GetMaximumBin()); 
  for(int i=1;i<=gh2->GetNbinsX();i++) { gh2->SetBinContent(i,gh2->GetBinContent(i)/max); gh2->SetBinError(i,gh2->GetBinError(i)/max); }

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","plots",0,0,1000,600);
  can->Divide(2,2);

  can->cd(1); gh1->SetStats(kFALSE);
  can->cd(4); gh2->SetStats(kFALSE);
  can->cd(2); gh3->SetStats(kFALSE);
  
  //***************** Fit the data *****************//
  //Build minimizer
  ROOT::Math::Minimizer* min = 
     ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  // set tolerance , etc...
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(0.001);
  min->SetStrategy(1);
  min->SetPrintLevel(1);
  min->Options().Print();
  
  const int NumPar = gNPAR_FIT+2*gAPAR_FIT1 +2*gAPAR_FIT2;
  
  // Create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor functor(&globChi2,NumPar);
  min->SetFunction(functor);
    
  const double step = 0.01;
  
  // Set the names
  string var_name[NumPar];
  prepareNames(gNPAR_FIT,gAPAR_FIT1,&var_name[0],"ch1_");
  prepareNames(0        ,gAPAR_FIT2,&var_name[gNPAR_FIT+2*gAPAR_FIT1],"ch2_");
  // Set the limites for the variables
  double limit[NumPar];
  for(int i=0;i<gNPAR_FIT;i++)       limit[i] = MRANGE_NPAR/abs(gT->getmax_value_of_integral(i));
  for(int i=0;i<2*gAPAR_FIT1;i++)    limit[i+gNPAR_FIT] = MRANGE_APAR*1.0;
  for(int i=0;i<2*gAPAR_FIT2;i++)    limit[i+gNPAR_FIT+2*gAPAR_FIT1] = MRANGE_APAR*1.0;
  // Set the variables to be minimized
  for(int i=0;i<NumPar;i++) min->SetVariable(i,var_name[i],0.0,step);//,-limit[i],limit[i] Limited

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
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,limit[i]/1e3);
    // starting point
    double start_pars[NumPar];
    int init_iterations = 0;
    chi2 = CUT_CHI2+1; double minimal_chi2 = 1e50; int best_seed = gRandom->GetSeed();
    while(init_iterations < MAX_INIT_ITERATION && chi2 > CUT_CHI2) {
      if(init_iterations>=MAX_INIT_ITERATION-1) {cerr<<"Warning: failed to find good starting point. Best seed "<<best_seed<<" will be used. chi2 = " << minimal_chi2 << endl; gRandom->SetSeed(best_seed);}
      for(int i=0;i<NumPar;i++) start_pars[i] = limit[i]*(2*gRandom->Rndm()-1);
      start_pars[gNPAR_FIT+1] = 0; min->FixVariable(gNPAR_FIT+1);
      min->SetVariableValues(start_pars);
      chi2 = globChi2(start_pars); if(chi2<minimal_chi2) {minimal_chi2=chi2; best_seed = gRandom->GetSeed();}
      init_iterations++;
    }

    //do the minimization
    //min->Minimize(); //min->Hesse();
    //for(int i=0;i<NumPar;i++)  min->SetVariableLimits(i,-1e10,1e10);
    min->Minimize(); //min->Hesse();
    bool fit_stat = 1;
    status = min->Status(); //if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);

    //Plot functions
    TF1 *fInt1 = new TF1("int1",intensity1,LEFT_INT1,RIGHT_INT1,NumPar); fInt1->SetParameters(final_pars);
    double final_cut_pars[gNPAR_FIT+2*gAPAR_FIT2];
    memcpy(&final_cut_pars[0]        ,&final_pars[0]                     ,  gNPAR_FIT *sizeof(double));
    memcpy(&final_cut_pars[gNPAR_FIT],&final_pars[gNPAR_FIT+2*gAPAR_FIT2],2*gAPAR_FIT2*sizeof(double));
    TF1 *fInt2 = new TF1("int2",intensity2,LEFT_INT2,RIGHT_INT2,NumPar); fInt2->SetParameters(final_cut_pars);
    TF1 *fPhi  = new TF1("phi" ,dphase,    LEFT_PHI ,RIGHT_PHI ,NumPar); fPhi ->SetParameters(final_pars);

    can->cd(1); gh1->Draw(); fInt1->Draw("same");
    can->cd(4); gh2->Draw(); fInt2->Draw("same");
    can->cd(2); gh3->Draw();  fPhi->Draw("same");
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
  vector<cd>     apar1(gAPAR_FIT2), apar2(gAPAR_FIT2);
  fillvectors(&par[0],apar1,npar);
  fillvectors(&par[gNPAR_FIT+2*gAPAR_FIT1],apar2);

  // *********************** Calculate chi2 ********************** //
  double chi2 = 0;
  //Intensity1
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT1 || mass>RIGHT_INT1) continue;
    double s = pow(mass,2);
    cd amp = Afull(gT,s,apar1,npar,S2);
    double intens = norm(amp)*gT->getrho(s);
    double err = gh1->GetBinError(i);
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //cout << "chi2 = " << chi2;
  //Intensity2
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double mass = gh2->GetBinCenter(i);
    if(mass<LEFT_INT2 || mass>RIGHT_INT2) continue;
    double s = pow(mass,2);
    cd amp = Afull(gT,s,apar2,npar,S2);
    double intens = norm(amp)*gT->getrho(s);
    double err = gh2->GetBinError(i);
    chi2 += pow(intens - gh2->GetBinContent(i),2)/(err*err);
  }
  //cout << ", second chi2 = " << chi2;
  //Phase
  for(int i=1;i<=gh3->GetNbinsX();i++) {
    double mass = gh3->GetBinCenter(i);
    if(mass<LEFT_PHI || mass>RIGHT_PHI) continue;
    double s = pow(mass,2);
    cd amp1 = Afull(gT,s,apar1,npar,S2);
    cd amp2 = Afull(gT,s,apar2,npar,S2);
    double phase = normPhase( arg(amp1*conj(amp2)) );
    double err = gh3->GetBinError(i);
    double dPhi = normPhase(phase-gh3->GetBinContent(i));
    chi2 += pow(dPhi,2)/(err*err);
  }
  //cout << ", final chi2 = " << chi2<<endl;
  return chi2;
}

cd Afull(NoverD *lT, double s, const vector<cd> &apar, const vector<double> &npar, double s2) {
  for(int i=0;i<npar.size();i++) lT->npar[i] = npar[i];
  cd T = lT->A(s);
  cd Aprod = NoverD::cexpand(s,s2,apar);
  return Aprod*T;
}

cd rho(cd s) {
  if(real(s)>1e3 && fabs(imag(s))<1e-5) return 1./(8*M_PI)*(1.0-(pow(F2_MASS,2)+pow(PI_MASS,2))/s);
  MIsobar RhoPi(F2_MASS,F2_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  cd value = RhoPi.rho3(s);
  return value;
}


////////////////////////////////////////////////////////////////////////
//                     /////////////////////////////////////////////////
////////////////////                   ///////////////            //////
/////////////////////////////////////                    ///////////////
////////////////////////////////////////////////////////////////////////

double intensity1(double *x, double *par) {
  vector<double> npar(gNPAR_FIT);
  vector<cd> apar(gAPAR_FIT1);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(gT,x[0]*x[0],apar,npar,S2);
  return norm(amp)*gT->getrho(x[0]*x[0]);
}

double intensity2(double *x, double *par) {
  vector<double> npar(gNPAR_FIT);
  vector<cd> apar(gAPAR_FIT2);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(gT,x[0]*x[0],apar,npar,S2);
  return norm(amp)*gT->getrho(x[0]*x[0]);
}

double dphase(double *x, double *par) {
  vector<double> npar(gNPAR_FIT);
  vector<cd>     apar1(gAPAR_FIT2), apar2(gAPAR_FIT2);
  fillvectors(&par[0],apar1,npar);
  fillvectors(&par[gNPAR_FIT+2*gAPAR_FIT1],apar2);
  cd amp1 = Afull(gT,x[0]*x[0],apar1,npar,S2);
  cd amp2 = Afull(gT,x[0]*x[0],apar2,npar,S2);
  return normPhase( arg(amp1*conj(amp2)) );
}

void fillvectors(const double *par, vector<cd> &apar, vector<double> &npar) {
  //Constract Numerator
  int np = npar.size();
  for(int i=0;i<np;i++) npar[i] = par[i];
  //Constract Alpha
  for(int i=0;i<apar.size();i++) apar[i] = cd(par[np+2*i],par[np+2*i+1]);
}
void fillvectors(const double *par, vector<cd> &apar) {for(int i=0;i<apar.size();i++) apar[i] = cd(par[2*i],par[2*i+1]);}

double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);
  if(a<-M_PI) return normPhase(a+2*M_PI);
  return a;
}

void prepareNames(int npar, int apar, string *var_name,const char *pref) {
  for(int i=0;i<npar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rN"<<i; var_name[i] = rIO.str();}
  for(int i=0;i<apar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rA"<<i; iIO<<pref<<"iA"<<i; var_name[npar+2*i] = rIO.str(); var_name[npar+2*i+1] = iIO.str();}
}

void printToFile(cd *arr, TH2D &t2d, ofstream &f) { 
  int Nbx = t2d.GetNbinsX();
  int Nby = t2d.GetNbinsY();
  for(int i=1;i<=Nbx;i++)
    for(int j=1;j<=Nby;j++) {
      int bin = t2d.GetBin(i,j);
      cd s(t2d.GetXaxis()->GetBinCenter(i),t2d.GetYaxis()->GetBinCenter(j));
      f << real(s) << " "
	<< imag(s) << " "
	<< real(arr[(i-1)*Nby+(j-1)]) << " "
	<< imag(arr[(i-1)*Nby+(j-1)]) << endl;
    }
}
 
