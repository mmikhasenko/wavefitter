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
#include <MStructureHolder.h>

#include <constants.h>
#include <deflib.h>

//******Parameters of the model********//
#define S1 0.157
#define S2 0.0

#define MRANGE_NPAR 1.
#define MRANGE_APAR 1.

#define CUT_CHI2 1e5
#define MAX_INIT_ITERATION 1

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
double normPhase(double a);
void prepareNames(int npar, int apar, string *var_name,const char *pref);
void printToFile(cd *arr, TH2D &t2d, ofstream &f);

//functions main
int fit_data(int nAttempt, const char *fout_name, int pars_all);
int plot_best_result(const char *fin_name);
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

NoverD *gT1, *gT2;
TH1D *gh1, *gh2, *gh3;
int gNPAR_FIT1, gAPAR_FIT1, gNPAR_FIT2, gAPAR_FIT2;
vector< pair<double,double> > glookup_rho_F2;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int plot_best_result(const char *fin_name) {
  
  TFile *fin = TFile::Open(fin_name);           if(!fin) {cerr<<"No file!"<<endl;return 1;}
  TTree *t = (TTree*)gDirectory->Get("mins");   if(!t)   {cerr<<"No tree!"<<endl;return 1;}

  //Get parameters first
  TObjArray *br = t->GetListOfBranches(); br->ls();
  std::vector<TString> vName1; std::vector<TString> vName2;
  std::vector<TString> rName1; std::vector<TString> rName2;
  std::vector<TString> iName1; std::vector<TString> iName2;
  //real N pars
  int ipar=0;
  ipar=0; while(1) { TString svName = TString::Format("ch1_rN%d",ipar); if( br->FindObject(svName.Data()) ) vName1.push_back(svName); else break; std::cout << ipar << ": found ch1 N" << std::endl; ipar++;}
  ipar=0; while(1) { TString svName = TString::Format("ch2_rN%d",ipar); if( br->FindObject(svName.Data()) ) vName2.push_back(svName); else break; std::cout << ipar << ": found ch2 N" << std::endl; ipar++;}
  //real A pars
  ipar=0; while(1) { TString svName = TString::Format("ch1_rA%d",ipar); if( br->FindObject(svName.Data()) ) rName1.push_back(svName); else break; std::cout << ipar << ": found ch1 rA" << std::endl; ipar++;}
  ipar=0; while(1) { TString svName = TString::Format("ch2_rA%d",ipar); if( br->FindObject(svName.Data()) ) rName2.push_back(svName); else break; std::cout << ipar << ": found ch2 rA" << std::endl; ipar++;}
  //real A pars
  ipar=0; while(1) { TString svName = TString::Format("ch1_iA%d",ipar); if( br->FindObject(svName.Data()) ) iName1.push_back(svName); else break; std::cout << ipar << ": found ch1 iA" << std::endl; ipar++;}
  ipar=0; while(1) { TString svName = TString::Format("ch2_iA%d",ipar); if( br->FindObject(svName.Data()) ) iName2.push_back(svName); else break; std::cout << ipar << ": found ch2 iA" << std::endl; ipar++;}

  //set result
  gNPAR_FIT1 = vName1.size(); gNPAR_FIT2 = vName2.size();
  gAPAR_FIT1 = rName1.size(); gAPAR_FIT2 = rName2.size();
  if(gAPAR_FIT1!=iName1.size()&&gAPAR_FIT2!=iName2.size()) {cerr<<"Error: Something wrong!"<<endl; return 1;}
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
  cout << "\n-------------> Best chi2 = " << chi2min << " ("<<emin<<")" << endl;

  double vPar1[gNPAR_FIT1]; double vPar2[gNPAR_FIT2];
  for(int i=0;i<vName1.size();i++) t->SetBranchAddress(vName1[i].Data(),&vPar1[i]);
  for(int i=0;i<vName2.size();i++) t->SetBranchAddress(vName2[i].Data(),&vPar2[i]);
  double rPar1[gAPAR_FIT1]; double rPar2[gAPAR_FIT2];
  for(int i=0;i<rName1.size();i++) t->SetBranchAddress(rName1[i].Data(),&rPar1[i]);
  for(int i=0;i<rName2.size();i++) t->SetBranchAddress(rName2[i].Data(),&rPar2[i]);
  double iPar1[gAPAR_FIT1]; double iPar2[gAPAR_FIT2];
  for(int i=0;i<iName1.size();i++) t->SetBranchAddress(iName1[i].Data(),&iPar1[i]);
  for(int i=0;i<iName2.size();i++) t->SetBranchAddress(iName2[i].Data(),&iPar2[i]);

  cout << "-------------> Creating models: " << gNPAR_FIT1 << "," << gNPAR_FIT2  << endl;
  const double eTH = 3*PI_MASS, sTH = pow(eTH,2);
  gT1 = new NoverD(gNPAR_FIT1,0,S1,rho,sTH,2000,500);
  gT2 = new NoverD(gNPAR_FIT2,0,S1,rho,sTH,2000,500);

  //fill the model
  t->GetEntry(emin);
  for(int i=0;i<vName1.size();i++) {gT1->npar[i] = vPar1[i]; cout << vPar1[i] << " ";} cout<<endl;
  for(int i=0;i<vName2.size();i++) {gT2->npar[i] = vPar2[i]; cout << vPar2[i] << " ";} cout<<endl;
  
  const int NumPar = gNPAR_FIT1+2*gAPAR_FIT1 + gNPAR_FIT2+2*gAPAR_FIT2;
  double final_pars[NumPar];
  for(int i=0;i<gNPAR_FIT1;i++) final_pars[i] = vPar1[i];
  for(int i=0;i<gAPAR_FIT1;i++) final_pars[gNPAR_FIT1+2*i] = 0;//rPar1[i];
  for(int i=0;i<gAPAR_FIT1;i++) final_pars[gNPAR_FIT1+2*i+1] = 0;//iPar1[i];
  for(int i=0;i<gNPAR_FIT2;i++) final_pars[gNPAR_FIT1+2*gAPAR_FIT1+i] = vPar2[i];
  for(int i=0;i<gAPAR_FIT2;i++) final_pars[gNPAR_FIT1+2*gAPAR_FIT1+gNPAR_FIT2+2*i] = 0;//rPar2[i];
  for(int i=0;i<gAPAR_FIT2;i++) final_pars[gNPAR_FIT1+2*gAPAR_FIT1+gNPAR_FIT2+2*i+1] = 0;//iPar2[i];
  final_pars[gNPAR_FIT1]=1.0;
  final_pars[gNPAR_FIT1+2*gAPAR_FIT1+gNPAR_FIT2]=1.0;
  for(int i=0;i<NumPar;i++) cout << final_pars[i] << " "; cout<<endl;

  cout << "------------------> Calculation of amplitude" << endl;
  TF1 *fInt1 = new TF1("int1",intensity1,0.5,2.5,NumPar); fInt1->SetParameters(final_pars);
  TF1 *fInt2 = new TF1("int2",intensity2,0.5,2.5,NumPar); fInt2->SetParameters(&final_pars[gNPAR_FIT1+2*gAPAR_FIT1]);
  TF1 *fPhi  = new TF1("phi" ,dphase    ,0.5,2.5,NumPar); fPhi ->SetParameters(final_pars);

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","plots",0,0,1200,400);
  can->Divide(3,1);
  can->cd(1); fInt1->Draw();
  can->cd(2); fInt2->Draw();
  can->cd(3);  fPhi->Draw();

  can->SaveAs("/tmp/1.pdf");

  return 0;
}

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
  ipar=0;
  while(1) {
    TString svName = TString::Format("ch2_rN%d",ipar);
    if( br->FindObject(svName.Data()) ) vName2.push_back(svName); else break;
    std::cout << ipar << ": found ch2 N" << std::endl;
    ipar++;
  }
  //set result
  gNPAR_FIT1 = vName1.size(); gNPAR_FIT2 = vName2.size();
  
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

  double vPar1[gNPAR_FIT1];
  double vPar2[gNPAR_FIT2];
  for(int i=0;i<vName1.size();i++) t->SetBranchAddress(vName1[i].Data(),&vPar1[i]);
  for(int i=0;i<vName2.size();i++) t->SetBranchAddress(vName2[i].Data(),&vPar2[i]);

  cout << "-------------> Creating models: " << gNPAR_FIT1 << "," << gNPAR_FIT2  << endl;
  const double eTH = 3*PI_MASS, sTH = pow(eTH,2);
  gT1 = new NoverD(gNPAR_FIT1,0,S1,rho,sTH,5000,500);
  gT2 = new NoverD(gNPAR_FIT2,0,S1,rho,sTH,5000,500);

  //fill the model
  t->GetEntry(emin);
  for(int i=0;i<vName1.size();i++) {gT1->npar[i] = vPar1[i]; cout << vPar1[i] << endl;}
  for(int i=0;i<vName2.size();i++) {gT2->npar[i] = vPar2[i]; cout << vPar2[i] << endl;}

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
      cd DI1   = gT1->DI(s);          fs1[i-1][j-1] = DI1;
      cd DII1  = DI1 - gT1->Disc(s);  ss1[i-1][j-1] = DII1;
      cd DI2   = gT2->DI(s);	      fs2[i-1][j-1] = DI2;
      cd DII2  = DI2 - gT2->Disc(s);  ss2[i-1][j-1] = DII2;
      if(i==61&&j==50) cout<<"s = "<<s<<": DI1 = "<<DI1<<", DI2 = "<<DI2<<endl;
    }
  }

  ofstream tfout;
  tfout.open ("/tmp/sheets.out");
  tfout << "-----------> FIRST  SHEET OF 1" << endl;
  printToFile(&fs1[0][0],t2d,tfout);
  tfout << "-----------> SECOND SHEET OF 1" << endl;
  printToFile(&ss1[0][0],t2d,tfout);
  tfout << "-----------> FIRST  SHEET OF 2" << endl;
  printToFile(&fs2[0][0],t2d,tfout);
  tfout << "-----------> SECOND SHEET OF 2" << endl;
  printToFile(&ss2[0][0],t2d,tfout);
  tfout.close();  

}


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int fit_data(int nAttempt, const char *fout_name, int pars_all) {

  gNPAR_FIT1 = (pars_all % 10000) / 1000;
  gAPAR_FIT1 = (pars_all % 1000 ) / 100;
  gNPAR_FIT2 = (pars_all % 100  ) / 10;
  gAPAR_FIT2 = (pars_all % 10   ) / 1;
  cout<<"Fit with model "<<gNPAR_FIT1<<gAPAR_FIT1<<gNPAR_FIT2<<gAPAR_FIT2<<endl;
  
  //**********************Seed**********************//
  gRandom->SetSeed(12); cout << "The seed " << gRandom->GetSeed() << " is used." << endl;

  //*****************Create the model***************//
  const double eTH = 3*PI_MASS, sTH = pow(eTH,2);
  gT1 = new NoverD(gNPAR_FIT1,0,S1,rho,sTH,10000,500);
  gT2 = new NoverD(gNPAR_FIT2,0,S1,rho,sTH,10000,500);

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
  //gh1->Scale(1./max);
  for(int i=1;i<=gh1->GetXaxis()->GetNbins();i++) { gh1->SetBinContent(i,gh1->GetBinContent(i)*(1./max)); gh1->SetBinError(i,gh1->GetBinError(i)*(1./max)); }

  max = gh2->GetBinContent(gh2->GetMaximumBin()); 
  //gh2->Scale(1./max);
  for(int i=1;i<=gh2->GetXaxis()->GetNbins();i++) { gh2->SetBinContent(i,gh2->GetBinContent(i)*(1./max)); gh2->SetBinError(i,gh2->GetBinError(i)*(1./max)); }

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","plots",0,0,1000,600);
  can->Divide(2,2);

  can->cd(1); gh1->SetStats(kFALSE);
  can->cd(4); gh2->SetStats(kFALSE);
  can->cd(2); gh3->SetStats(kFALSE);
  
  /* phase space */

  double gMth = 3*PI_MASS;
  double gSth = gMth*gMth;

  const int Nrho = 1000;
  const double Mhg = 3.0;
  for(int i=0;i<Nrho;i++) {
    double mass = gMth + (Mhg - gMth)/(Nrho-1)*i;
    glookup_rho_F2 .push_back( make_pair( mass, real(rho (mass*mass) )) );
  }

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
  
  const int NumPar = gNPAR_FIT1+2*gAPAR_FIT1 + gNPAR_FIT2+2*gAPAR_FIT2;
  
  // Create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor functor(&globChi2,NumPar);
  min->SetFunction(functor);
    
  const double step = 0.01;
  
  // Set the names
  string var_name[NumPar];
  prepareNames(gNPAR_FIT1,gAPAR_FIT1,&var_name[0],"ch1_");
  prepareNames(gNPAR_FIT2,gAPAR_FIT2,&var_name[gNPAR_FIT1+2*gAPAR_FIT1],"ch2_");
  // Set the limites for the variables
  double limit[NumPar];
  for(int i=0;i<gNPAR_FIT1;i++)      limit[i] = MRANGE_NPAR/abs(gT1->getmax_value_of_integral(i));
  for(int i=0;i<2*gAPAR_FIT1;i++)    limit[i+gNPAR_FIT1] = MRANGE_APAR*1.0;
  for(int i=0;i<gNPAR_FIT2;i++)      limit[i+gNPAR_FIT1+2*gAPAR_FIT1] = MRANGE_NPAR/abs(gT2->getmax_value_of_integral(i));
  for(int i=0;i<2*gAPAR_FIT2;i++)    limit[i+gNPAR_FIT2+gNPAR_FIT1+2*gAPAR_FIT1] = MRANGE_APAR*1.0;
  // Set the variables to be minimized
  for(int i=0;i<NumPar;i++) min->SetLimitedVariable(i,var_name[i],0.0,step,-limit[i],limit[i]);//
  for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,limit[i]/1e3);

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
    // starting point
    double start_pars[NumPar];
    int init_iterations = 0;
    chi2 = CUT_CHI2+1; double minimal_chi2 = 1e50; int best_seed = gRandom->GetSeed();
    while(init_iterations < MAX_INIT_ITERATION && chi2 > CUT_CHI2) {
      //if(init_iterations>=MAX_INIT_ITERATION-1) {cerr<<"Warning: failed to find good starting point. Best seed "<<best_seed<<" will be used. chi2 = " << minimal_chi2 << endl; gRandom->SetSeed(best_seed);}
      for(int i=0;i<NumPar;i++) start_pars[i] = limit[i]*(2*gRandom->Rndm()-1);
      start_pars[gNPAR_FIT1+1] = 0; min->FixVariable(gNPAR_FIT1+1);
      //double new_pars[] = {10,10,50,500,10,0,10,0,-10,0,1,0};
      //memcpy(start_pars,new_pars,12*sizeof(double));
      min->SetVariableValues(start_pars);
      chi2 = globChi2(start_pars); //if(chi2<minimal_chi2) {minimal_chi2=chi2; best_seed = gRandom->GetSeed();}
      init_iterations++;
    }
    cout<<"Initual parameters:\n";for(int i=0;i<NumPar;i++) {cout<<var_name[i]<<" = "<<start_pars[i]<<";\n";}

    //do the minimization
    min->Minimize(); //min->Hesse();
    bool fit_stat = 1;
    status = min->Status(); //if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);

    TF1 *fInt1 = new TF1("int1",intensity1,LEFT_INT1,RIGHT_INT1,NumPar); fInt1->SetParameters(final_pars);
    TF1 *fInt2 = new TF1("int2",intensity2,LEFT_INT2,RIGHT_INT2,NumPar); fInt2->SetParameters(&final_pars[gNPAR_FIT1+2*gAPAR_FIT1]);
    TF1 *fPhi  = new TF1("phi" ,dphase,    LEFT_PHI ,RIGHT_PHI ,NumPar); fPhi ->SetParameters(final_pars);

//    can->cd(1); gh1->Draw(); fInt1->Draw("same");
//    can->cd(4); gh2->Draw(); fInt2->Draw("same");
//    can->cd(2); gh3->Draw();  fPhi->Draw("same");
    can->cd(1); fInt1->Draw();
    can->cd(4); fInt2->Draw();
    can->cd(2);  fPhi->Draw();
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

  vector<double> npar1(gNPAR_FIT1), npar2(gNPAR_FIT2);  
  vector<cd>     apar1(gAPAR_FIT2), apar2(gAPAR_FIT2);
  fillvectors(&par[0],apar1,npar1);
  fillvectors(&par[gNPAR_FIT1+2*gAPAR_FIT1],apar2,npar2);

  // *********************** Calculate chi2 ********************** //
  //Intensity1
  double chi0 = 0;
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT1 || mass>RIGHT_INT1) continue;
    double s = pow(mass,2);
    cd amp = Afull(gT1,s,apar1,npar1,S2);
    double ph = MStructureHolder::getphsp(s,glookup_rho_F2);
    double intens = norm(amp)*ph;//gT1->getrho(s);
    double err = gh1->GetBinError(i);
//    cout << "M = " << mass << ", A = " << amp 
//	 << ", C = " << setprecision(10) << gh1->GetBinContent(i) 
//	 << ", E = " << setprecision(10) << gh1->GetBinError(i) 
//	 << ", H = " << std::setprecision(10) << ph
//	 << ", P = " << setprecision(10) << pow(intens - gh1->GetBinContent(i),2)/(err*err) << endl;
    chi0 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //Intensity2
  double chi1 = 0;
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double mass = gh2->GetBinCenter(i);
    if(mass<LEFT_INT2 || mass>RIGHT_INT2) continue;
    double s = pow(mass,2);
    cd amp = Afull(gT2,s,apar2,npar2,S2);
    double ph = MStructureHolder::getphsp(s,glookup_rho_F2);
    double intens = norm(amp)*ph;//gT2->getrho(s);
    double err = gh2->GetBinError(i);
//    cout << "M = " << mass << ", A = " << amp 
//	 << ", C = " << setprecision(10) << gh2->GetBinContent(i) 
//	 << ", E = " << setprecision(10) << gh2->GetBinError(i)
//	 << ", H = " << std::setprecision(10) << ph
//	 << ", P = " << setprecision(10) << pow(intens - gh2->GetBinContent(i),2)/(err*err) << endl;
    chi1 += pow(intens - gh2->GetBinContent(i),2)/(err*err);
  }
  //Phase
  double chiPhi = 0;
  for(int i=1;i<=gh3->GetNbinsX();i++) {
    double mass = gh3->GetBinCenter(i);
    if(mass<LEFT_PHI || mass>RIGHT_PHI) continue;
    double s = pow(mass,2);
    cd amp1 = Afull(gT1,s,apar1,npar1,S2);
    cd amp2 = Afull(gT2,s,apar2,npar2,S2);
    double phase = normPhase( arg(amp1*conj(amp2)) );
    double err = gh3->GetBinError(i);
    double dPhi = normPhase(phase-gh3->GetBinContent(i));
    //std::cout << value << std::endl;    
    chiPhi += pow(dPhi,2)/(err*err);
  }
  double chi2 = chi0+chi1+chiPhi;
  cout << "chi2 = " << setw(16) << setprecision(15) << chi0;
  cout << ", second chi2 = " << setw(16) << setprecision(15) << chi1;
  cout << ", last chiPhi = " << setw(16) << setprecision(15) << chiPhi;
  cout << ", final chi2 = "  << setw(16) << setprecision(15) << chi2 << endl;
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
  vector<double> npar(gNPAR_FIT1);
  vector<cd> apar(gAPAR_FIT1);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(gT1,x[0]*x[0],apar,npar,S2);
  return norm(amp)*gT1->getrho(x[0]*x[0]);
}

double intensity2(double *x, double *par) {
  vector<double> npar(gNPAR_FIT2);
  vector<cd> apar(gAPAR_FIT2);
  fillvectors(&par[0],apar,npar);
  cd amp = Afull(gT1,x[0]*x[0],apar,npar,S2);
  return norm(amp)*gT2->getrho(x[0]*x[0]);
}

double dphase(double *x, double *par) {
  vector<double> npar1(gNPAR_FIT1), npar2(gNPAR_FIT1);  
  vector<cd>     apar1(gAPAR_FIT2), apar2(gAPAR_FIT2);
  fillvectors(&par[0],apar1,npar1);
  fillvectors(&par[gNPAR_FIT1+2*gAPAR_FIT1],apar2,npar2);
  cd amp1 = Afull(gT1,x[0]*x[0],apar1,npar1,S2);
  cd amp2 = Afull(gT2,x[0]*x[0],apar2,npar2,S2);
  return normPhase( arg(amp1*conj(amp2)) );
}

void fillvectors(const double *par, vector<cd> &apar, vector<double> &npar) {
  //Constract Numerator
  int np = npar.size();
  for(int i=0;i<np;i++) npar[i] = par[i];
  //Constract Alpha
  for(int i=0;i<apar.size();i++) apar[i] = cd(par[np+2*i],par[np+2*i+1]);
}

double normPhase(double a) {
  if(a> M_PI)  return normPhase(a-2*M_PI);
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
 
