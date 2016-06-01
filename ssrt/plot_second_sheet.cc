#include <iostream>
#include <sstream>
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
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

#define NPAR_FIT  3
#define APAR_FIT  1

//#define S1 -10.0
#define S2 0.5

#define LEFT_INT   RHO_MASS+PI_MASS
#define RIGHT_INT  1.5

#define VFUNC_NPOINT 100
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
cd AfullI(cd s, vector<cd> apar, vector<cd> npar, double s1, double s2);
cd DfullI (cd s, vector<cd> npar, double s1);
cd DfullII(cd s, vector<cd> npar, double s1);
double globChi2(const double *par);

void fillvectors(const double *par, vector<cd> &apar, vector<cd> &npar);
double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);//*int(( a+M_PI)/(2*M_PI));
  if(a<-M_PI) return normPhase(a+2*M_PI);//a+2*M_PI*int((-a-M_PI)/(2*M_PI));
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
double normDII(const double *pars);
////////////////////////////////////////////

vector<vector<pair<double,cd> > > *garr;
TH1D *gh1, *gh2;
double gs1;
vector<cd> *gnpar;

int fit_data(int nAttempt, const char *fout_name);
int find_poles(int nAttempt, const char *fin_name, const char *fout_name);
int plot_best_result(const char *fin_name);

int main(int ac, char **av) {

  if(ac<2) { cerr << "first arg is FIT or FIND or PLOT" << endl; return 1; }

  const int nAttempt = (ac > 2) ? atoi(av[2]) : 1;
  const char *file_name = (ac > 3) ? av[3] : "/tmp/test.root";

  if(strcmp(av[1],"FIT")==0) {
    return fit_data(nAttempt,file_name);
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


int find_poles(int nAttempt, const char *fin_name, const char *fout_name) {

  const int NumPar = 2*NPAR_FIT + 2*APAR_FIT + 1;

  TFile *fin = TFile::Open(fin_name);
  TTree *t = (TTree*)gDirectory->Get("mins");
  if(!t) {cerr<<"no tree"<<endl; return 1;}

  //*********** Find entry with smallest chi2 **********//
  const int Nentries = t->GetEntries(0);
  double chi2; t->SetBranchAddress("chi2",&chi2);
  t->GetEntry(0); double minChi2=chi2; int minEntry=0;
  for(int i=1;i<Nentries;i++) {
    t->GetEntry(i);
    if(chi2<minChi2) { minChi2=chi2; minEntry=i;}
  }
  cout<<"min chi2 is "<<minChi2<<", entry is "<<minEntry<<endl;

  //******************* Set branches ********************//
  string var_name[NumPar];
  for(int i=0;i<NPAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rN"<<i; iIO<<"iN"<<i; var_name[2*i] = rIO.str(); var_name[2*i+1] = iIO.str();}
  for(int i=0;i<APAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rA"<<i; iIO<<"iA"<<i; var_name[2*NPAR_FIT+2*i] = rIO.str(); var_name[2*NPAR_FIT+2*i+1] = iIO.str();}
  var_name[NumPar-1] = "sn1"; 
  double par[NumPar];
  for(int i=0;i<NumPar;i++) t->SetBranchAddress(var_name[i].c_str(),&par[i]);
  t->GetEntry(minEntry);

  vector<cd> npar(NPAR_FIT), apar(APAR_FIT);
  fillvectors(&par[0],apar,npar);  
  double s1 = par[NumPar-1];
  gs1 = s1;
  gnpar = &npar;
  //cout << "gs1 = " << gs1 << endl;

  //**************** Calculate integrals ****************//
  vector<vector<pair<double,cd> > > data;
  build_integrals_table(data, s1);
  garr=&data;
  cout << "test access " << (*garr)[0][5].second << endl;
  
  //Build minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  //min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.01);
  min->SetPrintLevel(1);
  
  // create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor f(&normDII,2); 
  min->SetFunction(f);
    
  double step = 0.01;
  // Set the free variables to be minimized!
  min->SetVariable(0,"rS",0, step); 
  min->SetVariable(1,"iS",0, step); min->SetVariableUpperLimit(1,0);

  //Do minimisation several times
  TFile fout(fout_name,"RECREATE");
  TTree tp("pols","several_pols");
  double final_pos[2];
  tp.Branch("rS",&final_pos[0]);
  tp.Branch("iS",&final_pos[1]);
  double v_normDII; tp.Branch("normDII",&v_normDII);
  int status; tp.Branch("status",&status);
  for(int e=0;e<nAttempt;e++) {
    //reload step size
    min->SetVariableStepSize(0,step);
    min->SetVariableStepSize(1,step);
    //set parameters
    min->SetVariableValue(0,gRandom->Uniform(0.5,2.5));
    min->SetVariableValue(1,gRandom->Uniform(-0.5,0 ));
    cout<<"------start in ("<<min->X()[0]<<", "<<min->X()[1]<<")------"<<endl;

    //do minimization
    min->SetStrategy(1); min->Minimize(); 
    min->SetStrategy(2);  bool fit_stat = min->Minimize(); 
    status = min->Status(); if(status > 1) continue;
    memcpy(final_pos,min->X(),2*sizeof(double));
    v_normDII = normDII(final_pos);

    //store results
    tp.Fill();
  }  

  tp.Write();
  fout.Close();

  return 0;
}

double normDII(const double *pars) {
  cd s(pars[0],pars[1]);
  cd DII = DfullII(s,*gnpar,gs1);
  double nDII = norm(DII);
  //cout <<"s = "<<pars[0]<<"+I "<<pars[1]<<", rN0 = "<<(*gnpar)[0]<<", DII = " << nDII << endl;
  return nDII;
}

int plot_best_result(const char *fin_name) {
  const int NumPar = 2*NPAR_FIT + 2*APAR_FIT + 1;

  TFile *fin = TFile::Open(fin_name);
  TTree *t = (TTree*)gDirectory->Get("mins");
  if(!t) {cerr<<"no tree"<<endl; return 1;}

  //*********** Find entry with smallest chi2 **********//
  const int Nentries = t->GetEntries();
  double chi2; t->SetBranchAddress("chi2",&chi2);
  t->GetEntry(0); double minChi2=chi2; int minEntry=0;
  for(int i=1;i<Nentries;i++) {
    t->GetEntry(i);
    if(chi2<minChi2) { minChi2=chi2; minEntry=i;}
  }
  cout<<"min chi2 is "<<minChi2<<", entry is "<<minEntry<<endl;

  //******************* Set branches ********************//
  string var_name[NumPar];
  for(int i=0;i<NPAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rN"<<i; iIO<<"iN"<<i; var_name[2*i] = rIO.str(); var_name[2*i+1] = iIO.str();}
  for(int i=0;i<APAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rA"<<i; iIO<<"iA"<<i; var_name[2*NPAR_FIT+2*i] = rIO.str(); var_name[2*NPAR_FIT+2*i+1] = iIO.str();}
  var_name[NumPar-1] = "sn1"; 
  double par[NumPar];
  for(int i=0;i<NumPar;i++) t->SetBranchAddress(var_name[i].c_str(),&par[i]);
  t->GetEntry(minEntry);

  vector<cd> npar(NPAR_FIT), apar(APAR_FIT);
  fillvectors(&par[0],apar,npar);  
  double s1 = par[NumPar-1];

  //**************** Calculate integrals ****************//
  vector<vector<pair<double,cd> > > data;
  build_integrals_table(data, s1);
  garr=&data;
  cout << "test access " << (*garr)[0][5].second << endl;

  //**************** Something else *********************//
  const int nPoint(100);
  double mass[nPoint];
  double re[nPoint];
  double im[nPoint];
  double sigma[nPoint];
  double mf = RHO_MASS+PI_MASS;
  double ml = 1.4;
  double step = (ml-mf)/(nPoint-1);
  for(int i=0;i<nPoint;i++) {
    mass[i] = mf+step*i;
    cd amp = Afull(mass[i]*mass[i],apar,npar,s1,S2);
    re[i] = real(amp);
    im[i] = imag(amp);
    sigma[i] = norm(amp)*rho(mass[i],RHO_MASS,PI_MASS);
  }

  TGraph gre(nPoint,mass,re); gre.SetTitle("Real part;mass");
  TGraph gim(nPoint,mass,im); gim.SetTitle("Imag part;mass");
  TGraph gsig(nPoint,mass,sigma); gsig.SetTitle("Cross sectrion;mass");

  TCanvas c1("c1","canvas",0,0,1000,1000);
  c1.Divide(3,3);
  c1.cd(1); gre.Draw("apl");
  c1.cd(2); gim.Draw("apl");
  c1.cd(3); gsig.Draw("apl");

  //Plot first and second sheets
  const int Nx(100); const double xf(0.5) ,xl(1.5);
  const int Ny(100); const double yf(-0.5),yl(0.5);

  TH2D hreI  ("hreI"  ,"Real part, I sheet;mass;Im s"    ,Nx,xf,xl,Ny,yf,yl);
  TH2D himI  ("himI"  ,"Imag part, I sheet;mass;Im s"    ,Nx,xf,xl,Ny,yf,yl);
  TH2D habsI ("habsI" ,"Norm |A|^{2}, I sheet;mass;Im s" ,Nx,xf,xl,Ny,yf,yl);
  TH2D hreII ("hreII" ,"Real part, II sheet;mass;Im s"   ,Nx,xf,xl,Ny,yf,yl);
  TH2D himII ("himII" ,"Imag part, II sheet;mass;Im s"   ,Nx,xf,xl,Ny,yf,yl);
  TH2D habsII("habsII","Norm |A|^{2}, II sheet;mass;Im s",Nx,xf,xl,Ny,yf,yl);
  
  for(int i=1;i<=Nx;i++) {
    for(int j=1;j<=Ny;j++) {
      cd s(pow(hreI.GetXaxis()->GetBinCenter(i),2),hreI.GetYaxis()->GetBinCenter(j));
      cd DI  = DfullI (s,npar,s1);
      cd DII = DfullII(s,npar,s1);
      int bin = hreI.GetBin(i,j);
      hreI  .SetBinContent(bin,real(DI ));
      himI  .SetBinContent(bin,imag(DI ));
      habsI .SetBinContent(bin, abs(DI ));
      hreII .SetBinContent(bin,real(DII));
      himII .SetBinContent(bin,imag(DII));
      habsII.SetBinContent(bin, abs(DII));
    }
  }

  c1.cd(4); hreI.  Draw("colz");
  c1.cd(5); himI.  Draw("colz");
  c1.cd(6); habsI. Draw("colz");
  c1.cd(7); hreII. Draw("colz");
  c1.cd(8); himII. Draw("colz");
  c1.cd(9); habsII.Draw("colz");

  c1.SaveAs("c2.png");

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int fit_data(int nAttempt, const char *fout_name) {

  //**************** Calculate integrals ****************//
  vector<vector<pair<double,cd> > > data;

  gs1 = -1;
  build_integrals_table(data,gs1);

  garr=&data;
  cout << "test access " << (*garr)[0][5].second << endl; 
  
  //**************** Load the data *****************//
  double dataLeft = RHO_MASS+PI_MASS;  double dataRight = 1.5; double dataNbin = 100;
  TH1D *hist_data = new TH1D("data","Intensity;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  TH1D *hist_data_phase = new TH1D("phase","Phase;M_{3#pi}^{2}",dataNbin,dataLeft,dataRight);
  double sigma = 0.02*norm(BreitWignerA(A1_MASS,A1_MASS,A1_WIDTH));
  for(int i=1;i<=hist_data->GetNbinsX();i++) {
    hist_data->SetBinContent(i,
			     gRandom->Gaus(norm(BreitWignerA(hist_data->GetBinCenter(i),
						       A1_MASS,A1_WIDTH)),
				     sigma)
			     );
    hist_data->SetBinError(i,sigma);
  }
  double sigmaphi = 0.02*2*M_PI;
  for(int i=1;i<=hist_data->GetNbinsX();i++) {
    hist_data_phase->SetBinContent(i,
				   gRandom->Gaus(arg(BreitWignerA(hist_data->GetBinCenter(i),
								  A1_MASS,A1_WIDTH)),
						 sigmaphi)
				   );
    hist_data_phase->SetBinError(i,sigmaphi);

  }
  
  gh1 = hist_data;
  gh2 = hist_data_phase; //gh2->GetYaxis()->SetRangeUser(-M_PI,M_PI);

  TCanvas *can = new TCanvas("c1","plots",0,0,1000,600);
  can->Divide(2,1);
 
  //***************** Fit the data *****************//
  //Build minimizer
  ROOT::Math::Minimizer* min = 
     ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");

  // set tolerance , etc...
  min->SetMaxFunctionCalls(100000); // for Minuit/Minuit2 
  min->SetTolerance(0.1);
  ROOT::Math::MinimizerOptions mopt; 
  min->SetStrategy(2);
  min->SetPrintLevel(1);
  
  const int NumPar = 2*NPAR_FIT + 2*APAR_FIT + 1;
  
  // create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor f(&globChi2,NumPar);
  min->SetFunction(f);
    
  const double step = 0.01;
  const double limit = 30000;
  //form names
  string var_name[NumPar];
  for(int i=0;i<NPAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rN"<<i; iIO<<"iN"<<i; var_name[2*i] = rIO.str(); var_name[2*i+1] = iIO.str();}
  for(int i=0;i<APAR_FIT;i++) { ostringstream rIO,iIO; rIO<<"rA"<<i; iIO<<"iA"<<i; var_name[2*NPAR_FIT+2*i] = rIO.str(); var_name[2*NPAR_FIT+2*i+1] = iIO.str();}
  var_name[NumPar-1] = "sn1"; 
  // Set the free variables to be minimized!
  for(int i=0;i<NumPar;i++) min->SetLimitedVariable(i,var_name[i],0.0,step,-limit,limit);
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
    for(int i=0;i<NumPar;i++) { min->SetVariableStepSize(i,step);/*min->SetVariableUpperLimit(i,10000); min->SetVariableLowerLimit(i,-10000);*/}

    min->SetVariableStepSize(NumPar-1,0.1);
    // starting point
    double start_par[NumPar];
    gRandom->SetSeed(0);
    for(int i=0;i<2*NPAR_FIT;i++) if(i%2==0) start_par[i] = 10*(2*gRandom->Rndm()-1); else start_par[i] = 0;
    for(int i=0;i<2*APAR_FIT;i++) if(i%2==0) start_par[2*NPAR_FIT+i] = gRandom->Rndm(); else start_par[2*NPAR_FIT+i] = 0;
    //start_par[0]=0.194764; start_par[2]=8.36437; start_par[4]=46.3398; start_par[6]=1.;
    start_par[0]=10; start_par[2]=0; start_par[4]=0; start_par[6]=1.;
    //start_par[NumPar-1] = gRandom->Uniform(-1000,pow(RHO_MASS+PI_MASS,2));//pow(RHO_MASS+PI_MASS,2)-1./
    start_par[NumPar-1] = gs1;
    min->SetVariableValues(start_par);
    min->FixVariable(NumPar-1);
    min->FixVariable(NumPar-2);
    for(int i=0;i<NPAR_FIT;i++) min->FixVariable(2*i+1);

    //do the minimization
    bool fit_stat = 1;//min->Minimize(); //min->Hesse();
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
    double err = gh1->GetBinError(i);
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //Phase
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT || mass>RIGHT_INT) continue;
    double s = pow(mass,2);
    cd amp = Afull(s,apar,npar,par[2*NPAR_FIT+2*APAR_FIT],S2);
    double phase = arg(amp);
    double err = gh2->GetBinError(i);
    double dPhi = normPhase(phase-gh2->GetBinContent(i));
    chi2 += pow(dPhi,2)/(err*err);
  }
  //  cout << "par[2] = "<< par[2] <<", chi2 = " << chi2<<endl;
  return chi2;
}

void build_integrals_table(vector<vector<pair<double,cd> > > &data, double s1) {
  data.clear();
  for(int i=0;i<NPAR_FIT;i++) {
    vector<pair<double,cd> > vj;
    double step = 1.0*(VFUNC_RB - VFUNC_LB)/(VFUNC_NPOINT-1);
    for(int j=0;j<VFUNC_NPOINT;j++) {
      double s = VFUNC_LB + step*j;
      cd value = vf(s,i,s1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));
      cout << "s = " << s << ", v = " << value << endl;
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
  double p  = sqrt(LAMBDA(x*x,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*x);
  double p0 = sqrt(LAMBDA(m*m,PI_MASS*PI_MASS,PI_MASS*PI_MASS))/(2*m);
  double R = 5; double BW = (1+pow(R*p0,2))/(1+pow(R*p,2));
  double G = G0*rho(x,RHO_MASS,PI_MASS)/rho(m,RHO_MASS,PI_MASS);//*pow(p/p0,2)*BW

  cd A = sqrt(2*m*G)/(m*m-x*x-cd(0,1)*m*G);
  return A;
} 

cd expand(cd s, cd s1, vector<cd> par) {
  //cout << real(par[0]) << " " << imag(par[0]) << " " << norm(par[0]) << endl; 
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
  if(s1!=gs1) { gs1=s1; build_integrals_table(*garr, s1); } /*cout<<"------------->--rebuild all integrals, s1="<<s1<<"--<---------------"<<endl;*/ 
  for(int i=0;i<npar.size();i++) dsum += npar[i]*vtable(real(s),(*garr)[i]);//vf(s,i,s1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));//
  cd Denomin = 1. - s / (2*M_PI) * dsum;

  //Constract N, Alpha and A.
  cd Numir = expand(s,s1,npar);
  cd Aprod = expand(s,s2,apar);
  cd Afull = Aprod*Numir/Denomin;

  return Afull;
}

cd AfullI(cd s, vector<cd> apar, vector<cd> npar, double s1, double s2) {

  //Constract D
  cd Denomin = DfullI(s,npar,s1);

  //Constract N, Alpha and A.
  cd Numir = expand(s,s1,npar);
  cd Aprod = expand(s,s2,apar);
  cd Afull = Aprod*Numir/Denomin;

  return Afull;
}

cd DfullI(cd s, vector<cd> npar, double s1) {
  cd dsum =0;
  for(int i=0;i<npar.size();i++) dsum += npar[i]*vf(s,i,s1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));
  cd DfullI = 1. - s / (2*M_PI) * dsum;
  return  DfullI;
}

cd DfullII(cd s, vector<cd> npar, double s1) {

  //Constract D
  cd DI = DfullI(s,npar,s1);
  //Construct N
  cd Numir = expand(s,s1,npar);
  //Construct Rho
  double sl = pow(RHO_MASS-PI_MASS,2);
  double sr = pow(RHO_MASS+PI_MASS,2);  
  cd rho = sqrt( (s - sl)*(s -sr)/(s*s) * exp(-cd(0,1)*M_PI) ) * exp(cd(0,1)*M_PI/2.);  
  
  cd DII = DI + cd(0,1)*rho*Numir;
  return DII;
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
  cd alog = 1.-(s+cd(0,1)*EPSILON)/sr;
  cd first_int = integral(s,k,s0,sl,sr);
  return first_int-mult*log(alog);
}


    //min->Minimize();
    //min->Hesse();
    //min->FixVariable(4); min->FixVariable(5);
    //for(int i=2;i<2*NPAR_FIT;i++) min->SetVariableValue(i,start_par[i]);
    //min->Minimize();
    //min->ReleaseVariable(4); min->ReleaseVariable(5);
    //for(int i=0;i<2*NPAR_FIT;i++) min->FixVariable(i);
    //min->Minimize();
    //for(int i=0;i<NPAR_FIT;i++) min->ReleaseVariable(2*i);
    //
    //for(int i=1;i<NPAR_FIT;i++) min->FixVariable(2*i);
    //min->Minimize();
    //for(int i=0;i<NPAR_FIT;i++) min->ReleaseVariable(2*i);

    //min->Minimize();
    
//    min->ReleaseVariable(NumPar-1); 
//    for(int i=0;i<NumPar-1;i++) min->FixVariable(i);     
//    min->Minimize(); 
//    
//    min->FixVariable(NumPar-1);     
//    for(int i=0;i<NumPar-1;i++) min->ReleaseVariable(i);
//    for(int i=0;i<NPAR_FIT;i++) min->FixVariable(2*i+1);
//    min->Minimize(); 
//
//    min->ReleaseVariable(NumPar-1);
    







//  /***********************************************************/
//  vector<cd> npar(NPAR_FIT), apar(APAR_FIT);
//  for(int i=0;i<NumPar;i++) final_pars[i] = 0;
//  final_pars[0]=10; final_pars[2]=0; final_pars[4]=0; final_pars[6]=1.;
//  fillvectors(&final_pars[0],apar,npar);
//  cd s(1.1,0);
//  cout << "A = " << Afull(s,apar,npar,gs1,0) << endl;
//  cout << "N = " << expand(s,gs1,npar) << endl;
//  /***********************************************************/
////cd Afull(cd s, vector<cd> apar, vector<cd> npar, 
////	 double s1, double s2) {
////
////  //Constract D
//  cd dsum =0;
////  if(s1!=gs1) { gs1=s1; build_integrals_table(*garr, s1); } /*cout<<"------------->--rebuild all integrals, s1="<<s1<<"--<---------------"<<endl;*/ 
//  for(int i=0;i<npar.size();i++) dsum += npar[i]*vtable(real(s),(*garr)[i]);//vf(s,i,s1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));//
//  cd Denomin = 1. - s / (2*M_PI) * dsum;
//  cout << "D = " << Denomin << endl;
////
////  //Constract N, Alpha and A.
////  cd Numir = expand(s,s1,npar);
////  cd Aprod = expand(s,s2,apar);
////  cd Afull = Aprod*Numir/Denomin;
////
////  return Afull;
////}
////



