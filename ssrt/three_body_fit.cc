//root
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
#include "Minuit2/Minuit2Minimizer.h"

//my
#include "NoverD.h"
#include "MIsobar.h"
#include <constants.h>

//model settings
#define NPAR_FIT1  3
#define APAR_FIT1  1
#define NPAR_FIT2  3
#define APAR_FIT2  1

#define S1 -1.0
#define S2 -2.0

#define LIMIT_TAKEN_AT 3.0
#define LIMIT_ACOEFF 5.
#define LIMIT_NCOEFF 30.

//region of fit
#define LEFT_INT1   (3*PI_MASS)
#define RIGHT_INT1  2.5
#define LEFT_INT2   (3*PI_MASS)
#define RIGHT_INT2  2.5
#define LEFT_PHI    1.5
#define RIGHT_PHI   2.5

using namespace std; 
typedef std::complex<double> cd;

NoverD *a1,*a2;
TH1D *gh1, *gh2, *gh3;
///////////////////////////////////////////
double globChi2(const double *par);
cd rho3_f2_S(cd s);
cd rho3_f2_D(cd s);
///////////////////////////////////////////
double intensity1(double *x, double *par);
double intensity2(double *x, double *par);
double dphase(double *x, double *par);
///////////////////////////////////////////
double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);
  if(a<-M_PI) return normPhase(a+2*M_PI);
  return a;
}
void prepareNames(int npar, int apar, string *var_name,const char *pref);
void fillvectors(const double *par, vector<double> &npar, vector<cd> &apar);
///////////////////////////////////////////

int main(int ac, char** av) {  

  const int nAttempt = (ac > 1) ? atoi(av[1]) : 1;
  const char *fout_name = (ac>2) ? av[2] : "/tmp/test.root";
  
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

  //Normalise
  double Smax = f2piS_int->GetBinContent( f2piS_int->GetMaximumBin() );
  double Dmax = f2piD_int->GetBinContent( f2piD_int->GetMaximumBin() );
  for(int i=1;i<=nBins;i++) {
    f2piS_int->SetBinContent(i,f2piS_int->GetBinContent(i)/Smax); f2piS_int->SetBinError(i,f2piS_int->GetBinError(i)/Smax);
    f2piD_int->SetBinContent(i,f2piD_int->GetBinContent(i)/Dmax); f2piD_int->SetBinError(i,f2piD_int->GetBinError(i)/Dmax);
  }


  gh1 = f2piS_int;
  gh2 = f2piD_int;
  gh3 = f2pi_dphi;

  TCanvas *can = new TCanvas("c1","plots",0,0,1000,600);
  can->Divide(2,2);
 
  can->cd(1); gh1->Draw();
  can->cd(4); gh2->Draw();
  can->cd(2); gh3->Draw();

  //**************** Build the model ***************//
  NoverD ndS(NPAR_FIT1,S1,rho3_f2_S,pow(3*PI_MASS,2),300,700);
  NoverD ndD(NPAR_FIT2,S1,rho3_f2_D,pow(3*PI_MASS,2),300,700);
  a1 = &ndS;
  a2 = &ndD;
  //cout << a1->A(pow(24,2)) << endl;
  //return 0;

  //***************** Fit the data *****************//
  //Build minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2","Combined");//Combined
  ROOT::Math::MinimizerOptions mopt = min->Options();
  mopt.SetMinimizerAlgorithm("Migrad");
  min->SetOptions(mopt);
  min->Options().Print();
//  ROOT::Minuit2::Minuit2Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer("scan");
//  min->Options().Print();

  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetTolerance(0.001);
  min->SetPrintLevel(2);
  min->SetStrategy(2);

  const int NumPar = 
    NPAR_FIT1 + 2*APAR_FIT1 +
    NPAR_FIT2 + 2*APAR_FIT2;
  
  // create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor f(&globChi2,NumPar);
  min->SetFunction(f);

  double step = 1.0;
  //form names
  string var_name[NumPar];
  prepareNames(NPAR_FIT1,APAR_FIT1,&var_name[0],"ch1_");
  prepareNames(NPAR_FIT2,APAR_FIT2,&var_name[NPAR_FIT1+2*APAR_FIT1],"ch2_");
  // Calculate limits
  double limit[NumPar];
  for(int i=0;i<NPAR_FIT1;i++) limit[i] = LIMIT_NCOEFF / abs(a1->getvi(LIMIT_TAKEN_AT,i));
  for(int i=0;i<2*APAR_FIT1;i++) limit[NPAR_FIT1+i] = LIMIT_ACOEFF * sqrt(gh1->GetBinContent(gh1->GetMaximumBin()));
  for(int i=0;i<NPAR_FIT2;i++) limit[NPAR_FIT1+2*APAR_FIT1+i] = LIMIT_NCOEFF / abs(a2->getvi(LIMIT_TAKEN_AT,i));
  for(int i=0;i<2*APAR_FIT2;i++) limit[NPAR_FIT1+2*APAR_FIT1+NPAR_FIT2+i] = LIMIT_ACOEFF * sqrt(gh2->GetBinContent(gh2->GetMaximumBin()));
  cout << "Limits: ";
  for(int i=0;i<NumPar;i++) cout << limit[i] << " ";
  cout << endl;		      
  // Set the free variables to be minimized
  for(int i=0;i<NumPar;i++) min->SetLimitedVariable(i,var_name[i],0, 
						    limit[i]/200.,
						    -limit[i], limit[i]);  

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

    //------------------- starting point -----------------------//
    
    double start_par[NumPar];
    gRandom->SetSeed(6); cout << "------------- Seed is set to "<<gRandom->GetSeed()<<" ----------------"<<endl;
    for(int i=0;i<NumPar;i++) if(i%2==0) start_par[i] = limit[i]*(2*gRandom->Rndm()-1); else start_par[i] = 0;
    start_par[NPAR_FIT1+1]=0;  //remove umbig phase  
    min->SetVariableValues(start_par);
    min->FixVariable(NPAR_FIT1+1);  //fix one phase
    /**/
    //setParameters
    /*
    double start_par2[NumPar] = {243.71186, 1082.8433, -1081.933, -5000, -26.79075, 0, 284.05362, -65.50888, 553.57364, -556.8757, -1765.791, -282.2118, 91.899669, 833.43312, 2007.2157, 675.12612, 242.49648, -214.0594, 1068.1972, -1068.197, 734.61093, -863.4845, -1068.197, 1068.1972};
    min->SetVariableValues(start_par);
    memcpy(start_par,start_par2,NumPar*sizeof(double));
    /**/

    //Scheam of fixation and releasation parameters    
//    for(int i=NPAR_FIT1+2*APAR_FIT1;i<NumPar;i++) min->FixVariable(i);
//    min->Minimize();
//    for(int i=NPAR_FIT1+2*APAR_FIT1;i<NumPar;i++) min->ReleaseVariable(i);
//    for(int i=0;i<NPAR_FIT1+2*APAR_FIT1;i++) min->FixVariable(i);
//    min->Minimize();
//    for(int i=0;i<NPAR_FIT1+2*APAR_FIT1;i++) if(i!=NPAR_FIT1+1) min->ReleaseVariable(i);
    bool fit_stat = min->Minimize();
    /**/
    //done
    
    status = min->Status();
    //if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);
    for(int i=0;i<NumPar;i++) cout << min->VariableName(i) << " = " << final_pars[i] << endl;
    
    TF1 fInt1("int1",intensity1,LEFT_INT1,RIGHT_INT1,NPAR_FIT1+2*APAR_FIT1); fInt1.SetParameters(final_pars);
    TF1 fInt2("int2",intensity2,LEFT_INT2,RIGHT_INT2,NPAR_FIT2+2*APAR_FIT2); fInt2.SetParameters(&final_pars[NPAR_FIT1+2*APAR_FIT1]);
    TF1 fPhi("dphi",dphase,     LEFT_PHI, RIGHT_PHI, NumPar);                 fPhi.SetParameters(final_pars);
    //fInt2.SetNpx(1000);
    //for(int i=0;i<100;i++) cout << abs(a2->A(pow(1+0.002*i,2))) << " " << abs(a2->AI(pow(1+0.002*i,2))) << endl;

    double mass_r = gh1->GetBinLowEdge(gh1->GetNbinsX())+gh1->GetBinWidth(1);
    //TF1 fInt1_fr(fInt1); fInt1_fr.SetRange(F2_MASS+PI_MASS,mass_r); fInt1_fr.SetLineColor(kOrange); fInt1_fr.SetParameters(final_pars);
    //TF1 fInt2_fr(fInt2); fInt2_fr.SetRange(F2_MASS+PI_MASS,mass_r); fInt2_fr.SetLineColor(kOrange); fInt2_fr.SetParameters(&final_pars[NPAR_FIT1+2*APAR_FIT1 ]);
    TF1  fPhi_fr("dphi_fr",dphase,(3*PI_MASS), RIGHT_PHI, NumPar);                 fPhi.SetParameters(final_pars); fPhi_fr.SetLineColor(kOrange);  fPhi_fr.SetParameters(final_pars);//(fPhi);   fPhi_fr.SetRange(F2_MASS+PI_MASS,mass_r);  
    can->cd(1);/* fInt1_fr.Draw("same");*/ fInt1.Draw("same");
    can->cd(2);  fPhi_fr.Draw("same");  fPhi.Draw("same");
    can->cd(4);/* fInt2_fr.Draw("same");*/ fInt2.Draw("same");
    
    if(nAttempt==1) can->SaveAs("c1.png");
    t.Fill();
    
  }

  t.Write();
  fout.Close();
  
  cout << "finished" << endl;
  return 0;
}



double globChi2(const double *par) {

  vector<cd> apar1(APAR_FIT1); fillvectors(&par[0]                    ,a1->npar,apar1);
  vector<cd> apar2(APAR_FIT2); fillvectors(&par[NPAR_FIT1+2*APAR_FIT1],a2->npar,apar2);

  // *********************** Calculate chi2 ********************** //
  double chi2 = 0;

  //Intensity1
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT1 || mass>RIGHT_INT1) continue;
    double s = pow(mass,2);
    cd amp = a1->ExtendA(s,S2,apar1);    
    double intens = norm(amp)*a1->getrho(mass*mass);
    double err = gh1->GetBinError(i);
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //cout << "chi2_1 = " << chi2 << endl;
  //Intensity2
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double mass = gh2->GetBinCenter(i);
    if(mass<LEFT_INT2 || mass>RIGHT_INT2) continue;
    double s = pow(mass,2);
    cd amp = a2->ExtendA(s,S2,apar2);    
    double intens = norm(amp)*a2->getrho(mass*mass);
    double err = gh2->GetBinError(i);
    chi2 += pow(intens - gh2->GetBinContent(i),2)/(err*err);
  }
  //cout << "chi2_2 = " << chi2 << endl;
  //dPhase
  for(int i=1;i<=gh3->GetNbinsX();i++) {
    double mass = gh3->GetBinCenter(i);
    if(mass<LEFT_PHI || mass>RIGHT_PHI) continue;
    double s = pow(mass,2);
    cd amp1 = a1->ExtendA(s,S2,apar1);
    cd amp2 = a2->ExtendA(s,S2,apar2);
    double phi = arg(amp1*conj(amp2));
    double err = gh3->GetBinError(i);
    double dPhi = normPhase(phi-gh3->GetBinContent(i));
    chi2 += pow(dPhi,2)/(err*err);
  }
  //cout << "chi2_3 = " << chi2 << endl;
  return chi2;
}

void prepareNames(int npar, int apar, string *var_name,const char *pref) {
  for(int i=0;i<npar;i++) { ostringstream rIO; rIO<<pref<<"rN"<<i; var_name[i] = rIO.str();}
  for(int i=0;i<apar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rA"<<i; iIO<<pref<<"iA"<<i; var_name[npar+2*i] = rIO.str(); var_name[npar+2*i+1] = iIO.str();}
}


cd rho3_f2_S(cd s) {
  MIsobar F2Pi(F2_MASS,F2_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  return F2Pi.rho3(s);
}
cd rho3_f2_D(cd s) {
  MIsobar F2Pi(F2_MASS,F2_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  return F2Pi.rho3(s);
}

void fillvectors(const double *par, vector<double> &npar, vector<cd> &apar) {
  //**************** Convert parameters to format **************** //
  //Constract parameters N
  for(int i=0;i<npar.size();i++) npar[i] = par[i];
  //Constract Alpha
  for(int i=0;i<apar.size();i++) { 
    cd normi(par[npar.size()+2*i],par[npar.size()+2*i+1]);  
    apar[i] = normi;
  }
}

//plot
double intensity1(double *x, double *par) {
  vector<cd> apar(APAR_FIT1);
  fillvectors(&par[0],a1->npar,apar);
  cd amp = a1->ExtendA(x[0]*x[0],S2,apar);
  return norm(amp)*a1->getrho(x[0]*x[0]);
}
double intensity2(double *x, double *par) {
  vector<cd> npar(NPAR_FIT2), apar(APAR_FIT2);
  fillvectors(&par[0],a2->npar,apar);
  cd amp = a2->ExtendA(x[0]*x[0],S2,apar);
  return norm(amp)*a2->getrho(x[0]*x[0]);
}

double dphase(double *x, double *par) {
  vector<cd> apar1(APAR_FIT1);
  vector<cd> apar2(APAR_FIT2);
  fillvectors(&par[0],                    a1->npar,apar1);
  fillvectors(&par[NPAR_FIT1+2*APAR_FIT1],a2->npar,apar2);
  cd amp1 = a1->ExtendA(x[0]*x[0],S2,apar1);
  cd amp2 = a2->ExtendA(x[0]*x[0],S2,apar2);
  return arg(amp1*conj(amp2));
}
