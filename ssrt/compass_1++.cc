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
#include <MCoupledChannelIsobar.h>
#include <NoverD.h>

#include <constants.h>
#include <deflib.h>

//******Parameters of the data*********// 3*PI_MASS
#define LEFT_INT1   3*PI_MASS
#define LEFT_INT2   1.20
#define LEFT_PHI    1.20
#define RIGHT_INT1  2.1
#define RIGHT_INT2  2.3
#define RIGHT_PHI   1.95

using namespace std; 

//functions data
double dphase(double *x, double *par);
double intensity1(double *x, double *par);
double intensity2(double *x, double *par);
//
double intensityDeck1(double *x, double *par);
double intensityDeck2(double *x, double *par);
double intensityTr1(double *x, double *par);
double intensityBW(double *x, double *par);

//functions fit
double globChi2(const double *par);
double normPhase(double a);
//addition function
cd awave1(double s, double m1, double g1, cd c1, double b1, double sl1, cd d1_c);
cd awave2(double s, cd tr1_c, cd tr2_c, double b2, double sl2, cd d2_c,
	  /*wave1*/ double m1, double g1, cd c1, double b1, double sl1, cd d1_c);
cd a1_bw(double s, double m1, double g);
double deck1(double s, double b, double slope, vector<pair<double,double> > &table);
double ndeck1(double s, double b, double slope, vector<pair<double,double> > &table);
template <typename Type> Type getvalue(double M, vector< pair<double,Type> > &table);
//double getvalue(double M, vector<double> &table);
double rhoRHO(double s);
double rhoF0(double s);

//functions main
int fit_data(int nAttempt, const char *fout_name, int pars_all);
int plot_best_result(const char *fin_name) {return 0;};
int plot_best_sheets(const char *fin_name) {return 0;};

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

double gSth,gMth;
TH1D *gh1, *gh2, *gh3;
int gNPAR_FIT1, gAPAR_FIT1, gNPAR_FIT2, gAPAR_FIT2;
vector< pair<double,double> > glookup_rho_RHO, glookup_rho_F0;
vector< pair<double,cd> > glookup_rhopiTof0pi, glookup_KstarKTof0pi;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int fit_data(int nAttempt, const char *fout_name, int pars_all) {

  MIsobar test(RHO_MASS,RHO_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  double rho3 = test.rho3(1.27*1.27);
  cout << " g ~= " << sqrt(0.5*2*1.27/rho3) << endl;

  //**********************Seed**********************//
  gRandom->SetSeed(0);/*12314*/ cout << "The seed " << gRandom->GetSeed() << " is used." << endl;

  //**************** Load the data *****************//
  TGraph2D rhopiS("data/rhopiS.txt"); rhopiS.SetName("grhopiS");  
  TGraph2D f0piP("data/f0piP.txt"); f0piP.SetName("gf0piP"); 
  TGraph2D phi_f0piP_rhopiS("data/phi_f0piP_rhopiS.txt"); phi_f0piP_rhopiS.SetName("gphi_f0piP_rhopiS"); 
  const int nBins = 100;
  TH1D *rhopiS_int = new TH1D("rhopiS","#rho#pi S;M_{3#pi}",nBins,0.5,2.5);
  TH1D *f0piP_int = new TH1D("f0piP","f_{0}#pi D;M_{3#pi}",nBins,0.5,2.5);
  TH1D *onepp_dphi = new TH1D("phi_f0piP_rhopiS","#Delta(#rho#pi S and f_{0}#pi P);M_{3#pi}",nBins,0.5,2.5);
  if(rhopiS.GetN() != rhopiS_int->GetNbinsX()) {cout << "Very strange!" << endl; return 0;} 
  for(int i=1;i<=nBins;i++) {
    rhopiS_int->SetBinContent(i,rhopiS.GetY()[i-1]);
    rhopiS_int->SetBinError  (i,rhopiS.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    f0piP_int->SetBinContent(i,f0piP.GetY()[i-1]);
    f0piP_int->SetBinError  (i,f0piP.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    onepp_dphi->SetBinContent(i,normPhase(phi_f0piP_rhopiS.GetY()[i-1]/180.*M_PI));
    onepp_dphi->SetBinError  (i,phi_f0piP_rhopiS.GetZ()[i-1]/180.*M_PI);
  }

  gh1 = rhopiS_int;
  gh2 = f0piP_int;
  gh3 = onepp_dphi; gh3->GetYaxis()->SetRangeUser(-M_PI,M_PI); 

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

  //***************** look up tables ***************//
  gMth = 3*PI_MASS;
  gSth = gMth*gMth;

  const int Nrho = 1000;
  const double Mhg = 3.0;
  for(int i=0;i<Nrho;i++) {
    double mass = gMth + (Mhg - gMth)/(Nrho-1)*i;
    glookup_rho_RHO.push_back( make_pair( mass, rhoRHO(mass*mass) ) );
    glookup_rho_F0 .push_back( make_pair( mass, rhoF0 (mass*mass) ) );
  }

  TGraph2D gr_tr1("calculations/sclr_trgl_KstarKTof0pi.out"); gr_tr1.SetName("tr1");//_two_isospin_channels
  for(int i=0;i<gr_tr1.GetN();i++)
    glookup_rhopiTof0pi.push_back( make_pair( gr_tr1.GetX()[i], 
					      cd(gr_tr1.GetY()[i],gr_tr1.GetZ()[i]) ) );

  TGraph2D gr_tr2("calculations/sclr_trgl_rhopiTof0pi.out"); gr_tr2.SetName("tr2");
  for(int i=0;i<gr_tr2.GetN();i++)
    glookup_KstarKTof0pi.push_back( make_pair( gr_tr2.GetX()[i], 
					       cd(gr_tr2.GetY()[i],gr_tr2.GetZ()[i]) ) );
  
  //normalizetion
  max = abs( glookup_rhopiTof0pi[0].second );
  for(int i=1;i<glookup_rhopiTof0pi .size();i++) if(abs(glookup_rhopiTof0pi [i].second)>max) max = abs(glookup_rhopiTof0pi [i].second); 
  for(int i=0;i<glookup_rhopiTof0pi .size();i++) glookup_rhopiTof0pi[i].second *= 1./max;
  max = abs( glookup_KstarKTof0pi[0].second );
  for(int i=1;i<glookup_KstarKTof0pi.size();i++) if(abs(glookup_KstarKTof0pi[i].second)>max) max = abs(glookup_KstarKTof0pi[i].second); 
  for(int i=0;i<glookup_KstarKTof0pi.size();i++) glookup_KstarKTof0pi[i].second *= 1./max;  
  
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
  
  const int NumPar1 = 8;
  const int NumPar2 = 8;
  const int NumPar = NumPar1+NumPar2;

  // Create funciton wrapper for minmizer a IMultiGenFunction type 
  ROOT::Math::Functor functor(&globChi2,NumPar);
  min->SetFunction(functor);
    
  const double step = 0.01;
  
  // Set the limites for the variables
  string var_name[NumPar] = {"a1_mass","a1_g","a1_rc0","a1_ic0",
			     "b1","slope1","d1_rc0","d1_ic0",
			     "tr1_rc","tr1_ic","tr2_rc","tr2_ic",
			     "b2","slope2","d2_rc0","d2_ic0"};
//  double start_pars[NumPar] = {1.2,  6., 1.0, 0.0, 
//			       1.0, -2.0, 0.0, 0.0,
//			       0.5,  0.0, 0.0, 0.0,
//			       1.0, -2.0, 0.0, 0.0};
  double start_pars[NumPar] = {1.2,  7, 0.07, 0.0, 
			       1.0, -2.0, 0.0, 0.0,
			       1.0,  0.0, 0.0, 0.0,
			       1.0, -2.0, 0.0, 0.0};
  double  up_limit[NumPar] = {1.3, 7.0, 0.08, 0.0,
			      10.0, 0.0, .1, .1,
			      2.0, 1.0, 1.0, 1.0,
			      5.0, 0.0, .1, .1};
  double low_limit[NumPar] = {1.1, 3,  0.02, -0.0,
			      0.0, -10.0, -.1, -.1,
			      1.0, -1.0, -1.0, -1.0,
			      0.0, -10.0, -.1, -.1};
  //WAVE 1
  for(int i=0;i<NumPar;i++) min->SetVariable( i,var_name[i].c_str(), start_pars[i], step );
  for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,(up_limit[i]-low_limit[i])/1e3);
 
  vector<int> unused = {3,10,11};
  const int NstepsInfit = 3;
  vector<int> fix[NstepsInfit] = { {0,1,4,5,6,7,12,13,14,15},
				   {0,1,12,13,14,15},
				   {0,1,4,5,6,7}};
    
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
    //step
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,(up_limit[i]-low_limit[i])/1e3);
    // starting point
    for(int i=0;i<NumPar;i++)  
      start_pars[i] = low_limit[i] + (up_limit[i]-low_limit[i])*gRandom->Rndm();
    min->SetVariableValues(start_pars);
    for(int i=0;i<unused.size();i++) {
      min->SetVariableValue(unused[i],0);
      min->FixVariable(unused[i]);
    }

    //do the minimization
    for(int j=0;j<NstepsInfit;j++) {//
      //step J
      for(int i=0;i<fix[j].size();i++) min->FixVariable    (fix[j][i]);
      min->Minimize();
      for(int i=0;i<fix[j].size();i++) min->ReleaseVariable(fix[j][i]);
    }

    min->Minimize();
    //finish
    bool fit_stat = 1;
    status = min->Status(); //if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);

    TF1 *fInt1 = new TF1("int1",intensity1,LEFT_INT1,RIGHT_INT1,NumPar1); fInt1->SetParameters(final_pars);
    TF1 *fInt2 = new TF1("int2",intensity2,LEFT_INT2,RIGHT_INT2,NumPar ); fInt2->SetParameters(final_pars);
    TF1 *fPhi  = new TF1("phi" ,dphase,    LEFT_PHI ,RIGHT_PHI ,NumPar ); fPhi ->SetParameters(final_pars);

    TF1 *fInt1_out = new TF1("int1",intensity1,3*PI_MASS,2.5,NumPar1); fInt1_out->SetLineColor(kGray); fInt1_out->SetLineWidth(0.5); fInt1_out->SetParameters(final_pars);
    TF1 *fInt2_out = new TF1("int2",intensity2,3*PI_MASS,2.5,NumPar ); fInt2_out->SetLineColor(kGray); fInt2_out->SetLineWidth(0.5); fInt2_out->SetParameters(final_pars);
    TF1 * fPhi_out = new TF1("phi" ,dphase,    3*PI_MASS,2.5,NumPar );  fPhi_out->SetLineColor(kGray);  fPhi_out->SetLineWidth(0.5);  fPhi_out->SetParameters(final_pars);

    TF1 *fInt1_bw    = new TF1("int1_bw"  ,intensityBW   ,3*PI_MASS,2.5,4);   fInt1_bw  ->SetLineColor(kBlue);  fInt1_bw  ->SetLineWidth(1); fInt1_bw  ->SetParameters(final_pars);
    TF1 *fInt1_deck  = new TF1("int1_deck",intensityDeck1,3*PI_MASS,2.5,4);   fInt1_deck->SetLineColor(kGreen); fInt1_deck->SetLineWidth(1); fInt1_deck->SetParameters(&final_pars[4]);
    TF1 *fInt2_tr    = new TF1("int2_tr",  intensityTr1,  3*PI_MASS,2.5,10 ); fInt2_tr  ->SetLineColor(kBlue);  fInt2_tr  ->SetLineWidth(1); fInt2_tr->SetParameters(final_pars);
    TF1 *fInt2_deck  = new TF1("int2_deck",intensityDeck2,3*PI_MASS,2.5,4 );  fInt2_deck->SetLineColor(kGreen); fInt2_deck->SetLineWidth(1); fInt2_deck->SetParameters(&final_pars[12]);

    can->cd(1); gh1->Draw(); fInt1_out->Draw("same"); fInt1->Draw("same"); fInt1_bw->Draw("same"); fInt1_deck->Draw("same");
    can->cd(4); gh2->Draw(); fInt2_out->Draw("same"); fInt2->Draw("same"); fInt2_tr->Draw("same"); fInt2_deck->Draw("same");
    can->cd(2); gh3->Draw();  fPhi_out->Draw("same");  fPhi->Draw("same"); 
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
  //WAVE 1
  double m1  = par[0];
  double g1  = par[1];
  cd c1(par[2],par[3]);
  double b1  = par[4];
  double sl1 = par[5];
  cd d1_c(par[6],par[7]);
  //WAVE 2
  cd tr1_c(par[8],par[9]);
  cd tr2_c(par[10],par[11]);
  double b2     = par[12];
  double sl2 = par[13];
  cd d2_c(par[14],par[15]);

  // *********************** Calculate chi2 ********************** //
  double chi2 = 0;
  //Intensity1
  for(int i=1;i<=gh1->GetNbinsX();i++) {
    double mass = gh1->GetBinCenter(i);
    if(mass<LEFT_INT1 || mass>RIGHT_INT1) continue;
    double s = mass*mass;
    cd amp = awave1(mass*mass,m1,g1,c1,b1,sl1,d1_c);
    double intens = norm(amp)*getvalue(sqrt(s), glookup_rho_RHO);
    double err = gh1->GetBinError(i);
    chi2 += pow(intens - gh1->GetBinContent(i),2)/(err*err);
  }
  //cout << "chi2 = " << chi2;
  //Intensity2
  for(int i=1;i<=gh2->GetNbinsX();i++) {
    double mass = gh2->GetBinCenter(i);
    if(mass<LEFT_INT2 || mass>RIGHT_INT2) continue;
    double s = mass*mass;
    cd amp = awave2(s,tr1_c,tr2_c,b2,sl2,d2_c,/*wave1*/m1,g1,c1,b1,sl1,d1_c);
    double intens = norm(amp)*getvalue(sqrt(s), glookup_rho_F0);
    double err = gh2->GetBinError(i);
    chi2 += pow(intens - gh2->GetBinContent(i),2)/(err*err);
  }
  //cout << ", second chi2 = " << chi2;
  //Phase
  for(int i=1;i<=gh3->GetNbinsX();i++) {
    double mass = gh3->GetBinCenter(i);
    if(mass<LEFT_PHI || mass>RIGHT_PHI) continue;
    double s = pow(mass,2);
    cd amp1 = awave1(s,m1,g1,c1,b1,sl1,d1_c);
    cd amp2 = awave2(s,tr1_c,tr2_c,b2,sl2,d2_c,/*wave1*/m1,g1,c1,b1,sl1,d1_c);
    double phase = normPhase( arg(amp2*conj(amp1)) );
    double err = gh3->GetBinError(i);
    double dPhi = normPhase(phase-gh3->GetBinContent(i));
    chi2 += pow(dPhi,2)/(err*err);
  }
  //cout << ", final chi2 = " << chi2<<endl;
  return chi2;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

cd awave1(double s, double m1, double g1, cd c1, double b1, double sl1, cd d1_c) {
  cd bw = a1_bw(s,m1,g1);
  double d1 = deck1(s,b1,sl1,glookup_rho_RHO);
  return c1*bw + d1_c*d1;
}

cd awave2(double s, cd tr1_c, cd tr2_c, double b2, double sl2, cd d2_c,
	  /*wave1*/ double m1, double g1, cd c1, double b1, double sl1, cd d1_c) {
  cd bw = a1_bw(s,m1,g1);  
  cd tr1 = c1*bw * getvalue(sqrt(s),glookup_rhopiTof0pi);
  double d2  = deck1(s,b2,sl2,glookup_rho_F0);
  cd value = tr1_c*tr1 + d2_c*d2;
  if(tr2_c == 0.0) return value;
  cd aw1 = awave1(s,m1,g1,c1,b1,sl1,d1_c);
  cd tr2 = aw1 * getvalue(sqrt(s),glookup_KstarKTof0pi);
  value += tr2_c*tr2;
  return value;
}

cd a1_bw(double s, double m1, double g) {
  double rho = getvalue(sqrt(s), glookup_rho_RHO);
  cd value = g*g/(m1*m1-s- 0.5*g*g*rho*cd(0,1));
  return value;
}

double deck1(double s, double b, double slope, vector<pair<double,double> > &table) {
  double rho = getvalue(sqrt(s), table);//1./(8*M_PI)*sqrt((s-pow(gMth,2))/s);
  double p = sqrt(s)/2.0*(8*M_PI)*rho;
  double value = pow(sqrt(s)-gMth,b)*exp(slope*p);
  //norm
  return value;
}

double ndeck1(double s, double b, double slope, vector<pair<double,double> > &table) {
  const double Nloop = 100000;
  double sum = 0;
  double hM = 2.5;
  for(int i=0;i<Nloop;i++) {
    double si = gMth + (hM-gMth)*gRandom->Rndm();
    sum += deck1(si,b,slope,table);
  }
  sum *= (hM-gMth)/Nloop;
  //norm
  double value = deck1(s,b,slope,table)/sum;
  return value;//sqrt(8*M_PI)*
}

template <typename Type>
Type getvalue(double M, vector< pair<double,Type> > &table) {
  const int N = table.size();
  const double lft = table[0].first;
  const double rht = table[N-1].first;
  const double Mstep = table[1].first - lft;
  const int Nsteps = (M - lft)/Mstep;
  if(Nsteps<0 || Nsteps>=N-1) {cerr<<"Error!! in getvalue!"<<endl; return 1;}
  const Type value = table[Nsteps].second + 
    ( table[Nsteps+1].second - table[Nsteps].second ) / 
    ( table[Nsteps+1].first  - table[Nsteps].first  ) * (M - table[Nsteps].first); 
  return value;
}

double rhoRHO(double s) {
  //if(s>1e3) return 1./(8*M_PI)*(1.0-(pow(RHO_MASS,2)+pow(PI_MASS,2))/s);
  MIsobar RhoPi(RHO_MASS,RHO_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  double value = RhoPi.rho3(s);
  //cout << "double rhoRHO("<<s<<") " << endl;
  return value;
}

double rhoF0(double s) {
  //if(s>1e3) return 1./(8*M_PI)*(1.0-(pow(F0_MASS,2)+pow(PI_MASS,2))/s);
  MCoupledChannelIsobar F0Pi(F0_MASS,F0_COUPLING_PI,F0_COUPLING_K,
			     PI_MASS,PI_MASS,PI_MASS,K_MASS,K_MASS);
  double value = F0Pi.rho3(s);
  //cout << "double rhoF0("<<s<<") " << endl;
  return value;
}


////////////////////////////////////////////////////////////////////////
//                     /////////////////////////////////////////////////
////////////////////                   ///////////////            //////
/////////////////////////////////////                    ///////////////
////////////////////////////////////////////////////////////////////////

double intensity1(double *x, double *par) {  
  //WAVE 1
  double m1  = par[0];
  double g1  = par[1];
  cd c1(par[2],par[3]);
  double b1  = par[4];
  double sl1 = par[5];
  cd d1_c(par[6],par[7]);
  //amp
  cd amp = awave1(x[0]*x[0],m1,g1,c1,b1,sl1,d1_c);
  return norm(amp)*getvalue(x[0],glookup_rho_RHO);
}

double intensity2(double *x, double *par) {
  //WAVE 1
  double m1  = par[0];
  double g1  = par[1];
  cd c1(par[2],par[3]);
  double b1  = par[4];
  double sl1 = par[5];
  cd d1_c(par[6],par[7]);
  //WAVE 2
  cd tr1_c(par[8],par[9]);
  cd tr2_c(par[10],par[11]);
  double b2     = par[12];
  double sl2 = par[13];
  cd d2_c(par[14],par[15]);
  //amp
  cd amp = awave2(x[0]*x[0],tr1_c,tr2_c,b2,sl2,d2_c,/*wave1*/m1,g1,c1,b1,sl1,d1_c);
//  cout << getvalue(sqrt(s), glookup_rho_F0) 
//       << ", intens = " << intens
//       << ", norm = " << norm(amp)
//       << endl;
  return norm(amp)*getvalue(x[0],glookup_rho_F0);
}

double dphase(double *x, double *par) {
  //WAVE 1
  double m1  = par[0];
  double g1  = par[1];
  cd c1(par[2],par[3]);
  double b1  = par[4];
  double sl1 = par[5];
  cd d1_c(par[6],par[7]);
  //WAVE 2
  cd tr1_c(par[8],par[9]);
  cd tr2_c(par[10],par[11]);
  double b2     = par[12];
  double sl2 = par[13];
  cd d2_c(par[14],par[15]);
  //amps
  cd amp1 = awave1(x[0]*x[0],m1,g1,c1,b1,sl1,d1_c);
  cd amp2 = awave2(x[0]*x[0],tr1_c,tr2_c,b2,sl2,d2_c,/*wave1*/m1,g1,c1,b1,sl1,d1_c);
  return normPhase( arg(amp2*conj(amp1)) );
}

/////////////////////////////////////////////////////////////////////////////////////////

double intensityBW(double *x, double *par) {  
  //BW 1
  double m1  = par[0];
  double g1  = par[1];
  cd c1(par[2],par[3]);
  //amp
  cd amp = c1*a1_bw(x[0]*x[0],m1,g1);
  return norm(amp)*getvalue(x[0],glookup_rho_RHO);
}

double intensityDeck1(double *x, double *par) {  
  //DECK1 1
  double b1  = par[0];
  double sl1 = par[1];
  cd d1_c(par[2],par[3]);
  //amp
  cd amp = d1_c*deck1(x[0]*x[0],b1,sl1,glookup_rho_RHO);
  return norm(amp)*getvalue(x[0],glookup_rho_RHO);
}

double intensityTr1(double *x, double *par) {  
  //WAVE 1
  double m1  = par[0];
  double g1  = par[1];
  cd c1(par[2],par[3]);
  double b1  = par[4];
  double sl1 = par[5];
  cd d1_c(par[6],par[7]);
  //WAVE 2
  cd tr1_c(par[8],par[9]);
  //TRIANGLE 1
  cd aw1 = awave1(x[0]*x[0],m1,g1,c1,b1,sl1,d1_c);
  cd tr1 = aw1 * getvalue(x[0],glookup_rhopiTof0pi);  
  //amp
  cd amp = tr1_c*tr1;
  return norm(amp)*getvalue(x[0],glookup_rho_F0);
}

double intensityDeck2(double *x, double *par) {  
  //DECK 2
  double b2     = par[0];
  double sl2 = par[1];
  cd d2_c(par[2],par[3]);
  //amp
  cd amp = d2_c*deck1(x[0]*x[0],b2,sl2,glookup_rho_F0);
  return norm(amp)*getvalue(x[0],glookup_rho_F0);
}



/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);
  if(a<-M_PI) return normPhase(a+2*M_PI);
  return a;
}
