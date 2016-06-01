#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
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
#define MAX_INIT_ITERATION 1e3

//******Parameters of the data*********// 3*PI_MASS
#define LEFT_INT1   0.8
#define LEFT_INT2   0.8
#define LEFT_PHI    1.6
#define RIGHT_INT1  2.5
#define RIGHT_INT2  2.3
#define RIGHT_PHI   2.5

using namespace std; 

//functions fit
cd rho(cd s);
void prepareNames(int npar, int ppar, int apar, string *var_name,const char *pref);
TMultiGraph *split_data(const TH1D* h, double low, double up);
template <typename Type> void adjust_plot(Type &obj/*,const char *fulltitle*/);
void setPars(const double *par);
double globChi2(const double *par);

double normPhase(double a);
double phase(cd a1,cd a2) { return normPhase(arg(a1*conj(a2))); };


cd model_wave1(int flag, double s, const double *pars);
cd model_wave2(int flag, double s, const double *pars);

//functions main
int fit_data(int nAttempt, const char *fout_name, int pars_all);

int main(int ac, char **av) {

  if(ac<2) { cerr << "first arg is FIT or PLOT or SHEET" << endl; return 1; }

  const int nAttempt    = (ac > 2) ? atoi(av[2]) : 1;
  const char *file_name = (ac > 3) ? av[3] : "/tmp/test.root";
  const int pars_all    = (ac > 4) ? atoi(av[4]) : 403403;  

  if(strcmp(av[1],"FIT")==0) {
    return fit_data(nAttempt,file_name,pars_all);
  } else if(strcmp(av[1],"PLOT")==0) {
    return 0;/* plot_polar(file_name,nAttempt) */
  } else if(strcmp(av[1],"BAND")==0) {
    return 0;/* get_error_band(nAttempt,file_name); */
  } else { cerr << "first arg is FIT or SHEET or PLOT" << endl; return 1; }
  
  return 0;
}

NoverD *gT1, *gT2;
int gNPAR_FIT1, gPPAR_FIT1, gAPAR_FIT1, gNPAR_FIT2, gPPAR_FIT2, gAPAR_FIT2;

double gNorm1,gNorm2,gNorm3;
double gSth,gMth;
TH1D *gh1, *gh2, *gh12;
TMultiGraph *gm1, *gm2, *gm12;
vector< pair<double,double> > glookup_rho_F2;

MStructureHolder *gStrHolder;
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int build_the_structure() {

  //**************** Load the data *****************//
  TGraph2D g1 ("data/f2piS.txt");  g1.SetName("g1"); 
  TGraph2D g2 ("data/f2piD.txt"); g2.SetName("g2");  
  TGraph2D g12("data/phi_f2piD_f2piS.txt");  g12.SetName("g12");
  const int nBins = 100;
  TH1D *h1 = new TH1D("h1","1^{++}0^{+} f_{0}#pi P;M_{3#pi}",nBins,0.5,2.5);
  TH1D *h2 = new TH1D("h2","1^{++}0^{+} #rho#pi S;M_{3#pi}",nBins,0.5,2.5);
  TH1D *h12 = new TH1D("phi_f0piP_rhopiS","#Delta(f_{0}#pi P, #rho#pi S);M_{3#pi}",nBins,0.5,2.5);
  TH1D *h13 = new TH1D("phi_f0piP_rhopiD","#Delta(f_{0}#pi P, #rho#pi D);M_{3#pi}",nBins,0.5,2.5);
  if(g1.GetN() != h1->GetNbinsX()) {cout << "Very strange!" << endl; return 0;} 
  for(int i=1;i<=nBins;i++) {
    h1->SetBinContent(i,g1.GetY()[i-1]);
    h1->SetBinError  (i,g1.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    h2->SetBinContent(i,g2.GetY()[i-1]);
    h2->SetBinError  (i,g2.GetZ()[i-1]);
  }
  //phase
  for(int i=1;i<=nBins;i++) {
    h12->SetBinContent(i,-normPhase(g12.GetY()[i-1]/180.*M_PI));
    h12->SetBinError  (i,g12.GetZ()[i-1]/180.*M_PI);
  }

  gh1  = h1;  gh1 ->SetStats(kFALSE);
  gh2  = h2;  gh2 ->SetStats(kFALSE);
  gh12 = h12; gh12->SetStats(kFALSE);

  //phase
  //gh12->Scale(M_PI/180.);

  //***************** look up tables ***************//
  gMth = 3*PI_MASS;
  gSth = gMth*gMth;

  const int Nrho = 1000;
  const double Mhg = 3.0;
  for(int i=0;i<Nrho;i++) {
    double mass = gMth + (Mhg - gMth)/(Nrho-1)*i;
    glookup_rho_F2 .push_back( make_pair( mass, real(rho (mass*mass) )) );
  }

  //Normalise intensities
  gNorm1 = gh1 ->GetBinContent(gh1->GetMaximumBin());
  for(int i=1;i<=gh1->GetXaxis()->GetNbins();i++) { gh1->SetBinContent(i,gh1->GetBinContent(i)/gNorm1); gh1->SetBinError(i,gh1->GetBinError(i)/gNorm1); }
  //gh1 ->Scale(1./gNorm1);
  gNorm2 = gh2 ->GetBinContent(gh2->GetMaximumBin());
  for(int i=1;i<=gh2->GetXaxis()->GetNbins();i++) { gh2->SetBinContent(i,gh2->GetBinContent(i)/gNorm2); gh2->SetBinError(i,gh2->GetBinError(i)/gNorm2); }
  //gh2 ->Scale(1./gNorm2);

  //Model
  gStrHolder = new MStructureHolder();
  //intensities
  gStrHolder->AddWave(*gh1,model_wave1,gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1,glookup_rho_F2 ); gStrHolder->SetWaveRange(0, LEFT_INT1,RIGHT_INT1);
  gStrHolder->AddWave(*gh2,model_wave2,gNPAR_FIT2+gPPAR_FIT2+2*gAPAR_FIT2,glookup_rho_F2 ); gStrHolder->SetWaveRange(1, LEFT_INT2,RIGHT_INT2);
  //interference  
  gStrHolder->AddInterference(*gh12,0,1,phase); gStrHolder->SetInterfRange(0, LEFT_PHI, RIGHT_PHI);

  //gStrHolder->JustToPlot(0,0);
  

  //separate data
  gm1  = split_data(gh1,gStrHolder->GetWaveLowRange(0),gStrHolder->GetWaveUpRange(0));
  gm2  = split_data(gh2,gStrHolder->GetWaveLowRange(1),gStrHolder->GetWaveUpRange(1));
  gm12 = split_data(gh12,gStrHolder->GetInterfLowRange(0),gStrHolder->GetInterfUpRange(0));

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////    ////   ///     //////////////////////////////////////////////////////////////////////////////
///// //////// ////// ////////////////////////////////////////////////////////////////////////////////
/////   ////// ////// ////////////////////////////////////////////////////////////////////////////////
///// ///////   ///// ////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int fit_data(int nAttempt, const char *fout_name, int pars_all) {

  cout << "------------ the model with " << pars_all << endl;

  gNPAR_FIT1 = (pars_all % 1000000) / 100000;
  gPPAR_FIT1 = (pars_all % 100000 ) / 10000;
  gAPAR_FIT1 = (pars_all % 10000  ) / 1000;
  gNPAR_FIT2 = (pars_all % 1000   ) / 100;
  gPPAR_FIT2 = (pars_all % 100    ) / 10;
  gAPAR_FIT2 = (pars_all % 10     ) / 1;

  //****************create the model****************//
  /***********/ build_the_structure(); /*************/ 

  const double eTH = 3*PI_MASS, sTH = pow(eTH,2);
  gT1 = new NoverD(gNPAR_FIT1,gPPAR_FIT1,S1,rho,sTH,10000,500);
  gT2 = new NoverD(gNPAR_FIT2,gPPAR_FIT2,S1,rho,sTH,10000,500);

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","plots",0,0,900,600);
  can->Divide(2,2,0.0001,0.0001);

  //**********************Seed**********************//
  gRandom->SetSeed(12);/*12314*/ cout << "The seed " << gRandom->GetSeed() << " is used." << endl;

  //***************** Fit the data *****************//
  const int NumPar1 = gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1;
  const int NumPar2 = gNPAR_FIT2+gPPAR_FIT2+2*gAPAR_FIT2;
  const int NumPar = NumPar1+NumPar2;
  const double step = 0.01;

  // Set the limites for the variables
  string var_name[NumPar];
  prepareNames(gNPAR_FIT1,gPPAR_FIT1,gAPAR_FIT1,&var_name[0],"ch1_");
  prepareNames(gNPAR_FIT2,gPPAR_FIT2,gAPAR_FIT2,&var_name[gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1],"ch2_");

  double  up_limit[NumPar];
  for(int i=0;i<NumPar;i++)  up_limit[i] =  10.;
  for(int i=0;i<gNPAR_FIT1;i++)      up_limit[i] = MRANGE_NPAR/abs(gT1->getmax_value_of_integral(i));
  for(int i=0;i<2*gAPAR_FIT1;i++)    up_limit[i+gNPAR_FIT1+gPPAR_FIT1] = MRANGE_APAR*1.0;
  for(int i=0;i<gNPAR_FIT2;i++)      up_limit[i+           gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1] = MRANGE_NPAR/abs(gT2->getmax_value_of_integral(i));
  for(int i=0;i<2*gAPAR_FIT2;i++)    up_limit[i+gNPAR_FIT2+gPPAR_FIT2+gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1] = MRANGE_APAR*1.0;
  double low_limit[NumPar];
  for(int i=0;i<NumPar;i++) low_limit[i] = -up_limit[i];

  double start_pars[NumPar];
  for(int i=0;i<gNPAR_FIT1;i++)  start_pars[i] = i;
  for(int i=0;i<gPPAR_FIT1;i++)  start_pars[gNPAR_FIT1+i] = 0;
  for(int i=0;i<gAPAR_FIT1;i++) {start_pars[gNPAR_FIT1+gPPAR_FIT1+2*i] = (i==0)? 1.:0.; start_pars[gNPAR_FIT1+gPPAR_FIT1+2*i+1] = 0;}
  for(int i=0;i<gNPAR_FIT2;i++)  start_pars[NumPar1+i] = i;
  for(int i=0;i<gPPAR_FIT2;i++)  start_pars[NumPar1+gNPAR_FIT2+i] = 0;
  for(int i=0;i<gAPAR_FIT2;i++) {start_pars[NumPar1+gNPAR_FIT2+gPPAR_FIT2+2*i] = (i==0)? 1.:0.; start_pars[NumPar1+gNPAR_FIT2+gPPAR_FIT2+2*i+1] = 0;}

  vector<int> unused = {gNPAR_FIT1+gPPAR_FIT1+1/*,
			//gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1,
			gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1+1,
			gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1+2*/};

//  const int NstepsInfit = 1;
//  vector<int> fix[NstepsInfit] = { {}};
//  vector<int> int_exclude[NstepsInfit] = { {},
//					   {0},
//					   {1,2}};
//  vector<int> phi_exclude[NstepsInfit] = { {0,2},
//					   {0},
//					   {2}};

  //Do minimisation several times
  TFile fout(fout_name,"RECREATE");
  TTree t("mins","several_mimimas");
  double final_pars[NumPar];
  for(int i=0;i<NumPar;i++) t.Branch(var_name[i].c_str(),&final_pars[i]);
  t.Branch("canva","TCanvas",&can);
  double chi2; t.Branch("chi2",   &chi2);
  int status; t.Branch("status",&status);
  TGraphErrors *gr[2][3], *gr_whole[2][3], *gphi[1], *gphi_whole[1];
  for(int w=0;w<2;w++) for(int i=0;i<3;i++) {gr[w][i] = 0; gr_whole[w][i] = 0; }
  for(int w=0;w<1;w++) {gphi[w] = 0; gphi_whole[w] = 0; }
  for(int e=0;e<nAttempt;e++) {
    cout << "---------- Attempt " << e << " -----------"<< endl;
    cout << "------------------------------------------"<< endl;

    /**************************************MINIMIZE******************************************/
    /****************************************************************************************/
    //Build minimizer
    ROOT::Math::Minimizer* min = 
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    
    // set tolerance , etc...
    min->SetMaxFunctionCalls(100000);
    min->SetTolerance(0.001);
    min->SetStrategy(1);
    min->SetPrintLevel(1);
    min->Options().Print();
    
    // Create funciton wrapper for minmizer a IMultiGenFunction type 
    ROOT::Math::Functor functor(&globChi2,NumPar);
    min->SetFunction(functor);

    for(int i=0;i<NumPar;i++) min->SetVariable( i,var_name[i].c_str(), 0.0, step);//,low_limit[i],up_limit[i] start_pars[i]
    
    //step
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,up_limit[i]/1e3);
    //starting point
    for(int i=0;i<NumPar;i++)  
      start_pars[i] = low_limit[i] + (up_limit[i]-low_limit[i])*gRandom->Rndm();
    for(int i=0;i<unused.size();i++) {
      start_pars[unused[i]] = 0;
      min->FixVariable(unused[i]);
    }
    min->SetVariableValues(start_pars); globChi2(start_pars);
    cout<<"Initual parameters:\n";for(int i=0;i<NumPar;i++) {cout<<var_name[i]<<" = "<<start_pars[i]<<";\n";}

    //do the set of minimizations with stategy
//    for(int j=0;j<NstepsInfit;j++) {//NstepsInfit
//      for(int i=0;i<fix[j].size();i++) min->FixVariable    (fix[j][i]);
//      for(int i=0;i<int_exclude[j].size();i++) gStrHolder->JustToPlot(int_exclude[j][i],true);
//      for(int i=0;i<phi_exclude[j].size();i++) gStrHolder->JustToPlot(0,phi_exclude[j][i],true);
//      gStrHolder->Print();
//      min->Minimize();
//      for(int i=0;i<fix[j].size();i++) min->ReleaseVariable(fix[j][i]);
//      for(int i=0;i<int_exclude[j].size();i++) gStrHolder->JustToPlot(int_exclude[j][i],false);
//      for(int i=0;i<phi_exclude[j].size();i++) gStrHolder->JustToPlot(0,phi_exclude[j][i],false);
//    }
    
    min->Minimize();
    //finish
    bool fit_stat = 1;
    status = min->Status(); //if(status > 1) continue;
    memcpy(final_pars,min->X(),NumPar*sizeof(double));
    chi2 = globChi2(final_pars);
    cout<<"Final parameters:\n";for(int i=0;i<NumPar;i++) {cout<<var_name[i]<<" = "<<final_pars[i]<<";\n";}

    /*****************************************PLOT*******************************************/
    /****************************************************************************************/
    //calculate functions to plot
    double arr_pars[NumPar];
    memcpy(arr_pars,min->X(),NumPar*sizeof(double));
    setPars(arr_pars);
    //intesities
    for(int w=0;w<2;w++) 
      for(int i=0;i<3;i++) {
	gr[w][i] = gStrHolder->GetWavePlot  (w,i, 0,0, gr[w][i] );
	if(i==0) gr[w][i]->SetLineColor(kRed);
	if(i==1) gr[w][i]->SetLineColor(kBlue);
	if(i==2) gr[w][i]->SetLineColor(kGreen);
      }
    for(int w=0;w<2;w++) 
      for(int i=0;i<3;i++) {
	gr_whole[w][i] = gStrHolder->GetWavePlot  (w,i, 0.5,2.5, gr_whole[w][i] );
	gr_whole[w][i]->SetLineStyle(2);
	gr_whole[w][i]->SetLineWidth(1);
	if(i==0) gr_whole[w][i]->SetLineColor(kRed);
	if(i==1) gr_whole[w][i]->SetLineColor(kBlue);
	if(i==2) gr_whole[w][i]->SetLineColor(kGreen);
      }
    //interferences
    gphi[0] = gStrHolder->GetInterfPlot(0, 0,0, gphi[0]); gphi[0]->SetLineColor(kRed);
    gphi_whole[0] = gStrHolder->GetInterfPlot(0, 0.5,2.5, gphi_whole[0]); gphi_whole[0]->SetLineStyle(2); gphi_whole[0]->SetLineColor(kRed);
    //plot
    //intensities
    TMultiGraph m[3];
    m[0].Add((TMultiGraph*)gm1->Clone()); m[1].Add((TMultiGraph*)gm2->Clone());
    m[0].SetTitle(gm1->GetTitle()); m[1].SetTitle(gm2->GetTitle());
    for(int w=0;w<2;w++) {
      for(int i=2;i>=0;i--) {//for(int i=2;i>=0;i--) {
	m[w].Add((TGraphErrors*)gr[w][i]->Clone(),"lz");
	m[w].Add((TGraphErrors*)gr_whole[w][i]->Clone(),"lz");	
      }
      int Npad = 1+w+w*2;
      can->cd(Npad); m[w].Draw("a");
      adjust_plot(*m[w].GetHistogram());
    }
    //interferences
    m[2].Add((TMultiGraph*)gm12->Clone());
    m[2].SetTitle(gm12->GetTitle()); 
    {
      m[2].Add((TGraphErrors*)gphi[0]->Clone(),"lz"); 
      m[2].Add((TGraphErrors*)gphi_whole[0]->Clone(),"lz");
      int Npad = 2; can->cd(Npad); m[2].Draw("a");
      adjust_plot(*m[2].GetHistogram());
    }
    //save
    if(nAttempt==1) can->SaveAs("c1.png");
    t.Fill();
    delete min;
  }

  t.Write();
  fout.Close();
  
  cout << "finished" << endl;

  return 0;
}

double globChi2(const double *par) {  
  setPars(par);
  double chi0 = gStrHolder->GetWaveChi2(0);
  double chi1 = gStrHolder->GetWaveChi2(1);
  double interfChi2 = gStrHolder->GetInterferenceChi2(0);
  double chi2 = chi0+chi1+interfChi2;//gStrHolder->GetChi2();
  cout << "chi2 = " << setw(16) << setprecision(15) << chi0;
  cout << ", second chi2 = " << setw(16) << setprecision(15) << chi1;
  cout << ", last chiPhi = " << setw(16) << setprecision(15) << interfChi2;
  cout << ", final chi2 = " << setw(16) << setprecision(15) << chi2 << endl;
  return chi2;
}

void setPars(const double *par) {
  gStrHolder->pars[0] = vector<double>(par,par+gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1);
  //cout << "P1:"; for(int i=0;i<gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1;i++) cout<<gStrHolder->pars[0][i]<<" "; cout << endl;
  gStrHolder->pars[1] = vector<double>(par+gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1,
				       par+gNPAR_FIT1+gPPAR_FIT1+2*gAPAR_FIT1
				          +gNPAR_FIT2+gPPAR_FIT2+2*gAPAR_FIT2);
  //cout << "P2:"; for(int i=0;i<gNPAR_FIT2+gPPAR_FIT2+2*gAPAR_FIT2;i++) cout<<gStrHolder->pars[1][i]<<" "; cout << endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
/* f2 pi S-wave */
cd model_wave1(int flag, double s, const double *pars) {

  for(int i=0;i<gNPAR_FIT1;i++) gT1->npar[i] = pars[i];
  for(int i=0;i<gPPAR_FIT1;i++) gT1->ppar[i] = pars[gNPAR_FIT1+i];
  cd T = gT1->A(s);

  vector<cd> apar(gAPAR_FIT1);
  for(int i=0;i<gAPAR_FIT1;i++) apar[i] = cd(pars[gNPAR_FIT1+gPPAR_FIT1+2*i],
					     pars[gNPAR_FIT1+gPPAR_FIT1+2*i+1]);
  
  cd Aprod = NoverD::cexpand(s,S2,apar);
  //cout << Aprod << " " << T << " " << Aprod*T << endl;
  //cout << Aprod*T << endl;
  if(flag==0) return Aprod*T;
  if(flag==1) return apar[0]*T;
  if(flag==2) return (Aprod-apar[0])*T;

  cerr << "Chech the code. Something is going wrong" << endl;  return 0;  
}

////////////////////////////////////////////////////////////////////////////////////////
/* rho pi D-wave */
cd model_wave2(int flag, double s, const double *pars) {  

  //for(int i=0;i<gNPAR_FIT2;i++) cout << pars[i] << " "; cout << endl;
  
  for(int i=0;i<gNPAR_FIT2;i++) gT2->npar[i] = pars[i];
  for(int i=0;i<gPPAR_FIT2;i++) gT2->ppar[i] = pars[gNPAR_FIT2+i];
  cd T = gT2->A(s);

  vector<cd> apar(gAPAR_FIT2);
  for(int i=0;i<gAPAR_FIT2;i++) apar[i] = cd(pars[gNPAR_FIT2+gPPAR_FIT2+2*i],
					     pars[gNPAR_FIT2+gPPAR_FIT2+2*i+1]);

  cd Aprod = NoverD::cexpand(s,S2,apar);
  //cout << Aprod << " " << T << " " << Aprod*T << endl;
  if(flag==0) return Aprod*T;
  if(flag==1) return apar[0]*T;
  if(flag==2) return (Aprod-apar[0])*T;
  
  cerr << "Chech the code. Something is going wrong" << endl;  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

cd rho(cd s) {
  if(real(s)>1e3 && fabs(imag(s))<1e-5) return 1./(8*M_PI)*(1.0-(pow(F2_MASS,2)+pow(PI_MASS,2))/s);
  MIsobar RhoPi(F2_MASS,F2_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  cd value = RhoPi.rho3(s);
  return value;
}


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

TMultiGraph *split_data(const TH1D* h, double low, double up) {
  TGraphErrors *gray = new TGraphErrors(0); gray->SetLineColor(kGray); gray->SetMarkerColor(kGray);
  TGraphErrors *dark = new TGraphErrors(0);

  int Nbins = h->GetXaxis()->GetNbins();
  double dm = h->GetBinCenter(2) - h->GetBinCenter(1);
  for(int i=0;i<Nbins;i++) {
    double mass = h->GetBinCenter(i+1);
    if(mass<low || mass >up) {
      if(i!=0&&h->GetBinContent(i+1)==0&&h->GetBinError(i+1)==0) continue;
      gray->Set(gray->GetN()+1);
      gray->SetPoint     (gray->GetN()-1,mass,h->GetBinContent(i+1));
      gray->SetPointError(gray->GetN()-1,dm  ,h->GetBinError  (i+1));
    } else {
      dark->Set(dark->GetN()+1);
      dark->SetPoint     (dark->GetN()-1,mass,h->GetBinContent(i+1));
      dark->SetPointError(dark->GetN()-1,dm  ,h->GetBinError  (i+1));
    }
  }
  TMultiGraph *m = new TMultiGraph();
  m->Add(gray,"pz");
  m->Add(dark,"pz");
  m->SetTitle(TString::Format("%s;%s;%s",h->GetTitle(),h->GetXaxis()->GetTitle(),h->GetYaxis()->GetTitle()));  
  return m;
}

template <typename Type> void adjust_plot(Type &obj/*,const char *fulltitle*/) {
  obj.SetLabelSize(0.06,"xy");
  obj.GetXaxis()->SetTitleSize(0.07);
  obj.GetXaxis()->SetTitleOffset(-0.5);
  //obj.SetTitle(fulltitle);
}


void prepareNames(int npar, int ppar, int apar, string *var_name,const char *pref) {
  for(int i=0;i<npar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rN"<<i; var_name[i] = rIO.str();}
  for(int i=0;i<ppar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rP"<<i; var_name[npar+i] = rIO.str();}
  for(int i=0;i<apar;i++) { ostringstream rIO,iIO; rIO<<pref<<"rA"<<i; iIO<<pref<<"iA"<<i; var_name[npar+ppar+2*i] = rIO.str(); var_name[npar+ppar+2*i+1] = iIO.str();}
}


double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);
  if(a<-M_PI) return normPhase(a+2*M_PI);
  return a;
}
