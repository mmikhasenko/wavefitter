// Copyright [09.2015 September 2015]
// Author: Misha Mikhasenko

#include <constants.h>
#include <deflib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <initializer_list>

#include "TGraph2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TRandom.h"
#include "TVirtualFitter.h"
#include "TMultiGraph.h"
#include "TArrow.h"
#include "TLine.h"

#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include <MIsobar.h>
#include <MCoupledChannelIsobar.h>
#include <MStructureHolder.h>

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
// #include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <MatrixInverse.h>

using namespace std;

// functions fit
void setPars(const double *par);
void scalePars(double *start_pars);
double globChi2(const double *par);
double normPhase(double a);
double phase(cd a1, cd a2) { return normPhase(arg(a1*conj(a2))); }
template <typename Type> Type getvalue(double M, vector< pair<double, Type> > & table);
// rho
double rhoRHO(double s);
double rhoF2(double s);

cd model_wave1(int flag, double s, const double *pars);
cd model_wave2(int flag, double s, const double *pars);
cd model_wave3(int flag, double s, const double *pars);
boost::numeric::ublas::matrix<cd> calculate_matrix(cd s, const double *par);
boost::numeric::ublas::vector<cd> calculate_vector(cd s, const double *par);

template<class T> bool InvertMatrix(const boost::numeric::ublas::matrix<T>& input, 
				      boost::numeric::ublas::matrix<T>& inverse);

template<class T> boost::numeric::ublas::matrix<T> gjinverse(const boost::numeric::ublas::matrix<T> &m, 
							       bool &singular);

// functions main
int fit_data(int nAttempt, const char *fout_name, int pars_all);
int plot_best_result(const char *fin_name) {return 0;}
int plot_best_sheets(const char *fin_name) {return 0;}

int main(int ac, char **av) {

  if(ac<2) { cerr << "first arg is FIT or PLOT or SHEET" << endl; return 1; }

  const int nAttempt    = (ac > 2) ? atoi(av[2]) : 1;
  const char *file_name = (ac > 3) ? av[3] : "/tmp/test.root";
  const int pars_all    = (ac > 4) ? atoi(av[4]) : 2400;

  if (strcmp(av[1], "FIT") == 0) {
    return fit_data(nAttempt, file_name, pars_all);
  } else if (strcmp(av[1], "PLOT") == 0) {
    return 0;  // plot_polar(file_name,nAttempt);
  } else if (strcmp(av[1], "BAND") == 0) {
    return 0;  // get_error_band(nAttempt,file_name);
  } else { cerr << "first arg is FIT or SHEET or PLOT" << endl; return 1; }
  
  return 0;
}

double gSth, gMth;
TH1D *gh1, *gh2, *gh3, *gh12, *gh13, *gh23;
TMultiGraph *gm1, *gm2, *gm3, *gm12, *gm13, *gm23;
vector< pair<double, double> > glookup_rho_RHO, glookup_rho_F2;
int gNPAR, gAPAR_FIT1, gAPAR_FIT2, gAPAR_FIT3;

MStructureHolder *gStrHolder;
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
double gNorm1, gNorm2, gNorm3;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
TMultiGraph *split_data(const TH1D* h, double low, double up);
template <typename Type> void adjust_plot(const Type &obj);

int build_the_structure() {

  //**************** Load the data *****************//
  TGraph2D g1("data/f2piS.txt");  g1.SetName("g1");
  TGraph2D g2("data/f2piD.txt");  g2.SetName("g2");
  TGraph2D g3("data/rhopiF.txt"); g3.SetName("g3");
  TGraph2D g21("data/phi_f2piD_f2piS.txt");  g21.SetName("g21");
  TGraph2D g31("data/phi_rhopiF_f2piS.txt"); g31.SetName("g31");
  TGraph2D g32("data/phi_rhopiF_f2piD.txt"); g32.SetName("g32");
  const int nBins = 100;
  TH1D *h1 = new TH1D("h1", "2^{-+}0^{+} f_{2}#pi S;M_{3#pi}", nBins, 0.5, 2.5);
  TH1D *h2 = new TH1D("h2", "2^{-+}0^{+} f_{2}#pi D;M_{3#pi}", nBins, 0.5, 2.5);
  TH1D *h3 = new TH1D("h3", "2^{-+}0^{+} #rho#pi F;M_{3#pi}", nBins, 0.5, 2.5);
  TH1D *h12 = new TH1D("phi_f2piD_f2piS", "#Delta(f_{2}#pi S,f_{2}#pi D);M_{3#pi}", nBins, 0.5, 2.5);
  TH1D *h13 = new TH1D("phi_rhopiF_f2piS", "#Delta(#rho#pi S,#rho#pi F);M_{3#pi}", nBins, 0.5, 2.5);
  TH1D *h23 = new TH1D("phi_rhopiF_f2piD", "#Delta(#rho#pi D,#rho#pi F);M_{3#pi}", nBins, 0.5, 2.5);
  if (g1.GetN() != h1->GetNbinsX()) {cout << "Very strange!" << endl; return 0;}
  for (int i=1; i <= nBins; i++) {
    h1->SetBinContent(i, g1.GetY()[i-1]);
    h1->SetBinError  (i, g1.GetZ()[i-1]);
  }
  for (int i=1; i <= nBins; i++) {
    h2->SetBinContent(i, g2.GetY()[i-1]);
    h2->SetBinError  (i, g2.GetZ()[i-1]);
  }
  for (int i=1; i <= nBins; i++) {
    h3->SetBinContent(i, g3.GetY()[i-1]);
    h3->SetBinError  (i, g3.GetZ()[i-1]);
  }
  // phase
  for (int i=1; i <= nBins; i++) {
    h12->SetBinContent(i, normPhase(-g21.GetY()[i-1]*M_PI/180));
    double err = g21.GetZ()[i-1]*M_PI/180;
    h12->SetBinError(i, err > M_PI ? M_PI : err);
  }
  for (int i=1; i <= nBins; i++) {
    h13->SetBinContent(i, normPhase(-g31.GetY()[i-1]*M_PI/180));
    double err = g31.GetZ()[i-1]*M_PI/180;
    h13->SetBinError(i, err > M_PI ? M_PI : err);
  }
  for (int i=1; i <= nBins; i++) {
    h23->SetBinContent(i, normPhase(-g32.GetY()[i-1]*M_PI/180));
    double err = g32.GetZ()[i-1]*M_PI/180;
    h23->SetBinError(i, err > M_PI ? M_PI : err);
  }

  gh1  = h1;  gh1 ->SetStats(kFALSE);
  gh2  = h2;  gh2 ->SetStats(kFALSE);
  gh3  = h3;  gh3 ->SetStats(kFALSE);
  gh12 = h12; gh12->SetStats(kFALSE);
  gh23 = h23; gh23->SetStats(kFALSE);
  gh13 = h13; gh13->SetStats(kFALSE);

  // phase
  gh12->Scale(M_PI/180.);
  gh13->Scale(M_PI/180.);
  gh23->Scale(M_PI/180.);

  //***************** look up tables ***************//
  gMth = 3*PI_MASS;
  gSth = gMth*gMth;

  const int Nrho = 1000;
  const double Mhg = 3.0;
  for (int i=0; i < Nrho; i++) {
    double mass = gMth + (Mhg - gMth)/(Nrho-1)*i;
    glookup_rho_RHO.push_back(make_pair(mass, rhoRHO(mass*mass) ) );
    glookup_rho_F2 .push_back(make_pair(mass, rhoF2 (mass*mass) ) );
  }

  // Model
  gStrHolder = new MStructureHolder();
  // Intensities
  gStrHolder->AddWave(*h1, model_wave1, 24, glookup_rho_F2);  gStrHolder->SetWaveRange(0, 1.4, 1.95);  // 2.4);
  gStrHolder->AddWave(*h2, model_wave2, 24, glookup_rho_F2);  gStrHolder->SetWaveRange(1, 1.6, 2.1);  // 2.3);
  gStrHolder->AddWave(*h3, model_wave3, 24, glookup_rho_RHO); gStrHolder->SetWaveRange(2, 1.3, 2.1);  // 2.3);
  // Interference
  gStrHolder->AddInterference(*h12, 0, 1, phase); gStrHolder->SetInterfRange(0, 1.6, 1.95);  // 2.4);
  gStrHolder->AddInterference(*h13, 0, 2, phase); gStrHolder->SetInterfRange(1, 1.4, 1.95);  // 2.4);
  gStrHolder->AddInterference(*h23, 1, 2, phase); gStrHolder->SetInterfRange(2, 1.6, 1.95);  // 2.4);

  gStrHolder->JustToPlot(0, 1);

  //Normalise intensities
  gNorm1 = gh1 ->GetBinContent(gh1->GetMaximumBin());
  gNorm2 = gh2 ->GetBinContent(gh2->GetMaximumBin());
  gNorm3 = gh3 ->GetBinContent(gh3->GetMaximumBin());
  gh1 ->Scale(1./gNorm1);
  gh2 ->Scale(1./gNorm2);
  gh3 ->Scale(1./gNorm3);
  // separate data
  // cout << " ************ " << endl;
  gm1  = split_data(gh1, gStrHolder->GetWaveLowRange(0), gStrHolder->GetWaveUpRange(0));
  gm2  = split_data(gh2, gStrHolder->GetWaveLowRange(1), gStrHolder->GetWaveUpRange(1));
  gm3  = split_data(gh3, gStrHolder->GetWaveLowRange(2), gStrHolder->GetWaveUpRange(2));
  gm12 = split_data(gh12, gStrHolder->GetInterfLowRange(0), gStrHolder->GetInterfUpRange(0));
  gm13 = split_data(gh13, gStrHolder->GetInterfLowRange(1), gStrHolder->GetInterfUpRange(1));
  gm23 = split_data(gh23, gStrHolder->GetInterfLowRange(2), gStrHolder->GetInterfUpRange(2));
  // cout << " ------------ " << endl;
  // gm1->Print();


  return 0;
}

template <typename Type> void adjust_plot(Type &obj/*,const char *fulltitle*/) {
  obj.SetLabelSize(0.06,"xy");
  obj.GetXaxis()->SetTitleSize(0.07);
  obj.GetXaxis()->SetTitleOffset(-0.5);
  //obj.SetTitle(fulltitle);
}

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

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////    ////   ///     //////////////////////////////////////////////////////////////////////////////
///// //////// ////// ////////////////////////////////////////////////////////////////////////////////
/////   ////// ////// ////////////////////////////////////////////////////////////////////////////////
///// ///////   ///// ////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int fit_data(int nAttempt, const char *fout_name, int pars_all) {

  cout << "------------ the model with " << pars_all << endl;

  gNPAR      = 24;//(pars_all % 100000 ) / 10000;
  gAPAR_FIT1 = (pars_all % 1000   ) / 100;
  gAPAR_FIT2 = (pars_all % 100    ) / 10;
  gAPAR_FIT3 = (pars_all % 10     ) / 1;

  //****************create the model****************//
  build_the_structure();

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","plots",0,0,900,600);
  can->Divide(3,3,0.0001,0.0001);
  //can->cd(9)->SetLogy();

  //**********************Seed**********************//
  gRandom->SetSeed(0);/*12314*/ cout << "The seed " << gRandom->GetSeed() << " is used." << endl;

  //***************** Fit the data *****************//
    
  const int NumPar = 24; 
  const double step = 0.01;
  
  // Set the limites for the variables
  string var_name[NumPar] = {"mass1","mass2","mass3",
			     "g11","g12","g13",
			     "g21","g22","g23",
			     "g31","g32","g33",
			     "a11","a12","a13","a14",
			     "a21","a22","a23","a24",
			     "a31","a32","a33","a34"};

  double start_pars[NumPar] = {1.655,  1.875, 2.105,//2 
			       3.0,  1.0, 3.0,  //5
			       1.0,  3.0, 3.0,  //8
			       0.1,  3.1, 5.1,  //11
			       0.15,  0.0, 0.0, 0.0,  //15
			       0.15,  0.0, 0.0, 0.0,  //19
			       0.1,  0.0, 0.0, 0.0}; //23
  double up_limit[NumPar] = {1.73, 1.9, 2.2, 
			     3.0,  1.0, 3.0,
			     1.0,  3.0, 3.0,
			     3.0,  3.0, 3.0,
			     0.5,  0.0, 0.0, 0.0,
			     0.5,  0.0, 0.0, 0.0,
			     0.5,  0.0, 0.0, 0.0};
  double low_limit[NumPar] = {1.55, 1.8, 2.0, //2
			      1.5,  0.0, 1.5, //5
			      0.0,  1.5, 1.5, //8
			      1.5,  1.5, 1.5, //11
			      0.3,  0.0, 0.0, 0.0,  //15
			      0.3,  0.0, 0.0, 0.0,  //19
			      0.3,  0.0, 0.0, 0.0}; //23
  
  vector<int> unused = {  2,           9,10,11,   14,15,     18,19,     22,23};
  const int NstepsInfit = 5;
  vector<int> fix[NstepsInfit] = { {0,1,  3,4,5,6,7,8,     13,  17,  21      },
				   {0,1,                   13,  17,  21      },
				   {      3,4,5,6,7,8                        },
				   {0,1,                                     }};
  vector<int> int_exclude[NstepsInfit] = { {},{},{},{},{}};
  vector<int> phi_exclude[NstepsInfit] = { {0,2},{0,2},{0,2},{},{}};

  //Do minimisation several times
  TFile fout(fout_name,"RECREATE");
  TTree t("mins","several_mimimas");
  double final_pars[NumPar];
  for(int i=0;i<NumPar;i++) t.Branch(var_name[i].c_str(),&final_pars[i]);
  t.Branch("canva","TCanvas",&can);
  double chi2; t.Branch("chi2",   &chi2);
  int status; t.Branch("status",&status);
  TGraphErrors *gr[3][3], *gr_whole[3][3], *gphi[3], *gphi_whole[3];
  for(int w=0;w<3;w++) for(int i=0;i<3;i++) {gr[w][i] = 0; gr_whole[w][i] = 0; }
  for(int w=0;w<3;w++) {gphi[w] = 0; gphi_whole[w] = 0; }
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
    min->SetPrintLevel(3);
    min->Options().Print();
    
    // Create funciton wrapper for minmizer a IMultiGenFunction type 
    ROOT::Math::Functor functor(&globChi2,NumPar);
    min->SetFunction(functor);
    
    for(int i=0;i<NumPar;i++) min->SetVariable( i,var_name[i].c_str(), start_pars[i], step );
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,(up_limit[i]-low_limit[i])/1e3);
    
    //step
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,(up_limit[i]-low_limit[i])/1e3);
    //starting point
//    for(int i=0;i<NumPar;i++)  
//      start_pars[i] = low_limit[i] + (up_limit[i]-low_limit[i])*gRandom->Rndm();
    for(int i=0;i<unused.size();i++) {
      start_pars[unused[i]] = 0;
      min->FixVariable(unused[i]);
    }
    min->SetVariableValues(start_pars);
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
    
    //min->Minimize();
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
    //scalePars(arr_pars);
    setPars(arr_pars);
    //intesities
    for(int w=0;w<3;w++) 
      for(int i=0;i<3;i++) {
	gr[w][i] = gStrHolder->GetWavePlot  (w,i, 0,0, gr[w][i] );
	if(i==0) gr[w][i]->SetLineColor(kRed);
	if(i==1) gr[w][i]->SetLineColor(kBlue);
	if(i==2) gr[w][i]->SetLineColor(kGreen);
      }
    for(int w=0;w<3;w++) 
      for(int i=0;i<3;i++) {
	gr_whole[w][i] = gStrHolder->GetWavePlot  (w,i, 0.5,2.5, gr_whole[w][i] );
	gr_whole[w][i]->SetLineStyle(2);
	gr_whole[w][i]->SetLineWidth(1);
	if(i==0) gr_whole[w][i]->SetLineColor(kRed);
	if(i==1) gr_whole[w][i]->SetLineColor(kBlue);
	if(i==2) gr_whole[w][i]->SetLineColor(kGreen);
      }
    //interferences
    for(int w=0;w<3;w++)  {
      gphi[w] = gStrHolder->GetInterfPlot(w, 0,0, gphi[w]); gphi[w]->SetLineColor(kRed);
      gphi_whole[w] = gStrHolder->GetInterfPlot(w, 0.5,2.5, gphi_whole[w]); gphi_whole[w]->SetLineStyle(2); gphi_whole[w]->SetLineColor(kRed);
    }
    //plot
    //intensities
    TMultiGraph m[6];
    m[0].Add((TMultiGraph*)gm1->Clone()); m[1].Add((TMultiGraph*)gm2->Clone()); m[2].Add((TMultiGraph*)gm3->Clone());
    m[0].SetTitle(gm1->GetTitle()); m[1].SetTitle(gm2->GetTitle()); m[2].SetTitle(gm3->GetTitle());
    for(int w=0;w<3;w++) {
      for(int i=2;i>=0;i--) {
	m[w].Add((TGraphErrors*)gr[w][i]->Clone(),"lz");
	m[w].Add((TGraphErrors*)gr_whole[w][i]->Clone(),"lz");	
      }
      int Npad = 1+w+w*3;
      can->cd(Npad); m[w].Draw("a");
      adjust_plot(*m[w].GetHistogram());
    }
    //interferences
    m[3].Add((TMultiGraph*)gm12->Clone()); m[4].Add((TMultiGraph*)gm13->Clone()); m[5].Add((TMultiGraph*)gm23->Clone());
    m[3].SetTitle(gm12->GetTitle()); m[4].SetTitle(gm13->GetTitle()); m[5].SetTitle(gm23->GetTitle());
    int count = 0;
    for(int wi=0;wi<3;wi++) {
      for(int wj=0;wj<3;wj++) {
	if(wj<=wi) continue;
	m[3+count].Add((TGraphErrors*)gphi[count]->Clone(),"lz"); 
	m[3+count].Add((TGraphErrors*)gphi_whole[count]->Clone(),"lz");
	int Npad = wi*3+wj+1;
	can->cd(Npad); m[3+count].Draw("a");
	adjust_plot(*m[3+count].GetHistogram());
	count++;
      }
    }
    //save
    if(nAttempt==1) {cout << "chi2 = " << chi2 << endl; can->SaveAs("c1.pdf");}
    t.Fill();
    delete min;
  }
  
  t.Write();
  fout.Close();
  
  cout << "finished for seed " << gRandom->GetSeed() << endl;

  return 0;
}  

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

double rhoF2(double s) {
  //if(s>1e3) return 1./(8*M_PI)*(1.0-(pow(RHO_MASS,2)+pow(PI_MASS,2))/s);
  MIsobar RhoPi(F2_MASS,F2_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  double value = RhoPi.rho3(s);
  //cout << "double rhoRHO("<<s<<") " << endl;
  return value;
}

double rhoRHO(double s) {
  //if(s>1e3) return 1./(8*M_PI)*(1.0-(pow(RHO_MASS,2)+pow(PI_MASS,2))/s);
  MIsobar RhoPi(RHO_MASS,RHO_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  double value = RhoPi.rho3(s);
  //cout << "double rhoRHO("<<s<<") " << endl;
  return value;
}

////////////////////////////////////////////////////////////////////////////////////////
double glob_last_s;
boost::numeric::ublas::vector<cd> glob_akeeper;
bool glob_recalc_flag;
////////////////////////////////////////////////////////////////////////////////////////

void setPars(const double *par) {
  //To write
  for(int i=0;i<gStrHolder->pars[0].size();i++)
    if(gStrHolder->pars[0][i]!=par[i]) {glob_recalc_flag = true; break;};

  if(!glob_recalc_flag) return;

  gStrHolder->pars[0] = vector<double>(par,par+gNPAR);
  gStrHolder->pars[1] = vector<double>(par,par+gNPAR);
  gStrHolder->pars[2] = vector<double>(par,par+gNPAR);
}

double globChi2(const double *par) {  
  setPars(par);
  double chi2 = gStrHolder->GetChi2();
  return chi2;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
/* f0 pi P-wave */
cd model_wave1(int flag, double s, const double *pars) {
  if(s!=glob_last_s) glob_recalc_flag=true;

  if(flag==0) if(!glob_recalc_flag) return glob_akeeper(0); else return calculate_vector(s,pars)(0);
  double reduced_pars[gNPAR]; memcpy(reduced_pars,pars,gNPAR*sizeof(double));
  if(flag==1) { for(int i=0;i<3;i++) {reduced_pars[6+i] = 0; reduced_pars[9+i] = 0;} cd A = calculate_vector(s,reduced_pars)(0); glob_recalc_flag = true; return A;}
  if(flag==2) { for(int i=0;i<3;i++) {reduced_pars[3+i] = 0; reduced_pars[9+i] = 0;} cd A = calculate_vector(s,reduced_pars)(0); glob_recalc_flag = true; return A;}
  if(flag==3) { for(int i=0;i<3;i++) {reduced_pars[3+i] = 0; reduced_pars[6+i] = 0;} cd A = calculate_vector(s,reduced_pars)(0); glob_recalc_flag = true; return A;}

  cerr << "Chech the code. Something is going wrong" << endl;  return 0;  
}

////////////////////////////////////////////////////////////////////////////////////////
/* rho pi S-wave */
cd model_wave2(int flag, double s, const double *pars) {  
  if(s!=glob_last_s) glob_recalc_flag=true;

  if(flag==0) if(!glob_recalc_flag) return glob_akeeper(1); else return calculate_vector(s,pars)(1);
  double reduced_pars[gNPAR]; memcpy(reduced_pars,pars,gNPAR*sizeof(double));
  if(flag==1) { for(int i=0;i<3;i++) {reduced_pars[6+i] = 0; reduced_pars[9+i] = 0;} cd A = calculate_vector(s,reduced_pars)(1); glob_recalc_flag = true; return A;}
  if(flag==2) { for(int i=0;i<3;i++) {reduced_pars[3+i] = 0; reduced_pars[9+i] = 0;} cd A = calculate_vector(s,reduced_pars)(1); glob_recalc_flag = true; return A;}
  if(flag==3) { for(int i=0;i<3;i++) {reduced_pars[3+i] = 0; reduced_pars[6+i] = 0;} cd A = calculate_vector(s,reduced_pars)(1); glob_recalc_flag = true; return A;}

  cerr << "Chech the code. Something is going wrong" << endl;  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
/* rho pi D-wave */
cd model_wave3(int flag, double s, const double *pars) {
  if(s!=glob_last_s) glob_recalc_flag=true;

  if(flag==0) if(!glob_recalc_flag) return glob_akeeper(2); else return calculate_vector(s,pars)(2);
  double reduced_pars[gNPAR]; memcpy(reduced_pars,pars,gNPAR*sizeof(double));
  if(flag==1) { for(int i=0;i<3;i++) {reduced_pars[6+i] = 0; reduced_pars[9+i] = 0;} cd A = calculate_vector(s,reduced_pars)(2); glob_recalc_flag = true; return A;}
  if(flag==2) { for(int i=0;i<3;i++) {reduced_pars[3+i] = 0; reduced_pars[9+i] = 0;} cd A = calculate_vector(s,reduced_pars)(2); glob_recalc_flag = true; return A;}
  if(flag==3) { for(int i=0;i<3;i++) {reduced_pars[3+i] = 0; reduced_pars[6+i] = 0;} cd A = calculate_vector(s,reduced_pars)(2); glob_recalc_flag = true; return A;}

  cerr << "Chech the code. Something is going wrong" << endl;  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

boost::numeric::ublas::vector<cd> calculate_vector(cd s, const double *pars) {

    boost::numeric::ublas::matrix<cd> T = calculate_matrix(s, pars);
    boost::numeric::ublas::vector<cd> alpha(3); 
    alpha(0) = cd(pars[12],pars[13])+s*cd(pars[14],pars[15]);
    alpha(1) = cd(pars[16],pars[17])+s*cd(pars[18],pars[19]);
    alpha(2) = cd(pars[20],pars[21])+s*cd(pars[22],pars[23]);
    boost::numeric::ublas::vector<cd> A = prod(T,alpha);

    glob_last_s = real(s);
    glob_akeeper = A; glob_recalc_flag = false;
    //cout << "Done " << endl;
    return A;
}

boost::numeric::ublas::matrix<cd> calculate_matrix(cd s, const double *par) {
  //std::cout << "s =  " << s << std::endl;

  using namespace boost::numeric::ublas;

  double m1 = par[0];
  double m2 = par[1];
  double m3 = par[2];
  //couplings to I
  double g1[3] = {par[3],par[4],par[5]};
  //couplings to II
  double g2[3] = {par[6],par[7],par[8]};
  //couplings to III
  double g3[3] = {par[9],par[10],par[11]};
  //std::cout << "m3:  " << m3 << std::endl;
  //std::cout << "g3:  " << g3[0] << g3[1] << g3[2] << std::endl;

  //fill the firsr matrix
  symmetric_matrix<cd, upper> km1 (3, 3);
  for (unsigned i = 0; i < km1.size1 (); ++ i)
    for (unsigned j = i; j < km1.size2 (); ++ j)
      km1 (i, j) = g1[i]*g1[j];
  km1 *= 1./(m1*m1-s);
  //std::cout << "km1, " << km1 << std::endl;

  //fill the second matrix
  symmetric_matrix<cd, upper> km2 (3, 3);
  for (unsigned i = 0; i < km2.size1 (); ++ i)
    for (unsigned j = i; j < km2.size2 (); ++ j)
      km2 (i, j) = g2[i]*g2[j];
  km2 *= 1./(m2*m2-s);
  //std::cout << "km2, " << km2 << std::endl;

  //fill the second matrix
  symmetric_matrix<cd, upper> km3 (3, 3);
  for (unsigned i = 0; i < km3.size1 (); ++ i)
    for (unsigned j = i; j < km3.size2 (); ++ j)
      km3 (i, j) = g3[i]*g3[j];
  km3 *= 1./(m3*m3-s);
  //std::cout << "km3, " << km3 << std::endl;

  //add to km
  symmetric_matrix<cd, upper> km = km1+km2+km3;
  //std::cout << "km, " << km << std::endl;
  //possibly add some background
  //for (unsigned i = 0; i < km.size1 (); ++ i)
  //  for (unsigned j = i; j < km.size2 (); ++ j)
  //    km (i, j) += bgrd;
  //std::cout << km << std::endl; 

  //ph.sp.matrix
  symmetric_matrix<cd, upper> mrho (3, 3);
  //possibly add some background
  for (unsigned i = 0; i < mrho.size1 (); ++ i)
    for (unsigned j = i; j < mrho.size2 (); ++ j)
      mrho (i, j) = 0.0;
  mrho(0,0) = (imag(s)<1e-6) ? getvalue(sqrt(real(s)), glookup_rho_F2) : rhoF2(real(s));
  mrho(1,1) = (imag(s)<1e-6) ? getvalue(sqrt(real(s)), glookup_rho_F2) : rhoF2(real(s));
  mrho(2,2) = (imag(s)<1e-6) ? getvalue(sqrt(real(s)), glookup_rho_F2) : rhoRHO(real(s));

  //inverse and and phas
  identity_matrix<cd> uni (3); 
  
  matrix<cd> irhoK = cd(0,1)*prod(mrho,km);
  //std::cout << "\nirhoK " << irhoK << std::endl;
  matrix<cd> din = uni - irhoK;
  //std::cout << "\nmatrix " << din << std::endl;
  //matrix<cd> din_inv(3,3); InvertMatrix(din,din_inv);
  bool sing = false;
  matrix<cd> din_inv = gjinverse(din,sing);
  if(sing) {cerr << "ERROR: SINGULAR" << endl;}// exit(); }
  //std::cout << "inv matrix " << din_inv << std::endl;
  matrix<cd> T = prod(km,din_inv);

  return T;
  //return km;
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);
  if(a<-M_PI) return normPhase(a+2*M_PI);
  return a;
}

template <typename Type>
Type getvalue(double M, vector< pair<double,Type> > &table) {
  const int N = table.size();
  const double lft = table[0].first;
  const double rht = table[N-1].first;
  const double Mstep = table[1].first - lft;
  const int Nsteps = (M - lft)/Mstep;
  if(Nsteps<0 || Nsteps>=N-1) {cerr<<"Error!! in getvalue! M = "<<M<<", "<<table[0].second<<endl; return 0;}
  const Type value = table[Nsteps].second + 
    ( table[Nsteps+1].second - table[Nsteps].second ) / 
    ( table[Nsteps+1].first  - table[Nsteps].first  ) * (M - table[Nsteps].first); 
  return value;
}

void scalePars(double *start_pars) { 
  for(int i=0;i<4;i++) start_pars[12+i] *= sqrt(gNorm1);
  for(int i=0;i<4;i++) start_pars[16+i] *= sqrt(gNorm2);
  for(int i=0;i<4;i++) start_pars[20+i] *= sqrt(gNorm3);
}


template<class T> bool InvertMatrix(const boost::numeric::ublas::matrix<T>& input, 
				    boost::numeric::ublas::matrix<T>& inverse) {

  using namespace boost::numeric::ublas;

  typedef permutation_matrix<std::size_t> pmatrix;
  
  // create a working copy of the input
  matrix<T> A(input);
  
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());
  
  // perform LU-factorization
  int res = lu_factorize(A, pm);
  if (res != 0) return false;
  
  // create identity matrix of "inverse"
  inverse.assign(identity_matrix<T> (A.size1()));
  
  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);
  
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
