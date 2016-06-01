//09.2015 //September 2015
//Author: Misha Mikhasenko

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <initializer_list>

#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TVirtualFitter.h>
#include <TMultiGraph.h>
#include <TArrow.h>
#include <TLine.h>

#include <Math/MinimizerOptions.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include <MIsobar.h>
#include <MCoupledChannelIsobar.h>
#include <MStructureHolder.h>

#include <constants.h>
#include <deflib.h>

using namespace std;

//functions fit
void setPars(const double *par);
void scalePars(double *start_pars);
double globChi2(const double *par);
double normPhase(double a);
double phase(cd a1,cd a2) { return normPhase(arg(a1*conj(a2))); };
template <typename Type> Type getvalue(double M, vector< pair<double,Type> > &table);
//rho
double rhoRHO(double s);
double rhoF0(double s);

cd model_wave1(int flag, double s, const double *pars);
cd model_wave2(int flag, double s, const double *pars);
cd model_wave3(int flag, double s, const double *pars);

cd a1_bw(double s, double m1, double g, vector<pair<double,double> > &table, int l=0, double R = 5.0);
double deck1(double s, double b, double slope, vector<pair<double,double> > &table);

//functions main
int fit_data(int nAttempt, const char *fout_name);
int plot_best_result(const char *fin_name) {return 0;};
int plot_best_sheets(const char *fin_name) {return 0;};
int get_error_band(int nAttempt, const char *fout_name);
int plot_with_error_band(const char *fin_name) {return 1;};
int plot_polar(const char *fin_name, int entry);

int main(int ac, char **av) {

  if(ac<2) { cerr << "first arg is FIT or PLOT or SHEET" << endl; return 1; }

  const int nAttempt    = (ac > 2) ? atoi(av[2]) : 1;
  const char *file_name = (ac > 3) ? av[3] : "/tmp/test.root";

  if(strcmp(av[1],"FIT")==0) {
    return fit_data(nAttempt,file_name);
  } else if(strcmp(av[1],"PLOT")==0) {
    return plot_polar(file_name,nAttempt);
  } else if(strcmp(av[1],"BAND")==0) {
    return get_error_band(nAttempt,file_name);
  } else { cerr << "first arg is FIT or SHEET or PLOT" << endl; return 1; }
  
  return 0;
}

double gSth,gMth;
TH1D *gh1, *gh2, *gh3, *gh12, *gh13, *gh23;
TMultiGraph *gm1, *gm2, *gm3, *gm12, *gm13, *gm23;
vector< pair<double,double> > glookup_rho_RHO, glookup_rho_F0;
vector< pair<double,cd> > glookup_rhopiTof0pi, glookup_KstarKTof0pi;

MStructureHolder *gStrHolder;
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
double gNorm1,gNorm2,gNorm3;
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
TMultiGraph *split_data(const TH1D* h, double low, double up);
template <typename Type> void adjust_plot(Type &obj);

int build_the_structure() {

  //**************** Load the data *****************//
  TGraph2D g1 ("data/oliver/f0piP.oliver.txt");  g1.SetName("g1"); 
  TGraph2D g2 ("data/oliver/rhopiS.oliver.txt"); g2.SetName("g2");  
  TGraph2D g3 ("data/oliver/rhopiD.oliver.txt"); g3.SetName("g3");  
  TGraph2D g12("data/oliver/phi_f0piP_rhopiS.oliver.txt");  g12.SetName("g12");
  TGraph2D g31("data/oliver/phi_rhopiD_f0piP.oliver.txt");  g31.SetName("g31");
  TGraph2D g32("data/oliver/phi_rhopiD_rhopiS.oliver.txt"); g32.SetName("g32");
  const int nBins = 100;
  TH1D *h1 = new TH1D("h1","1^{++}0^{+} f_{0}#pi P;M_{3#pi}",nBins,0.5,2.5);
  TH1D *h2 = new TH1D("h2","1^{++}0^{+} #rho#pi S;M_{3#pi}",nBins,0.5,2.5);
  TH1D *h3 = new TH1D("h3","2^{++}1^{+} #rho#pi D;M_{3#pi}",nBins,0.5,2.5);
  TH1D *h12 = new TH1D("phi_f0piP_rhopiS","#Delta(f_{0}#pi P, #rho#pi S);M_{3#pi}",nBins,0.5,2.5);
  TH1D *h13 = new TH1D("phi_f0piP_rhopiD","#Delta(f_{0}#pi P, #rho#pi D);M_{3#pi}",nBins,0.5,2.5);
  TH1D *h23 = new TH1D("phi_rhopiS_rhopiD","#Delta(#rho#pi S, #rho#pi D);M_{3#pi}",nBins,0.5,2.5);
  if(g1.GetN() != h1->GetNbinsX()) {cout << "Very strange!" << endl; return 0;} 
  for(int i=1;i<=nBins;i++) {
    h1->SetBinContent(i,g1.GetY()[i-1]);
    h1->SetBinError  (i,g1.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    h2->SetBinContent(i,g2.GetY()[i-1]);
    h2->SetBinError  (i,g2.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    h3->SetBinContent(i,g3.GetY()[i-1]);
    h3->SetBinError  (i,g3.GetZ()[i-1]);
  }
  //phase
  for(int i=1;i<=nBins;i++) {
    h12->SetBinContent(i,g12.GetY()[i-1]);
    h12->SetBinError  (i,g12.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    h13->SetBinContent(i,-g31.GetY()[i-1]);
    h13->SetBinError  (i,g31.GetZ()[i-1]);
  }
  for(int i=1;i<=nBins;i++) {
    h23->SetBinContent(i,-g32.GetY()[i-1]);
    h23->SetBinError  (i,g32.GetZ()[i-1]);
  }

  gh1  = h1;  gh1 ->SetStats(kFALSE);
  gh2  = h2;  gh2 ->SetStats(kFALSE);
  gh3  = h3;  gh3 ->SetStats(kFALSE);
  gh12 = h12; gh12->SetStats(kFALSE);
  gh23 = h23; gh23->SetStats(kFALSE);
  gh13 = h13; gh13->SetStats(kFALSE);

  //phase
  gh12->Scale(M_PI/180.);
  gh13->Scale(M_PI/180.);
  gh23->Scale(M_PI/180.);

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
    glookup_KstarKTof0pi.push_back( make_pair( gr_tr1.GetX()[i], 
					      cd(gr_tr1.GetY()[i],gr_tr1.GetZ()[i]) ) );

  TGraph2D gr_tr2("calculations/sclr_trgl_rhopiTof0pi.out"); gr_tr2.SetName("tr2");
  for(int i=0;i<gr_tr2.GetN();i++)
    glookup_rhopiTof0pi.push_back( make_pair( gr_tr2.GetX()[i], 
					       cd(gr_tr2.GetY()[i],gr_tr2.GetZ()[i]) ) );  
  //normalizetion
  double max;
  max = abs( glookup_rhopiTof0pi[0].second );
  for(int i=1;i<glookup_rhopiTof0pi .size();i++) if(abs(glookup_rhopiTof0pi [i].second)>max) max = abs(glookup_rhopiTof0pi [i].second); 
  for(int i=0;i<glookup_rhopiTof0pi .size();i++) glookup_rhopiTof0pi[i].second *= 1./max;
  max = abs( glookup_KstarKTof0pi[0].second );
  for(int i=1;i<glookup_KstarKTof0pi.size();i++) if(abs(glookup_KstarKTof0pi[i].second)>max) max = abs(glookup_KstarKTof0pi[i].second); 
  for(int i=0;i<glookup_KstarKTof0pi.size();i++) glookup_KstarKTof0pi[i].second *= 1./max;  

  //Model
  gStrHolder = new MStructureHolder();
  //intensities
  gStrHolder->AddWave(*h1,model_wave1,17,glookup_rho_F0 ); gStrHolder->SetWaveRange(0, 1.3,1.6);
  gStrHolder->AddWave(*h2,model_wave2,8 ,glookup_rho_RHO); gStrHolder->SetWaveRange(1, 0.8,2.0);
  gStrHolder->AddWave(*h3,model_wave3,4 ,glookup_rho_RHO); gStrHolder->SetWaveRange(2, 1.0,2.0);
  //interference  
  gStrHolder->AddInterference(*h12,0,1,phase); gStrHolder->SetInterfRange(0, 1.3,1.7);
  gStrHolder->AddInterference(*h13,0,2,phase); gStrHolder->SetInterfRange(1, 1.3,1.7);
  gStrHolder->AddInterference(*h23,1,2,phase); gStrHolder->SetInterfRange(2, 1.0,1.7);

  gStrHolder->JustToPlot(0,1);

  //separate data
  //cout << " ************ " << endl;
  gm1  = split_data(gh1,gStrHolder->GetWaveLowRange(0),gStrHolder->GetWaveUpRange(0));
  gm2  = split_data(gh2,gStrHolder->GetWaveLowRange(1),gStrHolder->GetWaveUpRange(1));
  gm3  = split_data(gh3,gStrHolder->GetWaveLowRange(2),gStrHolder->GetWaveUpRange(2));
  gm12 = split_data(gh12,gStrHolder->GetInterfLowRange(0),gStrHolder->GetInterfUpRange(0));
  gm13 = split_data(gh13,gStrHolder->GetInterfLowRange(1),gStrHolder->GetInterfUpRange(1));
  gm23 = split_data(gh23,gStrHolder->GetInterfLowRange(2),gStrHolder->GetInterfUpRange(2));
  //cout << " ------------ " << endl;
  //gm1->Print();

  //Normalise intensities
  gNorm1 = gh1 ->GetBinContent(gh1->GetMaximumBin());
  gNorm2 = gh2 ->GetBinContent(gh2->GetMaximumBin());
  gNorm3 = gh3 ->GetBinContent(gh3->GetMaximumBin());
  gh1 ->Scale(1./gNorm1);
  gh2 ->Scale(1./gNorm2);
  gh3 ->Scale(1./gNorm3);

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

int fit_data(int nAttempt, const char *fout_name) {

  //****************create the model****************//
  build_the_structure();

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","plots",0,0,900,600);
  can->Divide(3,3,0.0001,0.0001);
  can->cd(9)->SetLogy();

  //**********************Seed**********************//
  gRandom->SetSeed(0);/*12314*/ cout << "The seed " << gRandom->GetSeed() << " is used." << endl;

  //***************** Fit the data *****************//
    
  const int NumPar = 25;  
  const double step = 0.01;
  
  // Set the limites for the variables
  string var_name[NumPar] = {"a1_mass","a1_g","a1_rc0","a1_ic0",
			     "b1","slope1","d1_rc0","d1_ic0",
			     "tr1_rc","tr1_ic","tr2_rc","tr2_ic",
			     "b2","slope2","d2_rc0","d2_ic0",
			     "shift",
			     "a2_mass","a2_g","a2_rc0","a2_ic0",
			     "b3","slope3","d3_rc0","d3_ic0"};

  double start_pars[NumPar] = {1.2,  7, 0.07, 0.0, 
			       1.0, -2.0, 0.0, 0.0,
			       1.0,  0.0, 0.0, 0.0,
			       1.0, -2.0, 0.0, 0.0,
			       0.0,
			       1.32, 3.8, 0.07, 0.0,
			       1.0, -2.0, 0.0, 0.0};
  double  up_limit[NumPar] = {1.3, 7.0, 0.08, 0.0,
			      10.0, 0.0, .1, .1,
			      2.0, 1.0, 1.0, 1.0,
			      2.0, 0.0, 0.1, 0.1,
			      0.2,
			      1.4, 4.0, 0.08, 0.0,
			      2.0, 0.0, 0.1, 0.1};
  double low_limit[NumPar] = {1.1, 5.0,  0.02, -0.0,
			      0.0, -10.0, -.1, -.1,
			      1.0, -1.0, -1.0, -1.0,
			      0.0, -10.0, -.01, -0.1,
			      0.0,
			      1.3, 3.6, 0.02, 0.0,
			      0.0, -10.0,-0.1, -0.1};
   
  vector<int> unused = {3,12,13,14,15,16};//
  const int NstepsInfit = 3;
  vector<int> fix[NstepsInfit] = { {                  10,11,            21,22,23,24},
				   {0,1,  4,5,    8,9,10,11,                       },
				   {0,1,2,4,5,6,7,          17,18,19,20,21,22,23,24}};
  vector<int> int_exclude[NstepsInfit] = { {},
					   {0},
					   {1,2}};
  vector<int> phi_exclude[NstepsInfit] = { {0,2},
					   {0},
					   {2}};
    
  //Do minimisation several times
  TFile fout(fout_name,"RECREATE");
  TTree t("mins","several_mimimas");
  double final_pars[NumPar];
  for(int i=0;i<NumPar;i++) t.Branch(var_name[i].c_str(),&final_pars[i]);
  t.Branch("canva","TCanvas",&can);
  double chi2; t.Branch("chi2",   &chi2);
  int status; t.Branch("status",&status);
  TGraphErrors *gr[3][3], *gr_whole[3][3], *gphi[3], *gphi_whole[3];
  TGraphErrors *gr_second_tr_whole = 0, *gr_second_tr = 0;
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
    min->SetPrintLevel(1);
    min->Options().Print();
    
    // Create funciton wrapper for minmizer a IMultiGenFunction type 
    ROOT::Math::Functor functor(&globChi2,NumPar);
    min->SetFunction(functor);
    
    for(int i=0;i<NumPar;i++) min->SetVariable( i,var_name[i].c_str(), start_pars[i], step );
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,(up_limit[i]-low_limit[i])/1e3);
    
    //step
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,(up_limit[i]-low_limit[i])/1e3);
    //starting point
    for(int i=0;i<NumPar;i++)  
      start_pars[i] = low_limit[i] + (up_limit[i]-low_limit[i])*gRandom->Rndm();
    for(int i=0;i<unused.size();i++) {
      start_pars[unused[i]] = 0;
      min->FixVariable(unused[i]);
    }
    min->SetVariableValues(start_pars);
    cout<<"Initual parameters:\n";for(int i=0;i<NumPar;i++) {cout<<var_name[i]<<" = "<<start_pars[i]<<";\n";}

    //do the set of minimizations with stategy
    for(int j=0;j<NstepsInfit;j++) {//NstepsInfit
      for(int i=0;i<fix[j].size();i++) min->FixVariable    (fix[j][i]);
      for(int i=0;i<int_exclude[j].size();i++) gStrHolder->JustToPlot(int_exclude[j][i],true);
      for(int i=0;i<phi_exclude[j].size();i++) gStrHolder->JustToPlot(0,phi_exclude[j][i],true);
      gStrHolder->Print();
      min->Minimize();
      for(int i=0;i<fix[j].size();i++) min->ReleaseVariable(fix[j][i]);
      for(int i=0;i<int_exclude[j].size();i++) gStrHolder->JustToPlot(int_exclude[j][i],false);
      for(int i=0;i<phi_exclude[j].size();i++) gStrHolder->JustToPlot(0,phi_exclude[j][i],false);
    }
    
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
    scalePars(arr_pars);
    setPars(arr_pars);
    //intesities
    for(int w=0;w<3;w++) 
      for(int i=0;i<3;i++) {
	gr[w][i] = gStrHolder->GetWavePlot  (w,i, 0,0, gr[w][i] );
	if(i==0) gr[w][i]->SetLineColor(kRed);
	if(i==1) gr[w][i]->SetLineColor(kBlue);
	if(i==2) gr[w][i]->SetLineColor(kGreen);
      }
    gr_second_tr = gStrHolder->GetWavePlot  (0,3, 0,0, gr_second_tr );
    gr_second_tr->SetLineColor(9);
    for(int w=0;w<3;w++) 
      for(int i=0;i<3;i++) {
	gr_whole[w][i] = gStrHolder->GetWavePlot  (w,i, 0.5,2.5, gr_whole[w][i] );
	gr_whole[w][i]->SetLineStyle(2);
	gr_whole[w][i]->SetLineWidth(1);
	if(i==0) gr_whole[w][i]->SetLineColor(kRed);
	if(i==1) gr_whole[w][i]->SetLineColor(kBlue);
	if(i==2) gr_whole[w][i]->SetLineColor(kGreen);
      }
    gr_second_tr_whole = gStrHolder->GetWavePlot  (0,3, 0.5,2.5, gr_second_tr_whole );
    gr_second_tr_whole->SetLineStyle(2);
    gr_second_tr_whole->SetLineColor(9);
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
    m[0].Add((TGraphErrors*)gr_second_tr->Clone());
    m[0].Add((TGraphErrors*)gr_second_tr_whole->Clone());
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
/////    //////  ///// /// ///    ///////////////////////////////////////////////////////////////////
///// // ///// // ////   / /// //  //////////////////////////////////////////////////////////////////
///// /// ///     //// /   /// /// //////////////////////////////////////////////////////////////////
/////     // ////  /// //  ///    ///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

int get_error_band(int nAttempt, const char *fin_name) {

  const char* fout_name = "/tmp/test_band.root";

  TFile *fin = TFile::Open(fin_name);           if(!fin) {cerr<<"No file!"<<endl;return 1;}
  TTree *tin = (TTree*)gDirectory->Get("mins");   if(!tin)   {cerr<<"No tree!"<<endl;return 1;}

  //Find entries number with minimal chi2
  const int N = tin->GetEntries();
  double chi2; tin->SetBranchAddress("chi2",&chi2); tin->GetEntry(0);
  int emin=0;
  double chi2min = chi2;
  for(int i=1;i<N;i++) {
    if(i%100==0) std::cout << i << ", ";
    tin->GetEntry(i);
    if(chi2<chi2min) { chi2min=chi2; emin=i;}
  } 
  cout << "\n-------------> Best chi2 = " << chi2min << " ("<<emin<<")" << endl;
  
  //Get parameters first
  const int NumPar = 25;
  //List of parameters
  string var_name[NumPar] = {"a1_mass","a1_g","a1_rc0","a1_ic0",
			     "b1","slope1","d1_rc0","d1_ic0",
			     "tr1_rc","tr1_ic","tr2_rc","tr2_ic",
			     "b2","slope2","d2_rc0","d2_ic0",
			     "shift",
			     "a2_mass","a2_g","a2_rc0","a2_ic0",
			     "b3","slope3","d3_rc0","d3_ic0"};
  vector<int> unused = {3,12,13,14,15,16};//
  const int NstepsInfit = 3;
  vector<int> fix[NstepsInfit] = { {                  10,11,            21,22,23,24},
				   {0,1,  4,5,    8,9,10,11,                       },
				   {0,1,2,4,5,6,7,          17,18,19,20,21,22,23,24}};
  vector<int> int_exclude[NstepsInfit] = { {},
					   {0},
					   {1,2}};
  vector<int> phi_exclude[NstepsInfit] = { {0,2},
					   {0},
					   {2}};
  //Take them
  double start_pars[NumPar];
  for(int i=0;i<NumPar;i++) tin->SetBranchAddress(var_name[i].c_str(),&start_pars[i]);
  tin->GetEntry(emin);
  build_the_structure();

  //**********************Seed**********************//
  gRandom->SetSeed(0);/*12314*/ cout << "The seed " << gRandom->GetSeed() << " is used." << endl;
  //****************create the model****************//

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","plots",0,0,900,600);
  can->Divide(3,3,0.0001,0.0001);
  can->cd(9)->SetLogy();

  //Do minimisation several times
  TFile fout(fout_name,"RECREATE");
  TTree t("mins","bootstrup");
  double final_pars[NumPar];
  for(int i=0;i<NumPar;i++) t.Branch(var_name[i].c_str(),&final_pars[i]);
  t.Branch("canva","TCanvas",&can);
  t.Branch("chi2",   &chi2);
  int status; t.Branch("status",&status);
  TGraphErrors *gr[3][3], *gr_whole[3][3], *gphi[3], *gphi_whole[3];
  TGraphErrors *gr_second_tr_whole = 0, *gr_second_tr = 0;
  for(int w=0;w<3;w++) for(int i=0;i<3;i++) {gr[w][i] = 0; gr_whole[w][i] = 0; }
  for(int w=0;w<3;w++) {gphi[w] = 0; gphi_whole[w] = 0; }
  for(int e=0;e<nAttempt;e++) {
    cout << "---------- Attempt " << e << " -----------"<< endl;
    cout << "------------------------------------------"<< endl;
    //step

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
    const double step = 0.001;
    for(int i=0;i<NumPar;i++) min->SetVariable( i,var_name[i].c_str(), start_pars[i], step );
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,step);    

    //set unused variables
    min->SetVariableValues(start_pars);
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,step);    
    for(int i=0;i<unused.size();i++) {
      min->SetVariableValue(unused[i],0);
      min->FixVariable(unused[i]);
    }
    cout<<"Initual parameters:\n";for(int i=0;i<NumPar;i++) {cout<<var_name[i]<<" = "<<start_pars[i]<<";\n";}

    //generate data
    cout<<"..............reGenerating data............." << endl;
    gStrHolder->reGenerateData();
    cout<<"..................Minimizing................" << endl;
    //do the minimization
    for(int j=0;j<NstepsInfit;j++) {//NstepsInfit
      for(int i=0;i<fix[j].size();i++) min->FixVariable    (fix[j][i]);
      for(int i=0;i<int_exclude[j].size();i++) gStrHolder->JustToPlot(int_exclude[j][i],true);
      for(int i=0;i<phi_exclude[j].size();i++) gStrHolder->JustToPlot(0,phi_exclude[j][i],true);
      gStrHolder->Print();
      min->Minimize();
      for(int i=0;i<fix[j].size();i++) min->ReleaseVariable(fix[j][i]);
      for(int i=0;i<int_exclude[j].size();i++) gStrHolder->JustToPlot(int_exclude[j][i],false);
      for(int i=0;i<phi_exclude[j].size();i++) gStrHolder->JustToPlot(0,phi_exclude[j][i],false);
    }
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
    scalePars(arr_pars);
    setPars(arr_pars);
    //intesities
    for(int w=0;w<3;w++) 
      for(int i=0;i<3;i++) {
	gr[w][i] = gStrHolder->GetWavePlot  (w,i, 0,0, gr[w][i] );
	if(i==0) gr[w][i]->SetLineColor(kRed);
	if(i==1) gr[w][i]->SetLineColor(kBlue);
	if(i==2) gr[w][i]->SetLineColor(kGreen);
      }
    gr_second_tr = gStrHolder->GetWavePlot  (0,3, 0,0, gr_second_tr );
    gr_second_tr->SetLineColor(9);
    for(int w=0;w<3;w++) 
      for(int i=0;i<3;i++) {
	gr_whole[w][i] = gStrHolder->GetWavePlot  (w,i, 0.5,2.5, gr_whole[w][i] );
	gr_whole[w][i]->SetLineStyle(2);
	gr_whole[w][i]->SetLineWidth(1);
	if(i==0) gr_whole[w][i]->SetLineColor(kRed);
	if(i==1) gr_whole[w][i]->SetLineColor(kBlue);
	if(i==2) gr_whole[w][i]->SetLineColor(kGreen);
      }
    gr_second_tr_whole = gStrHolder->GetWavePlot  (0,3, 0.5,2.5, gr_second_tr_whole );
    gr_second_tr_whole->SetLineStyle(2);
    gr_second_tr_whole->SetLineColor(9);
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
    m[0].Add((TGraphErrors*)gr_second_tr->Clone());
    m[0].Add((TGraphErrors*)gr_second_tr_whole->Clone());
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
////     ///  /////      //      ////////////////////////////////////////////////////////////////////
//// ///  //  /////  //  ////  //////////////////////////////////////////////////////////////////////
////     ///  // //  //  ////  //////////////////////////////////////////////////////////////////////
////  //////     //      ////  //////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void rotate(int N,double *ax1,double *ay1,double *ax2,double *ay2);
TMultiGraph *mgm(std::initializer_list<TGraph*> a_args);
////////////////////////////////////////////////////////////////////////////////////
int plot_polar(const char *fin_name, int entry) {

  TFile *fin = TFile::Open(fin_name);           if(!fin) {cerr<<"No file!"<<endl;return 1;}
  TTree *tin = (TTree*)gDirectory->Get("mins");   if(!tin)   {cerr<<"No tree!"<<endl;return 1;}

  //Find entries number with minimal chi2
  const int N = tin->GetEntries();

  //Find entries number with minimal chi2
  int _entry=entry;
  if(entry < 0 || entry>=N) {
    double chi2; tin->SetBranchAddress("chi2",&chi2); tin->GetEntry(0);
    int emin=0;
    double chi2min = chi2;
    for(int i=1;i<N;i++) {
      if(i%100==0) std::cout << i << ", ";
      tin->GetEntry(i);
      if(chi2<chi2min) { chi2min=chi2; emin=i;}
    } 
    cout << "\n-------------> Best chi2 = " << chi2min << " ("<<emin<<")" << endl;
    _entry = emin; tin->ResetBranchAddress(tin->GetBranch("chi2"));
  }

  //Get parameters first
  const int NumPar = 25;
  //List of parameters
  string var_name[NumPar] = {"a1_mass","a1_g","a1_rc0","a1_ic0",
			     "b1","slope1","d1_rc0","d1_ic0",
			     "tr1_rc","tr1_ic","tr2_rc","tr2_ic",
			     "b2","slope2","d2_rc0","d2_ic0",
			     "shift",
			     "a2_mass","a2_g","a2_rc0","a2_ic0",
			     "b3","slope3","d3_rc0","d3_ic0"};
  //Take them
  double start_pars[NumPar];
  for(int i=0;i<NumPar;i++) tin->SetBranchAddress(var_name[i].c_str(),&start_pars[i]);

  /*..................*/
  build_the_structure();
  /*..................*/

  const int Nentries = tin->GetEntries();
  double all_pars[Nentries][NumPar];
  vector<double*> vpar(Nentries);
  for(int i=0;i<Nentries;i++) {//Nentries
    tin->GetEntry(i);
    //set correct normanization
    scalePars(start_pars);
    memcpy(all_pars[i],start_pars,NumPar*sizeof(double));
    vpar[i] = &all_pars[i][0];
  }
  //****************create the model****************//
  setPars(vpar[_entry]);

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","Argon plot",0,0,600,600);

  //***************************plot
  //convert data to polar plot
  //amp 1
  gh1->Scale(gNorm1);
  const int Nbins = gh1->GetXaxis()->GetNbins();
  TGraphErrors *gDt1 = new TGraphErrors(Nbins);
  for(int i=0;i<Nbins;i++) {
    double mass = gh1->GetBinCenter(i+1);
    double intensity = gh1->GetBinContent(i+1);
    double phasespace = getvalue(mass, glookup_rho_F0);
    double amp = sqrt(intensity/phasespace);
    double phi = gh12->GetBinContent(i+1);
    gDt1->GetX()[i] = amp*cos(phi);
    gDt1->GetY()[i] = amp*sin(phi);      
    //errors
    double damp1x = 1./sqrt(phasespace*intensity)/2.*cos(phi)*gh1 ->GetBinError(i+1);
    double damp2x = sqrt(intensity/phasespace)/2.*sin(phi)   *gh12->GetBinError(i+1);
    double dampx = sqrt(damp1x*damp1x+damp2x*damp2x);
    double damp1y = 1./sqrt(phasespace*intensity)/2.*sin(phi)*gh1 ->GetBinError(i+1);
    double damp2y = sqrt(intensity/phasespace)/2.*cos(phi)   *gh12->GetBinError(i+1);
    double dampy = sqrt(damp1y*damp1y+damp2y*damp2y);
    gDt1->GetEX()[i] = dampx;
    gDt1->GetEY()[i] = dampy;
  }

  //change the phase of amp1
  TGraphErrors *tg1[2],*tg1_whole[2],*tg2,*tg2_whole,*tg2_100;
  //RED line
  const double mlow = gStrHolder->GetWaveLowRange(0);
  const double mup = gStrHolder->GetWaveUpRange(0);
  const double mlow_full = gh1->GetBinLowEdge(1);
  //const int Nbins = h1->GetXaxis()->GetNbins();
  const double mup_full  = gh1->GetBinLowEdge(Nbins+1);

  tg1[0] = gStrHolder->GetPolarPlot(0,0, mlow,mup, 0 );
  tg2 = gStrHolder->GetPolarPlot(1,0, mlow,mup, 0 );
  tg1[0]->SetLineColor(kRed);
  tg1_whole[0] = gStrHolder->GetPolarPlot(0,0, mlow_full,mup_full, 0 );
  tg2_whole = gStrHolder->GetPolarPlot(1,0, mlow_full,mup_full, 0 );
  tg2_100 = gStrHolder->GetPolarPlot(1,0, mlow_full,mup_full, 0, Nbins);
  tg1_whole[0]->SetLineColor(kRed);  tg1_whole[0]->SetLineStyle(2);
  //BLUE line
  tg1[1] = gStrHolder->GetPolarPlot(0,1, mlow,mup, 0 );
  tg1[1]->SetLineColor(kBlue);
  tg1_whole[1] = gStrHolder->GetPolarPlot(0,1, mlow_full,mup_full, 0 );
  tg1_whole[1]->SetLineColor(kBlue);  tg1_whole[1]->SetLineStyle(2);
  //rotate data
  rotate(gDt1->GetN(),
	 gDt1->GetX(),gDt1->GetY(),
	 tg2_100->GetX(),tg2_100->GetY());
  //separate data
  TGraphErrors *gDt1_gray = new TGraphErrors(0);
  TGraphErrors *gDt1_dark = new TGraphErrors(0);
  for(int i=0;i<Nbins;i++) {
    double mass = gh1->GetBinCenter(i+1);
    if(mass<mlow||mass>mup) {
      gDt1_gray->Set(gDt1_gray->GetN()+1);
      gDt1_gray->SetPoint(gDt1_gray->GetN()-1,gDt1->GetX()[i],gDt1->GetY()[i]);
      gDt1_gray->SetPointError(gDt1_gray->GetN()-1,gDt1->GetEX()[i],gDt1->GetEY()[i]);
    } else {
      gDt1_dark->Set(gDt1_dark->GetN()+1);
      gDt1_dark->SetPoint(gDt1_dark->GetN()-1,gDt1->GetX()[i],gDt1->GetY()[i]);
      gDt1_dark->SetPointError(gDt1_dark->GetN()-1,gDt1->GetEX()[i],gDt1->GetEY()[i]);
    }
  }
  gDt1_dark->SetLineColor(kBlack); 
  gDt1_gray->SetLineColor(kGray); gDt1_gray->SetMarkerColor(kGray);

  can->cd(1); 
  TMultiGraph *mt = mgm({tg1[0],tg1[1],tg1_whole[0],tg1_whole[1]});
  mt->Add(gDt1_gray,"p");
  mt->Add(gDt1_dark,"p");
  mt->Draw("a");// gDt1_gray->Draw("same p"); gDt1_dark->Draw("same p");
  const int Nlines = 10;
  for(int i=0;i<Nlines;i++) {
    int index( tg1[0]->GetN()/(Nlines-1)*i );
    TArrow *l = new TArrow(tg1[1]->GetX()[index],tg1[1]->GetY()[index],
			   tg1[0]->GetX()[index],tg1[0]->GetY()[index],0.01,">");
    l->SetLineColor(9); l->SetAngle(60); 
    l->Draw();
  }
  mt->GetHistogram()->SetTitle("1^{++}0^{+} f_{0}#pi P: triangle (K*KK) + triangle ( #rho#pi#pi);real(A);imag(A)");
  TLine *lx = new TLine(mt->GetHistogram()->GetXaxis()->GetXmin(),0,mt->GetHistogram()->GetXaxis()->GetXmax(),0); lx->SetLineStyle(2); lx->Draw();
  TLine *ly = new TLine(0,mt->GetHistogram()->GetYaxis()->GetXmin(),0,mt->GetHistogram()->GetYaxis()->GetXmax()); ly->SetLineStyle(2); ly->Draw();
  can->SaveAs("c1.pdf");
  
  return 0;

}

void rotate(int N,double *ax1,double *ay1,double *ax2,double *ay2) {
  for(int i=0;i<N;i++) {
    double x1 = ax1[i];
    double y1 = ay1[i];
    double phi2 = -atan2(ay2[i],ax2[i]);
    double x1_sub =  x1*cos(phi2) + y1*sin(phi2);
    double y1_sub = -x1*sin(phi2) + y1*cos(phi2) ;
    ax1[i] = x1_sub;
    ay1[i] = y1_sub;
  }
}

TMultiGraph *mgm(std::initializer_list<TGraph*> a_args) {
  TMultiGraph *_mg = new TMultiGraph();
  for (auto i: a_args) _mg->Add(i,"lz");
  return _mg;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

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

void setPars(const double *par) {
  gStrHolder->pars[0] = vector<double>(par,par+17);
  gStrHolder->pars[1] = vector<double>(par,par+8);
  gStrHolder->pars[2] = vector<double>(par+17,par+17+8);
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
  double m2 = pars[0];
  double g2 = pars[1];
  cd c2(pars[2],pars[3]);
  cd bw = c2*a1_bw(s,m2,g2,glookup_rho_RHO);

  double b2  = pars[4];
  double sl2 = pars[5];
  cd d2_c(pars[6],pars[7]);  
  double shift = pars[16];

  cd tr1_c(pars[8],pars[9]);
  cd tr1 = tr1_c*bw * getvalue(sqrt(s)-shift,glookup_KstarKTof0pi);
  if(flag==1) return tr1;

  cd tr2_c(pars[10],pars[11]);
  cd amp2 = model_wave2(0, s, pars);
  cd tr2 = tr2_c*amp2 * getvalue(sqrt(s),glookup_rhopiTof0pi);
  if(flag==3) return tr2;

  double b1  = pars[12];
  double sl1 = pars[13];
  cd d1_c(pars[14],pars[15]);  
  cd d1 = d1_c*deck1(s,b1,sl1,glookup_rho_F0);
  if(flag==2) return d1;
  if(flag==0) return tr1 + tr2 + d1;

  cerr << "Chech the code. Something is going wrong" << endl;  return 0;  
}

////////////////////////////////////////////////////////////////////////////////////////
/* rho pi S-wave */
cd model_wave2(int flag, double s, const double *pars) {  
  double m2 = pars[0];
  double g2 = pars[1];
  cd c2(pars[2],pars[3]);
  cd bw = c2*a1_bw(s,m2,g2,glookup_rho_RHO);
  if(flag==1) return bw;

  double b2  = pars[4];
  double sl2 = pars[5];
  cd d2_c(pars[6],pars[7]);  
  cd d2 = d2_c*deck1(s,b2,sl2,glookup_rho_RHO);
  if(flag==2) return d2;

  if(flag==0) return bw + d2;
  cerr << "Chech the code. Something is going wrong" << endl;  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
/* rho pi D-wave */
cd model_wave3(int flag, double s, const double *pars) {
  double m3 = pars[0];
  double g3 = pars[1];
  cd c3(pars[2],pars[3]);
  cd bw = c3*a1_bw(s,m3,g3,glookup_rho_RHO,2,5.0);
  if(flag==1) return bw;

  double b3  = pars[4];
  double sl3 = pars[5];
  cd d3_c(pars[6],pars[7]);  
  cd d3 = d3_c*deck1(s,b3,sl3,glookup_rho_RHO);
  if(flag==2) return d3;

  if(flag==0) return bw + d3;
  cerr << "Chech the code. Something is going wrong" << endl;  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double normPhase(double a) {
  if(a>M_PI)  return normPhase(a-2*M_PI);
  if(a<-M_PI) return normPhase(a+2*M_PI);
  return a;
}

cd a1_bw(double s, double m1, double g, vector<pair<double,double> > &table, int l, double R) {
  double rho  = getvalue(sqrt(s), glookup_rho_RHO);
  if(m1<gMth) return 0;
  double rho0 = (m1<2.5) ? getvalue( m1, table) : 0.1/(8*M_PI);
  double p  = 4*M_PI*sqrt(s)*rho;
  double p0 = 4*M_PI*     m1*rho0;
  double factor = pow(p*p/(p0*p0)*((1+p0*p0*R*R)/(1+p*p*R*R)),l); 
  cd value = g*g*sqrt(factor)/(m1*m1-s- 0.5*g*g*rho*cd(0,1)*factor);
  return value;
}

double deck1(double s, double b, double slope, vector<pair<double,double> > &table) {
  double rho = getvalue(sqrt(s), table);//1./(8*M_PI)*sqrt((s-pow(gMth,2))/s);
  double p = sqrt(s)/2.0*(8*M_PI)*rho;
  double value = pow(sqrt(s)-gMth,b)*exp(slope*p);
  //norm
  return value;
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
  start_pars[2] *= sqrt(gNorm2); start_pars[3] *= sqrt(gNorm2); start_pars[6] *= sqrt(gNorm2); start_pars[7] *= sqrt(gNorm2);
  start_pars[8] *= sqrt(gNorm1/gNorm2); start_pars[9] *= sqrt(gNorm1/gNorm2); start_pars[10] *= sqrt(gNorm1/gNorm2); start_pars[11] *= sqrt(gNorm1/gNorm2);
  start_pars[14] *= sqrt(gNorm1); start_pars[15] *= sqrt(gNorm1);
  start_pars[19] *= sqrt(gNorm3); start_pars[20] *= sqrt(gNorm3); start_pars[23] *= sqrt(gNorm3); start_pars[24] *= sqrt(gNorm3);
}
