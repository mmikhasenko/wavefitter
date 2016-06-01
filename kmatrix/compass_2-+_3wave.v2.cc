//09.2015 //Oktober 2015
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
#include <MmatrixK.h>

#include <constants.h>
#include <deflib.h>

#include <boost/math/special_functions/legendre.hpp>

using std::cout;
using std::endl;

// functions fit
void setPars(const double *par);
void scalePars(double *start_pars);
double globChi2(const double *par);
double normPhase(double a);
double phase(cd a1,cd a2) { return normPhase(arg(a1*conj(a2))); };
template <typename Type> Type getvalue(double M, vector< pair<double,Type> > &table);

//rho
double rhoRHO(double s);
double rhoF2(double s);

cd model_wave_template(int channel, int flag, double s, const double *pars);

//functions main
int fit_data(int nAttempt, const char *fout_name, int pars_all);

int plot_best_result(const char *fin_name, int entry);

int main(int ac, char **av) {

  if(ac<2) { cerr << "first arg is FIT or PLOT or SHEET" << endl; return 1; }

  const int nAttempt    = (ac > 2) ? atoi(av[2]) : 1;
  const char *file_name = (ac > 3) ? av[3] : "/tmp/test.root";
  const int pars_all    = (ac > 4) ? atoi(av[4]) : 2444;  

  if(strcmp(av[1],"FIT")==0) {
    return fit_data(nAttempt,file_name,pars_all);
  } else if(strcmp(av[1],"PLOT")==0) {
    return plot_best_result(file_name,nAttempt);
  } else if(strcmp(av[1],"BAND")==0) {
    return 0;//get_error_band(nAttempt,file_name);
  } else { cerr << "first arg is FIT or SHEET or PLOT" << endl; return 1; }
  
  return 0;
}

TMultiGraph *gm[6];
vector< pair<double,double> > glookup_rho_RHO, glookup_rho_F2;

double rhoF2lookup(double s) { return getvalue(sqrt(s),glookup_rho_F2); }
double rhoRHOlookup(double s) { return getvalue(sqrt(s),glookup_rho_RHO); }

cd prod_template(int nCh, double s, const double *pars);

MStructureHolder *gStrHolder;
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//double gNorm1,gNorm2,gNorm3;
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
template <typename Type> void adjust_plot(Type &obj);
TCanvas *plotAll(bool is_data = true);
int build_the_structure() {

  //**************** Load the data *****************//
  TGraph2D g1 ("data/f2piS.txt");  g1.SetName("g1"); 
  TGraph2D g2 ("data/f2piD.txt");  g2.SetName("g2");  
  TGraph2D g3 ("data/rhopiF.txt"); g3.SetName("g3");  
  TGraph2D g21("data/phi_f2piD_f2piS.txt");  g21.SetName("g21");
  TGraph2D g31("data/phi_rhopiF_f2piS.txt"); g31.SetName("g31");
  TGraph2D g32("data/phi_rhopiF_f2piD.txt"); g32.SetName("g32");
  const int nBins = 100;
  TH1D *h1 = new TH1D("h1","2^{-+}0^{+} f_{2}#pi S;M_{3#pi}",nBins,0.5,2.5);
  TH1D *h2 = new TH1D("h2","2^{-+}0^{+} f_{2}#pi D;M_{3#pi}",nBins,0.5,2.5);
  TH1D *h3 = new TH1D("h3","2^{-+}0^{+} #rho#pi F;M_{3#pi}",nBins,0.5,2.5);
  TH1D *h12 = new TH1D("phi_f2piS_f2piD","#Delta(f_{2}#pi S,f_{2}#pi D);M_{3#pi}",nBins,0.5,2.5);
  TH1D *h13 = new TH1D("phi_f2piS_rhopiF","#Delta(#rho#pi S,#rho#pi F);M_{3#pi}",nBins,0.5,2.5);
  TH1D *h23 = new TH1D("phi_f2piD_rhopiF","#Delta(#rho#pi D,#rho#pi F);M_{3#pi}",nBins,0.5,2.5);
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
    h12->SetBinContent(i,normPhase(-g21.GetY()[i-1]*M_PI/180));
    double err = g21.GetZ()[i-1]*M_PI/180;
    h12->SetBinError  (i, err > M_PI ? M_PI : err);
  }
  for(int i=1;i<=nBins;i++) {
    h13->SetBinContent(i,normPhase(-g31.GetY()[i-1]*M_PI/180));
    double err = g31.GetZ()[i-1]*M_PI/180;
    h13->SetBinError  (i, err > M_PI ? M_PI : err);
  }
  for(int i=1;i<=nBins;i++) {
    h23->SetBinContent(i,normPhase(-g32.GetY()[i-1]*M_PI/180));
    double err = g32.GetZ()[i-1]*M_PI/180;
    h23->SetBinError  (i, err > M_PI ? M_PI : err);
  }

  //phase
  //h12->Scale(M_PI/180.);
  //h13->Scale(M_PI/180.);
  //h23->Scale(M_PI/180.);

  //***************** look up tables ***************//
  double gMth = 3*PI_MASS;
  double gSth = gMth*gMth;

  const int Nrho = 1000;
  const double Mhg = 3.0;
  for(int i=0;i<Nrho;i++) {
    double mass = gMth + (Mhg - gMth)/(Nrho-1)*i;
    glookup_rho_RHO.push_back( make_pair( mass, rhoRHO(mass*mass) ) );
    glookup_rho_F2 .push_back( make_pair( mass, rhoF2 (mass*mass) ) );
  }

  //Model
  gStrHolder = new MStructureHolder();
  //intensities
  int defaultNumPar = 0;
  gStrHolder->AddWave(*h1,
		      [](int,double,const double*)->cd {return 0.0;},
		      defaultNumPar,glookup_rho_F2 ); 
  gStrHolder->AddWave(*h2,
		      [](int,double,const double*)->cd {return 0.0;},
		      defaultNumPar,glookup_rho_F2);  
  gStrHolder->AddWave(*h3,
		      [](int,double,const double*)->cd {return 0.0;},
		      defaultNumPar,glookup_rho_RHO);
  //
  gStrHolder->SetWaveRange(0, 1.4, 2.4);//2.4);
  gStrHolder->SetWaveRange(1, 1.6, 2.3 );//2.3);
  gStrHolder->SetWaveRange(2, 1.3, 2.1 );//2.3);

  //interference  
  gStrHolder->AddInterference(*h12,0,1,phase); gStrHolder->SetInterfRange(0, 1.6, 2.3);//2.4);
  gStrHolder->AddInterference(*h13,0,2,phase); gStrHolder->SetInterfRange(1, 1.4, 2.3);//2.4);
  gStrHolder->AddInterference(*h23,1,2,phase); gStrHolder->SetInterfRange(2, 1.6, 2.3);//2.4);

  gStrHolder->JustToPlot(0,1);

  //separate data
  //cout << " ************ " << endl;
  gm[0]  = split_data(h1,gStrHolder->GetWaveLowRange(0),gStrHolder->GetWaveUpRange(0));
  gm[1]  = split_data(h2,gStrHolder->GetWaveLowRange(1),gStrHolder->GetWaveUpRange(1));
  gm[2]  = split_data(h3,gStrHolder->GetWaveLowRange(2),gStrHolder->GetWaveUpRange(2));
  gm[3] = split_data(h12,gStrHolder->GetInterfLowRange(0),gStrHolder->GetInterfUpRange(0));
  gm[4] = split_data(h13,gStrHolder->GetInterfLowRange(1),gStrHolder->GetInterfUpRange(1));
  gm[5] = split_data(h23,gStrHolder->GetInterfLowRange(2),gStrHolder->GetInterfUpRange(2));
  //cout << " ------------ " << endl;

  return 0;
}

template <typename Type> void adjust_plot(Type &obj/*,const char *fulltitle*/) {
  obj.SetLabelSize(0.06,"xy");
  obj.GetXaxis()->SetTitleSize(0.07);
  obj.GetXaxis()->SetTitleOffset(-0.5);
  //obj.SetTitle(fulltitle);
}

MmatrixK *ampA;
int gAPAR[3];

//////////////////////////////////////////////////////////////////////////////////////////////////////
///     ///  /////      //      //////////////////////////////////////////////////////////////////////
///  // ///  /////  //  ////  ////////////////////////////////////////////////////////////////////////
///     ///  // //  //  ////  ////////////////////////////////////////////////////////////////////////
///  //////     //      ////  ////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int plot_best_result(const char *fin_name, int entry) {
  
  TFile *fin = TFile::Open(fin_name);           if(!fin) {cerr<<"No file!"<<endl;return 1;}
  TTree *t = (TTree*)gDirectory->Get("mins");   if(!t)   {cerr<<"No tree!"<<endl;return 1;}

  //Get parameters first
  TObjArray *br = t->GetListOfBranches(); br->ls();
  const int Nch = 3;
  std::vector<TString> pName;
  std::vector<TString> gName[Nch];
  std::vector<TString> aName[Nch];

  //number of poles
  int ipar=0;
  ipar=0; while(1) { TString svName = TString::Format("mass%d",ipar); if( br->FindObject(svName.Data()) ) pName.push_back(svName); else break; std::cout << ipar << ": found pole" << std::endl; ipar++;}
  //couplings
  for(int i=0;i<pName.size();i++) 
    for(int j=0;j<Nch;j++) {
      TString svName = TString::Format("gr%d_to%d",i,j); 
      if( br->FindObject(svName.Data()) ) gName[j].push_back(svName);
      else { std::cerr << "Error: couplings are not found in the tree. Something is wrong!" << std::endl; return 1; }
    }
  //production function
  for(int j=0;j<Nch;j++) { 
    ipar=0; while(1) { TString svName = TString::Format("alpha%d_pr%d",j,ipar); if( br->FindObject(svName.Data()) ) aName[j].push_back(svName); else break; std::cout << svName.Data() << ": found production coefficient" << std::endl; ipar++; } 
  }

  //set result
  const int gNPLS = pName.size();
  for(int j=0;j<Nch;j++) gAPAR[j] = aName[j].size()/2;
  //Find entries number with minimal chi2
  const int N = t->GetEntries();
  int mentry = entry;
  if(entry<0) {
    double chi2; t->SetBranchAddress("chi2",&chi2); t->GetEntry(0);
    double chi2min = chi2;
    int emin=0;
    for(int i=1;i<N;i++) {
      if(i%100==0) std::cout << i << ", ";
      t->GetEntry(i);
      if(chi2<chi2min) { chi2min=chi2; emin=i;}
    } 
    cout << "\n-------------> Best chi2 = " << chi2min << " ("<<emin<<")" << endl;
    mentry = emin;
  }
  cout << "-------------> entry " << mentry <<" will be get <-----------------" << endl;
  
  int NumPar = (Nch + 1)*gNPLS; for(int j=0;j<Nch;j++) NumPar+=2*gAPAR[j];
  double final_pars[NumPar];
  //masses
  for(int i=0;i<pName.size();i++) t->SetBranchAddress(pName[i].Data(),&final_pars[i]);
  //couplings
  for(int i=0;i<gNPLS;i++) 
    for(int j=0;j<Nch;j++) 
      t->SetBranchAddress(gName[j][i].Data(),&final_pars[gNPLS+i*Nch+j]);
  //productions
  int counter = 0;
  for(int j=0;j<Nch;j++)
    for(int i=0;i<2*gAPAR[j];i++) {
      t->SetBranchAddress(aName[j][i].Data(),&final_pars[gNPLS*(Nch+1)+counter]);
      counter++;
    }
  
  cout << "-------------> Creating models: " << Nch << " channels, " 
       << gNPLS << " poles, " 
       << "(" << gAPAR[0] << "," << gAPAR[1] << "," << gAPAR[2] << ")  <-------------" << endl;
  

  //****************create the model****************//
  build_the_structure();

  ampA = new MmatrixK(3,gNPLS,{2*gAPAR[0],2*gAPAR[1],2*gAPAR[2]});
  ampA->setPhSp(0, rhoF2lookup);
  ampA->setPhSp(1, rhoF2lookup);
  ampA->setPhSp(2,rhoRHOlookup);
  //production
  ampA->setProd(0,[](double s,const double*pars)->cd {return prod_template(0,s,pars);});
  ampA->setProd(1,[](double s,const double*pars)->cd {return prod_template(1,s,pars);});
  ampA->setProd(2,[](double s,const double*pars)->cd {return prod_template(2,s,pars);});

  ampA->Print();

  const int Nchannels = ampA->getNch();
  if( NumPar != ampA->getNpar() ) { std::cerr << "Something wrong with NumPar" << endl; return 1; }

  gStrHolder->SetWaveModel(0,[](int flag,double s,const double*pars)->cd {return model_wave_template(0,flag,s,pars);},NumPar);
  gStrHolder->SetWaveModel(1,[](int flag,double s,const double*pars)->cd {return model_wave_template(1,flag,s,pars);},NumPar);
  gStrHolder->SetWaveModel(2,[](int flag,double s,const double*pars)->cd {return model_wave_template(2,flag,s,pars);},NumPar);

  //fill the model
  t->GetEntry(mentry);


  //////////////////////////////////////////////////////////////////////////////////////////////
  ///////       /////////////        ///////////////       ///////////////      ////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////

  //fit result
  setPars(final_pars);
  TCanvas *can = plotAll();
  can->SaveAs("/tmp/original.pdf");

  //only poles
  double reduced_pars[NumPar]; memcpy(reduced_pars,final_pars,NumPar*sizeof(double));
  counter = 0;
  for(int j=0;j<Nch;j++)
    for(int i=0;i<2*gAPAR[j];i++) {
      reduced_pars[gNPLS*(Nch+1)+counter] = (i==0) ? 1.0 : 0.0;
      counter++;
    }
  setPars(reduced_pars);
  TCanvas *can_poles = plotAll(0);
  can_poles->SaveAs("/tmp/poles.pdf");
  //reduced_shifted
  for(int i=0;i<gNPLS;i++) reduced_pars[i] = 100;
  setPars(reduced_pars);
  TCanvas *can_poles2 = plotAll(0);
  can_poles2->SaveAs("/tmp/poles_shifted.pdf");

  //only production
  double shifted_pars[NumPar]; memcpy(shifted_pars,final_pars,NumPar*sizeof(double));
  for(int i=0;i<gNPLS;i++) shifted_pars[i] = 100;
  setPars(shifted_pars);
  TCanvas *can_prod = plotAll(0);
  can_prod->SaveAs("/tmp/shifted.pdf");

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
  
  const int gNPLS    = (pars_all % 10000 ) / 1000;
  gAPAR[0] = (pars_all % 1000  ) / 100;
  gAPAR[1] = (pars_all % 100   ) / 10;
  gAPAR[2] = (pars_all % 10    ) / 1;

  //****************create the model****************//
  build_the_structure();

  ampA = new MmatrixK(3,gNPLS,{2*gAPAR[0],2*gAPAR[1],2*gAPAR[2]});
  ampA->setPhSp(0, rhoF2lookup);
  ampA->setPhSp(1, rhoF2lookup);
  ampA->setPhSp(2,rhoRHOlookup);
  //production
  ampA->setProd(0,[](double s,const double*pars)->cd {return prod_template(0,s,pars);});
  ampA->setProd(1,[](double s,const double*pars)->cd {return prod_template(1,s,pars);});
  ampA->setProd(2,[](double s,const double*pars)->cd {return prod_template(2,s,pars);});

  ampA->Print();

  const int Nchannels = ampA->getNch();
  const int NumPar    = ampA->getNpar();

  gStrHolder->SetWaveModel(0,[](int flag,double s,const double*pars)->cd {return model_wave_template(0,flag,s,pars);},NumPar);
  gStrHolder->SetWaveModel(1,[](int flag,double s,const double*pars)->cd {return model_wave_template(1,flag,s,pars);},NumPar);
  gStrHolder->SetWaveModel(2,[](int flag,double s,const double*pars)->cd {return model_wave_template(2,flag,s,pars);},NumPar);

  //preparation for drawing
  TCanvas *can = new TCanvas("c1","plots",0,0,900,600);
  can->Divide(3,3,0.0001,0.0001);

  //**********************Seed**********************//
  gRandom->SetSeed(0);/*12314*/ cout << "The seed " << gRandom->GetSeed() << " is used." << endl;
  //***************** Fit the data *****************//
  
  // Set the NAMES for the variables
  string   var_name[NumPar];
  double start_pars[NumPar]; 
  double   up_limit[NumPar];
  double  low_limit[NumPar];
  //masses
  for(int i=0;i<gNPLS;i++) { 
    ostringstream oname; oname<<"mass"<<i; var_name[i] = oname.str(); 
    if(i==0) { start_pars[i]=1.605;  low_limit[i]=1.55;  up_limit[i]=1.73; }
    if(i==1) { start_pars[i]=1.885; low_limit[i]=1.8;   up_limit[i]=1.9;  }
    if(i==2) { start_pars[i]=2.1;  low_limit[i]=2.0;   up_limit[i]=2.2;  }
  }
  //couplings
  for(int i=0;i<gNPLS;i++) 
    for(int j=0;j<Nchannels;j++) { 
      ostringstream oname; oname<<"gr"<<i<<"_to"<<j; var_name[gNPLS+i*Nchannels+j] = oname.str(); 
      start_pars[gNPLS+i*Nchannels+j] = 1.;
       low_limit[gNPLS+i*Nchannels+j] = 0.;
        up_limit[gNPLS+i*Nchannels+j] = 3.;
    }

  //productions
  int counter = 0;
  for(int j=0;j<Nchannels;j++)
    for(int i=0;i<2*gAPAR[j];i++) {
      ostringstream oname; oname<<"alpha"<<j<<"_pr"<<i; var_name[gNPLS*(Nchannels+1)+counter] = oname.str(); 
      start_pars[gNPLS*(Nchannels+1)+counter] = (i==0) ? 60. : 0. ;
       low_limit[gNPLS*(Nchannels+1)+counter] = -60.;
        up_limit[gNPLS*(Nchannels+1)+counter] = 70.;
      counter++;
    }

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
    
    const double step = 0.01;
    for(int i=0;i<NumPar;i++) min->SetVariable( i,var_name[i].c_str(), start_pars[i], step );
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,(up_limit[i]-low_limit[i])/1e3);
    
    //step
    for(int i=0;i<NumPar;i++) min->SetVariableStepSize(i,(up_limit[i]-low_limit[i])/1e3);
    //starting point
    for(int i=0;i<NumPar;i++)  
      start_pars[i] = low_limit[i] + (up_limit[i]-low_limit[i])*gRandom->Rndm();
    
    min->SetVariableValues(start_pars);
    cout<<"Initual parameters:\n";for(int i=0;i<NumPar;i++) {cout<<var_name[i]<<" = "<<start_pars[i]<<";\n";}

    gStrHolder->JustToPlot(0,0,true);
    gStrHolder->JustToPlot(0,2,true);
    min->Minimize();
    gStrHolder->JustToPlot(0,0,false);
    gStrHolder->JustToPlot(0,2,false);
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
    for(int w=0;w<3;w++) {m[w].Add((TMultiGraph*)gm[w]->Clone()); m[w].SetTitle(gm[w]->GetTitle());}
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
    for(int i=3;i<6;i++) {m[i].Add((TMultiGraph*)gm[i]->Clone()); m[i].SetTitle(gm[i]->GetTitle());}
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


void setPars(const double *par) {
  gStrHolder->pars[0] = vector<double>(par,par+ampA->getNpar());
  gStrHolder->pars[1] = vector<double>(par,par+ampA->getNpar());
  gStrHolder->pars[2] = vector<double>(par,par+ampA->getNpar());
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
cd model_wave_template(int channel, int flag, double s, const double *pars) {  
  //cout<<"ch"<<channel<<", flag"<<flag<<endl;

  if(flag==0) return ampA->getA(channel,s,pars);

  const int Nall = ampA->getNpar();
  const int Nch = ampA->getNch();
  const int Np  = ampA->getNp();
  double reduced_pars[ampA->getNpar()]; memcpy(reduced_pars,pars,ampA->getNpar()*sizeof(double));
  for(int t=1;t<=Np;t++)
    if(flag==t && Np>=flag) { 
      for(int i=0;i<Np;i++) 
	for(int j=0;j<Nch;j++) 
	  if(i!=flag-1) reduced_pars[Np+i*Nch+j] = 0; 
      //cout << endl;
      //for(int u=0;u<ampA->getNpar();u++) cout << reduced_pars[u] << " ";
      cd value = ampA->getA(channel,s,reduced_pars);
      //cout<<"ch"<<channel<<", flag"<<flag<<", s="<<s<<","<<value<<" ___ " << endl;
      return value;
    }

  cerr << "Chech the code. Something is going wrong" << endl;  return 0;  
}

#define gOmega(a,b,c) (sqrt(c)-sqrt(a-b))/(sqrt(c)+sqrt(a-b))
#define gS0 0.0
#define gB  9.0

cd prod_template(int nCh, double s, const double *pars) {
  
  int nPar = gAPAR[nCh];
  cd ca[ nPar ];
  for(int i=0;i<nPar;i++) ca[i] = cd(pars[2*i],pars[2*i+1]);

  double omega = gOmega(s,gS0,gB);
  
  cd value = 0.0; 
  for(int i=0;i<nPar;i++) value += ca[i]*boost::math::legendre_p(i,omega);
  return value;

}


TCanvas *plotAll(bool is_data) {

  TCanvas *can = new TCanvas("c1","plots",0,0,900,600);
  can->Divide(3,3,0.0001,0.0001);
  
  TGraphErrors *gr[3][3], *gr_whole[3][3], *gphi[3], *gphi_whole[3];
  //intesities
  for(int w=0;w<3;w++) 
    for(int i=0;i<3;i++) {
      gr[w][i] = gStrHolder->GetWavePlot  (w,i, 0,0);
      if(i==0) gr[w][i]->SetLineColor(kRed);
      if(i==1) gr[w][i]->SetLineColor(kBlue);
      if(i==2) gr[w][i]->SetLineColor(kGreen);
    }
  for(int w=0;w<3;w++) 
    for(int i=0;i<3;i++) {
      gr_whole[w][i] = gStrHolder->GetWavePlot  (w,i, 0.5,2.5);
      gr_whole[w][i]->SetLineStyle(2);
      gr_whole[w][i]->SetLineWidth(1);
      if(i==0) gr_whole[w][i]->SetLineColor(kRed);
      if(i==1) gr_whole[w][i]->SetLineColor(kBlue);
      if(i==2) gr_whole[w][i]->SetLineColor(kGreen);
    }
  //interferences
  for(int w=0;w<3;w++)  {
    gphi[w] = gStrHolder->GetInterfPlot(w, 0,0); gphi[w]->SetLineColor(kRed);
    gphi_whole[w] = gStrHolder->GetInterfPlot(w, 0.5,2.5); gphi_whole[w]->SetLineStyle(2); gphi_whole[w]->SetLineColor(kRed);
  }
  //plot
  //intensities
  TMultiGraph *m[6];
  for(int w=0;w<3;w++) {m[w]=new TMultiGraph(); if(is_data) m[w]->Add((TMultiGraph*)gm[w]->Clone()); m[w]->SetTitle(gm[w]->GetTitle());}
  for(int w=0;w<3;w++) {
    for(int i=2;i>=0;i--) {
      m[w]->Add((TGraphErrors*)gr[w][i]->Clone(),"lz");
      m[w]->Add((TGraphErrors*)gr_whole[w][i]->Clone(),"lz");	
    }
    int Npad = 1+w+w*3;
    can->cd(Npad); m[w]->Draw("a");
    adjust_plot(*m[w]->GetHistogram());
  }
  //interferences
  for(int i=3;i<6;i++) {m[i]=new TMultiGraph(); if(is_data) m[i]->Add((TMultiGraph*)gm[i]->Clone()); m[i]->SetTitle(gm[i]->GetTitle());}
  int count = 0;
  for(int wi=0;wi<3;wi++) {
    for(int wj=0;wj<3;wj++) {
      if(wj<=wi) continue;
      m[3+count]->Add((TGraphErrors*)gphi[count]->Clone(),"lz"); 
      m[3+count]->Add((TGraphErrors*)gphi_whole[count]->Clone(),"lz");
      int Npad = wi*3+wj+1;
      can->cd(Npad); m[3+count]->Draw("a");
      adjust_plot(*m[3+count]->GetHistogram());
      count++;
    }
  }

  return can;
}
