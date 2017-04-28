// Copyright [2016] Misha Mikhasenko
// Description:
//   The program calculates integrals numerically based on pregenerated MC data sample

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "MDeck.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "waves.h"
#include "mintegrate.h"
#include "M3bodyAngularBasis.h"
#include "TText.h"

typedef struct {
  uint index;
  uint J;
  bool parity;
  uint M;
  bool pos_refl;
  int S;
  uint L;
  double threshold;
  std::string title;
} wave;

void fill_wavepull(const char* wave_fname, std::vector<wave> *waves);

int main(int ac, char *av[]) {

  std::vector<wave> waves;
  fill_wavepull("/localhome/mikhasenko/results/pwa_3pi/wavelist_formated.txt", &waves);
  uint nWaves = waves.size();
  std::cout << "waves.size() = " << waves.size() << "\n";
  for (auto & w : waves)
    std::cout << w.index << ": " << w.title << " "
              << w.J << " " << (w.parity ? "+" : "-") << " " << w.M << " " << (w.pos_refl ? "+" : "-")
              << " " << w.S << " " << w.L << "\n";


  // for out mode;
  MIsobar rho(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.); rho.setIntU();
  MIsobar  f2(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.); f2.setIntU();
  MIsobar rho3(1.69, 0.16, PI_MASS, PI_MASS, 1, 5.); rho3.setIntU();
  MIsobarPiPiS pipiS; pipiS.setIntU();
  MIsobar *iso[] = {&pipiS, &rho, &f2, &rho3};
  MIsobar f980(0.99, 0.04, PI_MASS, PI_MASS, 0); f980.setIntU();
  MIsobar f1500(1.504, 0.11, PI_MASS, PI_MASS, 0); f1500.setIntU();
  MIsobar *iso_scalars[] = {&pipiS, &f980, &f1500};

  for (uint e = 0; e < 100; e++) {
    std::cout << "---> File #" << e << "\n";
    TString fin_name = TString::Format("/mnt/data/compass/2008/phase_space_MC/_with_deck_%d.root", e);  // _large1e6

    // open file and check
    TFile *f = TFile::Open(fin_name);
    if (!f) {std::cout << "Error: no file" << std::endl; return 0;}
    TTree *tin = 0; gDirectory->GetObject("angles", tin);
    if (!tin) {std::cout << "Error: no tree" << std::endl; return 0;}

    // out file, recreate
    TString fout_name = TString::Format("/mnt/data/compass/2008/phase_space_MC/_with_deck_and_PWs_%d.root", e);  // _large1e6
    TFile *fout = TFile::Open(fout_name, "recreate");
    TTree *tout = tin->CloneTree();
    // add a branch for every wave
    double amp_real[nWaves][2], amp_imag[nWaves][2];
    TBranch *br[nWaves][4];
    for (uint w = 0; w < nWaves; w++) {
      br[w][0] = tout->Branch(TString::Format("amp%d_frame1_real", waves[w].index), &amp_real[w][0]);
      br[w][1] = tout->Branch(TString::Format("amp%d_frame1_imag", waves[w].index), &amp_imag[w][0]);
      br[w][2] = tout->Branch(TString::Format("amp%d_frame3_real", waves[w].index), &amp_real[w][1]);
      br[w][3] = tout->Branch(TString::Format("amp%d_frame3_imag", waves[w].index), &amp_imag[w][1]);
    }
    
    double s;
    tin->SetBranchAddress("s", &s);

    // (23)-frame
    double s1, costheta1, phi1, costheta23, phi23;
    tin->SetBranchAddress("s1", &s1);
    tin->SetBranchAddress("costheta1", &costheta1);
    tin->SetBranchAddress("phi1", &phi1);
    tin->SetBranchAddress("costheta23", &costheta23);
    tin->SetBranchAddress("phi23", &phi23);

    // (12)-frame
    double s3, costheta3, phi3, costheta12, phi12;
    tin->SetBranchAddress("s3", &s3);
    tin->SetBranchAddress("costheta3", &costheta3);
    tin->SetBranchAddress("phi3", &phi3);
    tin->SetBranchAddress("costheta12", &costheta12);
    tin->SetBranchAddress("phi12", &phi12);

    // create and clean integral variables
    cd integrals[nWaves][nWaves];
    for (uint iw = 0; iw < waves.size(); iw++)
      for (uint jw = 0; jw < waves.size(); jw++)
        integrals[iw][jw] = 0.;

    // integration loop
    cd amp[nWaves][2];
    const int Nentries = tin->GetEntries();
    for (int i = 0; i < Nentries; i++) {
      if (i%1000000 == 0 && i != 0) std::cout << "Processing entry " << i << "\n";
      tin->GetEntry(i);

      for (uint bose = 0; bose < 2; bose++) {
        /*************************************/
        double sI        = bose?         s3 :         s1;
        double costhetaI = bose?  costheta3 :  costheta1;
        double phiI      = bose?       phi3 :       phi1;
        double costheta  = bose? costheta12 : costheta23;
        double phi       = bose?      phi12 :      phi23;

        double thetaI = acos(costhetaI);
        double theta  = acos(costheta);
        // loop over waves
        for (uint w = 0; w < nWaves; w++) {
          cd iso_shape = waves[w].S == -7 ? 1. :
            (waves[w].S > 0 ?
             iso        [ waves[w].S]->ToneVertex(sI) :
             iso_scalars[-waves[w].S]->ToneVertex(sI));

          double R = 5.;
          double qsq = LAMBDA(s, sI, POW2(PI_MASS))/(4*s);
          double BlattWeisskopf = pow(R*R*qsq/(1.+R*R*qsq), waves[w].L/2.);

          amp[w][bose] = Math::ZJMLS_refl(waves[w].J, waves[w].M,
                                          (waves[w].pos_refl == waves[w].parity),  // (-1)*(-1) = (+1)*(+1) = true, otherwise is false
                                          waves[w].L, (waves[w].S > 0 ? waves[w].S : 0),
                                          thetaI, phiI, theta, phi) * iso_shape * BlattWeisskopf;
          amp_real[w][bose] = real(amp[w][bose]);
          amp_imag[w][bose] = imag(amp[w][bose]);
        }
        /*************************************/
      }
      // tout->Fill();
      for (uint w = 0; w < nWaves; w++) for (uint b = 0; b < 4; b++) br[w][b]->Fill();
    }
    f->Close();
    fout->cd(); tout->Write();
    // save titles
    for (uint w = 0; w < nWaves; w++) {
      TText tx(0., 0., waves[w].title.c_str()); tx.SetName(TString::Format("t%d", waves[w].index));
      tx.Write();
    }
    std::cout << "File " << fout->GetName() << " have been completed!\n";
    fout->Close();
  }

  std::cout << "Done!\n";
  return 0;
}


void fill_wavepull(const char* wave_fname, std::vector<wave> *waves) {
  std::ifstream fin(wave_fname);
  if (!fin.is_open()) {
    std::cerr << "Error: can not find the file!\n";
    return;
  }
  std::string line;
  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    wave n;
    iss >> n.index;
    iss >> n.title;
    if (n.title == "FLAT") {
      n.J = 0; n.M = 0; n.S = -7; n.L = 0;
      n.parity = false;
      n.pos_refl = true;
      waves->push_back(n);
      continue;
    }
    iss >> n.J;
    std::string parity;
    iss >> parity;
    n.parity = (parity == "+");
    iss >> n.M;
    std::string epsilon;
    iss >> epsilon;
    n.pos_refl = (epsilon == "+");
    iss >> n.S;
    iss >> n.L;
    waves->push_back(n);
  }  
}

/***************************************************************************************/
/**** script to plot matrix ************************************************************/
/***************************************************************************************/
/*
TCanvas c1
c1.Divide(4,4)
for(int i=0; i<4; i++) {
for (int j=0; j<4; j++) {
c1.cd(1+i*4+j);
if(i==j) gDirectory->Get(TString::Format("h%d",i+1))->Draw();
if(i<j) gDirectory->Get(TString::Format("h1%03d%03d",i+1,j+1))->Draw();
if(i>j) gDirectory->Get(TString::Format("h2%03d%03d",j+1,i+1))->Draw();
}
}
*/
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/


/***************************************************************************************/
/**** script to plot in the matrix form and colors *************************************/
/***************************************************************************************/
/*
TMatrixD m(88,88);
for(int i=0; i<88;i++) {
  for (int j=0; j<88; j++) {
    TH1D *h = (TH1D*)gDirectory->Get((i==j) ? TString::Format("h%d",i+1) : ((i<j) ? TString::Format("h1%03d%03d",i+1,j+1) : TString::Format("h1%03d%03d",j+1,i+1)));
    m(i,j) = h->GetBinContent(81);
  }
}
TCanvas *c1 = new TCanvas("c1", "title", 1000, 1000);
m.Draw("colz");
TH2D *h = (TH2D*)gROOT->FindObject("TMatrixDBase");
if (h) {
  for(int i=0; i<88;i++) {
    TH1D *hdiag = (TH1D*)gROOT->FindObject(TString::Format("h%d", i+1));
    if (!hdiag) break;
    h->GetXaxis()->SetBinLabel(i+1, hdiag->GetTitle());
    h->GetYaxis()->SetBinLabel(i+1, hdiag->GetTitle());
  }
}
h->SetStats(kFALSE);
h->GetXaxis()->SetLabelSize(0.02);
h->GetYaxis()->SetLabelSize(0.02);
h->GetZaxis()->SetRangeUser(0.01, 1.0);
h->GetXaxis()->SetTickSize(0);
h->GetYaxis()->SetTickSize(0);
h->SetTitle("PW integral matrix");
h->Draw("colz");
*/
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
