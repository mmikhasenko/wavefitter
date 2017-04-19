// Copyright [2016] Misha Mikhasenko
#include <iostream>
#include <complex>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "Math/SpecFuncMathMore.h"

#include "MDeck.h"
#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "waves.h"

typedef struct {
  uint J;
  uint M;
  uint S;
  uint L;
  double threshold;
  std::string title;
} wave;

wave make_wave(uint _J, uint _M, uint _S, uint _L, double _threshold,
               const std::string &_title = "") {
  wave n; n.J = _J; n.M = _M; n.S = _S; n.L = _L; n.title = _title; n.threshold = _threshold;
  return n;
}

int main(int argc, char *argv[]) {

  std::vector<wave> waves;

  TH1D htest("text", "phase space 2^{++}0^{+}#rho#pi D-wave", 100, 0.5, 2.5);
  waves.push_back(make_wave(2, 0, 2, 2, 1.0, "2^{++}0^{+}#rho#pi D-wave"));
  const double iw = 0;  // wave index waves[iw];

  for (uint e = 0; e < 100; e++) {
    TString fin_name = TString::Format("/mnt/data/compass/2008/phase_space_MC/%d.root", e);

    // open file and check
    TFile *f = TFile::Open(fin_name);
    if (!f) {std::cout << "Error: no file" << std::endl; return 0;}
    TTree *tin = 0; gDirectory->GetObject("events", tin);
    if (!tin) {std::cout << "Error: no tree" << std::endl; return 0;}

    TLorentzVector *beam_lv    = new TLorentzVector();
    TLorentzVector *pi1_lv     = new TLorentzVector();
    TLorentzVector *pi2_lv     = new TLorentzVector();
    TLorentzVector *pi3_lv     = new TLorentzVector();
    TLorentzVector *recl_lv    = new TLorentzVector();
    TLorentzVector *reso_lv     = new TLorentzVector();

    tin->SetBranchAddress("beam", &beam_lv);
    tin->SetBranchAddress("ipion1", &pi1_lv);
    tin->SetBranchAddress("ipion2", &pi2_lv);
    tin->SetBranchAddress("bpion",  &pi3_lv);
    tin->SetBranchAddress("3pi", &reso_lv);
    tin->SetBranchAddress("recoil", &recl_lv);

    double mpisq = PI_MASS*PI_MASS;

    // for out mode;
    MIsobar rho_iso(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
    MIsobar  f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.);

    double integral = 0;
    const int Nentries = tin->GetEntries();
    for (int i = 0; i < Nentries; i++) {
      if (i%1000000 == 0 && i != 0) std::cout << "Processing entry " << i << "\n";
      tin->GetEntry(i);

      // introduce quantity which does not depend on subchannels
      TLorentzVector t_lv = *beam_lv - *reso_lv;

      cd amp_w0 = 0.;
      cd amp_w1 = 0.;
      for (uint bose = 0; bose < 2; bose++) {
        if (bose == 1) {
          TLorentzVector *tmp = pi1_lv;
          pi1_lv = pi3_lv;
          pi3_lv = tmp;
        }

        TLorentzVector iso_lv = *pi1_lv+*pi2_lv;
        // calculate important quantities
        TLorentzVector t_exch_lv = *beam_lv - iso_lv;
        double t_exch = t_exch_lv.M2();
        double s_ppi = (*recl_lv+*pi3_lv).M2();

        // calculate z between beam and pion from isobar in the isobar rest frame
        double s1 = iso_lv.M2();
        double mbsq = beam_lv->M2();
        double eb_iso_rf = (s1+mbsq-t_exch)/(2.*sqrt(s1));
        double epi_iso_rf = (s1+pi1_lv->M2()-pi2_lv->M2())/(2.*sqrt(s1));
        double pb_iso_rf = sqrt(LAMBDA(s1, mbsq, t_exch)/(4*s1));
        double ppi_iso_rf = sqrt(LAMBDA(s1, pi1_lv->M2(), pi2_lv->M2())/(4*s1));
        double z_iso = (eb_iso_rf*epi_iso_rf - (*beam_lv)*(*pi1_lv))/(pb_iso_rf*ppi_iso_rf);

        // other variables
        double misq = iso_lv.M2();

        /*************************************/
        // loop over waves
        integral /*[w]*/ += 0.0;  // should be an anglular function.
        /*************************************/
      }
    }
    
    f->Close();

    // waves[iw];
    htest.SetBinContent(e+1, integral);
  }

  TCanvas c1("c1");
  htest.Draw("hist");
  c1.Print("/tmp/weghted_plots.pdf", "pdf");



  return 0;
}
