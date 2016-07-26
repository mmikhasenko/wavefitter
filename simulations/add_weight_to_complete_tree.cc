// Copyright [2016] Misha Mikhasenko
#include <iostream>
#include <complex>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "MDeck.h"
#include "MIsobar.h"

int add_weight_to_tree(const char *fin_name, bool save_flag = false, const char* fout_name = "/tmp/updated_test.root");
int add_weight_to_tree(const char *fin_name, bool save_flag, const char* fout_name) {
  // open file and check
  TFile *f = TFile::Open(fin_name);
  if (!f) {std::cout << "Error: no file" << std::endl; return 0;}
  TTree *tin = 0; gDirectory->GetObject("events", tin);
  if (!tin) {std::cout << "Error: no tree" << std::endl; return 0;}

  TFile *fout = (!save_flag) ? 0 : TFile::Open(fout_name, "recreate");
  TTree *tout = (!save_flag) ? 0 : tin->CloneTree();

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
  double w0 = 1;
  double w1 = 1;
  TBranch *bpt0 = (!tout) ? 0 : tout->Branch("weight_ascoli_simplified", &w0);
  TBranch *bpt1 = (!tout) ? 0 : tout->Branch("weight_crossed_helicity", &w1);

  TH1D *his = new TH1D("invMassSquare", "Square of the invariant mass of system", 100, 0, 9.);
  TH1D *ht  = new TH1D("transfM", "t distribution;t(GeV^{2})", 100, -1, 0);
  TH1D *hz  = new TH1D("scatt_angle", "z distribution", 100, -1, 1);

  double mpisq = PI_MASS*PI_MASS;

  // for out mode;
  MIsobar rho_iso(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
  MIsobar  f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.);
  MIsobar  pipiS_iso(0.5, 0.5,  PI_MASS, PI_MASS, 0, 5.);
  
  const int Nentries = tin->GetEntries();
  for (int i = 0; i < Nentries; i++) {
    if (i%1000000 == 0 && i!=0) std::cout << "Processing entry " << i << "\n";
    tin->GetEntry(i);

    // introduce quantity which does not depend on subchannels
    TLorentzVector t_lv = *beam_lv - *reso_lv;


    cd amp_w0 = 0.;
    cd amp_w1 = 0.;
    for (uint bose = 0; bose < 2; bose++ ) {
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
    
      /************* simplified Ascoli **************/
      // calculate amplitude
      double exch_amp_w0 = s_ppi/(mpisq-t_exch);
      cd pipi_amp_w0 = 0.;
      MIsobar *iso;
      iso = &rho_iso;   pipi_amp_w0 += (2*iso->GetL()+1)*ROOT::Math::legendre(iso->GetL(), z_iso) * iso->T(misq);
      iso = &f2_iso;    pipi_amp_w0 += (2*iso->GetL()+1)*ROOT::Math::legendre(iso->GetL(), z_iso) * iso->T(misq);
      iso = &pipiS_iso; pipi_amp_w0 += (2*iso->GetL()+1)*ROOT::Math::legendre(iso->GetL(), z_iso) * iso->T(misq);
      // calculate weight
      amp_w0 += exch_amp_w0*pipi_amp_w0;
    
      /************* crossed helicity ***************/
      // calculate amplitude
      double mt1sq = t_lv.M2();
      double m4sq = pi3_lv->M2();

      double s3pi = reso_lv->M2();
      // calculate z between beam and isobar in 3pi rest frame
      double eb_3pi_rf = (s3pi+mbsq-mt1sq)/(2*sqrt(s3pi));
      double ei_3pi_rf = (s3pi+misq-m4sq)/(2*sqrt(s3pi));
      double pb_3pi_rf = sqrt(LAMBDA(s3pi, mbsq, mt1sq)/(4*s3pi));
      double pi_3pi_rf = sqrt(LAMBDA(s3pi, misq, m4sq)/(4*s3pi));
      double z = (eb_3pi_rf*ei_3pi_rf - (*beam_lv)*(iso_lv)) / (pb_3pi_rf*pi_3pi_rf);
      if (fabs(z)>1.) { std::cout << "Warning: |z|>1. Sc.pr" << (*beam_lv)*(iso_lv)
                                  << ", ee = " << eb_3pi_rf*ei_3pi_rf
                                  << ", pp = " << pb_3pi_rf*pi_3pi_rf
                                  << "\n"; continue; } 
      hz->Fill(z);
      double deck = MDeck::getDeck(mbsq, mt1sq, misq, m4sq, mpisq,
                                   s3pi, z,
                                   1, 0, 1, 0,
                                   5);
      
      amp_w1 += deck;
    }

    // calculate weight
    w0 = norm(amp_w0);
    w1 = norm(amp_w1);
  
    // fill tree branches
    if (bpt0) bpt0->Fill();
    if (bpt1) bpt1->Fill();

    // fil some hostograms
    his->Fill(reso_lv->M2(), w1);
    ht->Fill(t_lv.M2(), w1);
  }

  TCanvas c1("c1");
  his->Draw("hist");
  c1.Print("/tmp/weghted_plots.pdf(", "pdf");
  // put some plots here without parentheses
  //
  hz->Draw("hist");
  c1.Print("/tmp/weghted_plots.pdf", "pdf");
  //
  ht->Draw("hist");
  c1.Print("/tmp/weghted_plots.pdf)", "pdf");

  if (fout && tout) {
    fout->cd();
    tout->Write();
    fout->Close();
  }
  f->Close();

  return 0;
}

int main(int argc, char *argv[]) {
  if (argc == 2) { return add_weight_to_tree(argv[1]);
  } else if (argc == 3) { return add_weight_to_tree(argv[1], true, argv[2]);
  } else { std::cerr << "Usage: ./add_weight_to_tree fin_name fout_name\n"; }
  return 0;
}
