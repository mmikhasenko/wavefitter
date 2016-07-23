// Copyright [2016] Misha Mikhasenko
#include <iostream>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "MDeck.h"

int add_weight_to_tree(const char *fin_name, bool save_flag = false, const char* fout_name = "/tmp/updated_test.root");
int add_weight_to_tree(const char *fin_name, bool save_flag, const char* fout_name) {
  // open file and check
  TFile *f = TFile::Open(fin_name);
  if (!f) {std::cout << "Error: no file" << std::endl; return 0;}
  TTree *tin = 0; gDirectory->GetObject("events", tin);
  if (!tin) {std::cout << "Error: no tree" << std::endl; return 0;}

  TFile *fout = (!save_flag) ? 0 : TFile::Open(fout_name, "recreate");
  TTree *tout = (!save_flag) ? 0 : tin->CloneTree();

  TLorentzVector *beam_lv   = new TLorentzVector();
  TLorentzVector *pi_lv     = new TLorentzVector();
  TLorentzVector *iso_lv = new TLorentzVector();
  TLorentzVector *recl_lv = new TLorentzVector();

  tin->SetBranchAddress("beam", &beam_lv);
  tin->SetBranchAddress("pion", &pi_lv);
  tin->SetBranchAddress("isobar", &iso_lv);
  tin->SetBranchAddress("recoil", &recl_lv);
  double w0 = 1;
  double w1 = 1;
  TBranch *bpt0 = (!tout) ? 0 : tout->Branch("weight_ascoli_simplified", &w0);
  TBranch *bpt1 = (!tout) ? 0 : tout->Branch("weight_crossed_helicity", &w1);

  TH1D *his = new TH1D("invMassSquare", "Square of the invariant mass of system", 100, 0, 9.);
  TH1D *ht  = new TH1D("transfM", "t distribution", 100, -1, 0);

  double mpi = 0.139;

  // for out mode;
  
  const int Nentries = tin->GetEntries();
  for (int i = 0; i < Nentries; i++) {
    tin->GetEntry(i);

    // calculate important quantities
    TLorentzVector pi3_lv = *pi_lv+*iso_lv;
    TLorentzVector t_lv = *beam_lv - pi3_lv;
    TLorentzVector t_exch_lv = *beam_lv - *iso_lv;
    double t_exch = t_exch_lv.M2();
    double s_ppi = (*recl_lv+pi3_lv).M();
    
    { /* simplified Ascoli*/
      // calculate amplitude
      double amp = s_ppi/(mpi*mpi-t_exch);
      // calculate weight
      w0 = amp*amp;
    }
    { /* crossed hilicity */
      // calculate amplitude
      double z = beam_lv->Vect().Unit().Dot(iso_lv->Vect().Unit());
      double deck = MDeck::getDeck(beam_lv->M2(), t_lv.M2(), iso_lv->M2(), pi_lv->M2(), mpi*mpi,
                                   pi3_lv.M2(), z,
                                   1, 0, 1, 0,
                                   5);
      double amp = deck;
      // calculate weight
      w1 = amp*amp;
    }

    // fill tree branches
    if (bpt0) bpt0->Fill();
    if (bpt1) bpt1->Fill();

    // fil some hostograms
    his->Fill(pi3_lv.M2(), w1);
    ht->Fill(t_lv.M2(), w1);
  }

  TCanvas c1("c1");
  his->Draw("hist");
  c1.Print("/tmp/weghted_plots.pdf(", "pdf");
  // put some plots here without parentheses
  //
  // hist->Draw()
  // c1.Print("/tmp/weghted_plots.pdf", "pdf");
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
