// Copyright [2016] Misha Mikhasenko
#include <iostream>
#include <complex>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "constants.h"
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
  TBranch *bpt1 = (!tout) ? 0 : tout->Branch("weight_ascoli_simplified_no_bose", &w1);

  TH1D *his = new TH1D("invMassSquare", "Square of the invariant mass of system", 100, 0, 9.);
  TH1D *ht  = new TH1D("transfM", "t distribution;t(GeV^{2})", 100, -1, 0);
  TH1D *hz  = new TH1D("scatt_angle", "z distribution", 100, -1, 1);

  const int Nentries = tin->GetEntries();
  for (int i = 0; i < Nentries; i++) {
    if (i%1000000 == 0 && i != 0) std::cout << "Processing entry " << i << "\n";
    tin->GetEntry(i);

    TLorentzVector t_lv = *beam_lv - *reso_lv;
    TLorentzVector trgt_lv = *recl_lv - t_lv;

    //**************************************************************************
    //* First method, boost and call twice *************************************
    //**************************************************************************

    double stot = (*recl_lv+*reso_lv).M2();
    double t = (*beam_lv-*reso_lv).M2();

    MDeck::fromLabToGJ(*pi1_lv, *pi2_lv, *pi3_lv, *beam_lv, trgt_lv);

    // calculate amplitude
    cd amp_w0 = 0.;

    amp_w0 += MDeck::getAmplitude(pi1_lv->Px(), pi1_lv->Py(), pi1_lv->Pz(),  // pi1-
                                  pi2_lv->Px(), pi2_lv->Py(), pi2_lv->Pz(),  // pi+
                                  t, POW2(PI_MASS), stot,
                                  beam_lv->M2(), trgt_lv.M2(), recl_lv->M2(),
                                  pi1_lv->M2(), pi2_lv->M2(), pi3_lv->M2());
    cd amp_w1 = amp_w0;

    amp_w0 += MDeck::getAmplitude(pi3_lv->Px(), pi3_lv->Py(), pi3_lv->Pz(),  // pi3-
                                  pi2_lv->Px(), pi2_lv->Py(), pi2_lv->Pz(),  // pi+
                                  t, POW2(PI_MASS), stot,
                                  beam_lv->M2(), trgt_lv.M2(), recl_lv->M2(),
                                  pi1_lv->M2(), pi2_lv->M2(), pi3_lv->M2());

    amp_w0 /= sqrt(2.);
    //**************************************************************************
    //* Alternatively one can use direct call **********************************
    //**************************************************************************
    //* cd amp_w0 = MDeck::symmetriedFromLab(pi1_lv->Px(), pi1_lv->Py(), pi1_lv->Pz(), pi1_lv->M2(),
    //*                                      pi2_lv->Px(), pi2_lv->Py(), pi2_lv->Pz(), pi2_lv->M2(),
    //*                                      pi3_lv->Px(), pi3_lv->Py(), pi3_lv->Pz(), pi3_lv->M2(),
    //*                                      beam_lv->Px(), beam_lv->Py(), beam_lv->Pz(), beam_lv->M2(),
    //*                                      trgt_lv.Px(), trgt_lv.Py(), trgt_lv.Pz(), trgt_lv.M2());
    //* cd amp_w1 = MDeck::nonSymmetriedFromLab(pi1_lv->Px(), pi1_lv->Py(), pi1_lv->Pz(), pi1_lv->M2(),
    //*                                         pi2_lv->Px(), pi2_lv->Py(), pi2_lv->Pz(), pi2_lv->M2(),
    //*                                         pi3_lv->Px(), pi3_lv->Py(), pi3_lv->Pz(), pi3_lv->M2(),
    //*                                         beam_lv->Px(), beam_lv->Py(), beam_lv->Pz(), beam_lv->M2(),
    //*                                         trgt_lv.Px(), trgt_lv.Py(), trgt_lv.Pz(), trgt_lv.M2());
    //**************************************************************************

    // calculate weight
    w0 = norm(amp_w0);
    w1 = norm(amp_w1);

    // fill tree branches
    if (bpt0) bpt0->Fill();
    if (bpt1) bpt1->Fill();

    // fil some hostograms
    his->Fill(reso_lv->M2(), w0);
    ht->Fill(t_lv.M2(), w0);

    // temperary fix
    // if (i == 2) break;
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
