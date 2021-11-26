#include <iostream>
#include <vector>
#include "draw_distribution_from_tree.C"

#define SYMM true

int plot_a11420_raw_data(uint t_slice) {
  TString t_slices_name[] = {
    "0.100000-0.112853",
    "0.112853-0.127471",
    "0.127471-0.144385",
    "0.144385-0.164401",
    "0.164401-0.188816",
    "0.188816-0.219907",
    "0.219907-0.262177",
    "0.262177-0.326380",
    "0.326380-0.448588",
    "0.448588-0.724294",
    "0.724294-1.000000"
  };

  TString projection_file = TString("~/cernbox/tmp/pwa_results/amplitudes/mfit_")
    + t_slices_name[t_slice] + TString(".root");
  const char *phase_space = "/mnt/data/compass/2008/Deck.88waves.linal.misha/phsp.lowt/direct_%d_large1e6.root";
  const char *cut = "(((sqrt(s3)>0.95)&&(sqrt(s3)<0.995))||((sqrt(s1)>0.95)&&(sqrt(s1)<0.995)))&&((sqrt(s3)<0.65||sqrt(s3)>0.95)&&(sqrt(s1)<0.65||sqrt(s1)>0.95))";

  std::vector<uint> what;  // (87);
  for (uint i = 1; i <= 88; i++) if (i + 1 != 16) what.push_back(i);

  TString name;
  name = "h87_cut";
  TH1D *h87 = new TH1D(name, name, 200, 1.3, 1.7);
  for (int bin = 41; bin <= 60; bin++) {
    draw_distribution_from_tree(TString::Format("sqrt(s)>>+%s", name.Data()),
                                projection_file, phase_space, bin, what, SYMM, cut);
  }

  name = "h16_cut";
  TH1D *h16 = new TH1D(name, name, 200, 1.3, 1.7);
  for (int bin = 41; bin <= 60; bin++) {
    draw_distribution_from_tree(TString::Format("sqrt(s)>>+%s", name.Data()),
                                projection_file, phase_space, bin, {16}, SYMM, cut);
  }

  what.push_back(16);
  name = "hall_cut";
  TH1D *hall = new TH1D(name, name, 200, 1.3, 1.7);
  for (int bin = 41; bin <= 60; bin++) {
    draw_distribution_from_tree(TString::Format("sqrt(s)>>+%s", name.Data()),
                                projection_file, phase_space, bin, what, SYMM, cut);
  }

  TFile *f = new TFile(TString("/tmp/raw_data_a11420_")
                       + t_slices_name[t_slice]
                       + TString(".root"), "RECREATE");
  h16->Write();
  h87->Write();
  hall->Write();
  f->Close();

  return 0;
}

