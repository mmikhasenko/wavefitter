// Copyright [2017] Misha Mikhasenko
// Discription:
//     part of linear algebra method chain to project to PWs
//     script to invert matrix

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "deflib.h"
#include "MatrixInverse.h"

void fill_thresholds(std::vector<double> *thr);

int main(int ac, char **av) {
  if (ac < 4) { std::cerr << "Usage: ./invert_matrix [PH.SP file] [PROJ file] [OUT file]\n"; return -1; }
  const char *name_phsp = av[1];  // "/mnt/data/compass/2008/Deck.88waves.linal.misha/waves.calculate_phase_space_symm_with_BW2.root";
  const char *name_proj = av[2];  // "/mnt/data/compass/2008/Deck.88waves.linal.misha/waves.calculate_Deck_integrals_symm_with_BW2.root"
  const char *fout_name = av[3];  // "/tmp/result.invert_matrix_symm_BW2.root"

  TFile *fincr = new TFile(name_phsp);
  if (!fincr) { std::cerr << "Error: no file!\n"; return 1; }
  const uint Nwaves = 88;
  std::vector<uint> chs0;

  // get histograms from the file
  TH1D *hdiag[Nwaves];
  TH1D *hintf[Nwaves][Nwaves];
  for (uint i = 0; i < Nwaves; i++)
    for (uint j = 0; j < Nwaves; j++) {
      TH1D **h = &((i == j) ? hdiag[i] : hintf[i][j]); *h = 0;
      fincr->cd();
      if (i == j) gDirectory->GetObject(TString::Format("h%d", i+1), *h);
      if (j > i) gDirectory->GetObject(TString::Format("h1%03d%03d", i+1, j+1), *h);
      if (i > j) gDirectory->GetObject(TString::Format("h2%03d%03d", j+1, i+1), *h);
      if (*h == 0) {
        std::cerr << "Error: no hist (i = " << i << ",j = " << j << ")\n";
        return 1;
      }
    }

  TFile *fproj = new TFile(name_proj);
  if (!fproj) { std::cerr << "Error: no Deck file!\n"; return 1; }
  TH1D *hreal[Nwaves], *himag[Nwaves];
  for (uint i = 0; i < Nwaves; i++) {
    fproj->cd();
    gDirectory->GetObject(TString::Format("r%d", i+1), hreal[i]);
    gDirectory->GetObject(TString::Format("i%d", i+1), himag[i]);
    if (hreal[i] == 0 || himag[i] == 0) {
      std::cerr << "Error: no proj-hist (i = " << i << ")\n";
      return 1;
    }
  }

  const uint Nbins = 100;
  // result
  TH1D *hresr[Nwaves], *hresi[Nwaves], *hint[Nwaves];
  for (uint i = 0; i < Nwaves; i++) {
    hresr[i] = new TH1D(TString::Format("hr%d", i+1),
                        TString::Format("Expansion coeff, real part, %s", hdiag[i]->GetTitle()),
                        100, 0.5, 2.5);
    hresi[i] = new TH1D(TString::Format("hi%d", i+1),
                        TString::Format("Expansion coeff, imag part, %s", hdiag[i]->GetTitle()),
                        100, 0.5, 2.5);
    hint[i] = new TH1D(TString::Format("h%d", i+1),
                       TString::Format("Intensity %s", hdiag[i]->GetTitle()),
                       100, 0.5, 2.5);
  }

  // fix the standatd set, start from 1, FLAT is not included
  for (uint i = 1; i < 88; i++) chs0.push_back(i);

  std::vector<double> thresholds(Nwaves, 0.0); fill_thresholds(&thresholds);

  // const uint iBin = 40;  // loop over bins
  for(uint iBin = 1; iBin <= Nbins; iBin++) {

    // loop over waves to pich those which have non-zero content
    std::vector<uint> chs;
    for (uint c = 0; c < chs0.size(); c++) {
      // std::cout << "chs0[c] = " << chs0[c] << ", " << hdiag[chs0[c]]->GetBinContent(iBin) << "; ";
      if (hdiag[chs0[c]]->GetBinCenter(iBin) > thresholds[chs0[c]]) {
        if (hdiag[chs0[c]]->GetBinContent(iBin) != 0.0) {chs.push_back(chs0[c]);
        } else { std::cout << "Interesting! wave " << hdiag[chs0[c]]->GetName() << " has zero!\n"; }
      } else {
        std::cout << "wave " << hdiag[chs0[c]]->GetName() << ", "
                  << hdiag[chs0[c]]->GetTitle() << " is rejected by threshold!\n";
      }
    }
    std::cout << "Bin " << iBin << ": basis is selected, " << chs.size() << " waves\n";

    // finally, fill the matrix
    b::matrix<cd> crosses(chs.size(), chs.size());
    for (uint i = 0; i < chs.size(); i++)
      for (uint j = i; j < chs.size(); j++) {
        // std::cout << "(" << chs[i] << ", " << chs[j] << ")\n";
        if (i == j) { crosses(i, j) = hdiag[chs[i]]->GetBinContent(iBin);
        } else {
          crosses(i, j) = cd(hintf[chs[i]][chs[j]]->GetBinContent(iBin),
                             hintf[chs[j]][chs[i]]->GetBinContent(iBin));
          crosses(j, i) = crosses(i, j);
        }
      }
    // finally, invert
    bool sing = false;
    b::matrix<cd> crosses_inv = gjinverse(crosses, sing);

    // vector of B
    b::vector<cd> vBint(chs.size());
    for (uint i = 0; i < chs.size(); i++)
      vBint(i) = cd(hreal[chs[i]]->GetBinContent(iBin), himag[chs[i]]->GetBinContent(iBin));

    b::vector<cd> vBres = prod(crosses_inv, vBint);

    for (uint i = 0; i < chs.size(); i++) {
      hresr[chs[i]]->SetBinContent(iBin, real(vBres(i)));
      hresi[chs[i]]->SetBinContent(iBin, imag(vBres(i)));
      hint[chs[i]]->SetBinContent(iBin, norm(vBres(i)) * hdiag[chs[i]]->GetBinContent(iBin) *
                                  1./POW2(4*M_PI) * 1./(8*M_PI));  // because of the normalization of integrals
    }
  }

  // Save the result
  TFile fout(fout_name, "RECREATE");
  for (uint i = 0; i < Nwaves; i++) { hresr[i]->Write(); hresi[i]->Write(); hint[i]->Write(); }
  fincr->cd();
  TH1D *hnorm; gDirectory->GetObject("hnorm", hnorm);
  if (hnorm) { fout.cd(); hnorm->Write(); }  // just propagate histogram to the next file

  std::cout << "File " << fout.GetName() << " has been created!\n";

  return 0;
}


void fill_thresholds(std::vector<double> *thr) {
  std::vector<std::pair<uint, double> > Dindex_threshold =
    {
      std::make_pair(4 , 1.2 ),
      std::make_pair(6 , 1.6 ),
      std::make_pair(12, 1.1 ),
      std::make_pair(13, 1.22),
      std::make_pair(16, 1.18),
      std::make_pair(17, 1.14),
      std::make_pair(23, 1.0 ),
      std::make_pair(24, 1.4 ),
      std::make_pair(25, 0.8 ),
      std::make_pair(27, 1.1 ),
      std::make_pair(38, 1.0 ),
      std::make_pair(39, 1.3 ),
      std::make_pair(42, 1.16),
      std::make_pair(44, 1.34),
      std::make_pair(49, 0.96),
      std::make_pair(50, 1.14),
      std::make_pair(51, 1.38),
      std::make_pair(52, 1.38),
      std::make_pair(54, 1.38),
      std::make_pair(55, 1.38),
      std::make_pair(59, 1.6 ),
      std::make_pair(60, 1.4 ),
      std::make_pair(64, 1.7 ),
      std::make_pair(68, 1.36),
      std::make_pair(71, 0.98),
      std::make_pair(86, 1.18),
      std::make_pair(87, 1.3 )
    };
  for (auto && p : Dindex_threshold) thr->at(p.first - 1)  = p.second;
}
