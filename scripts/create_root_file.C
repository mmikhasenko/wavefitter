// Copyright [2017] Misha Mikhasenko

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TString.h"

int create_root_file() {

  TFile fout("integrals.root", "RECREATE");
  for (uint i = 1; i < 88; i++) {
    TString fname = TString::Format("/mnt/data/compass/2008/integrals_stefan/h%d.n.out", i+1);
    std::ifstream fstr(fname.Data());
    if (!fstr.is_open()) {std::cerr << "Error with fstream!\n"; return 1;}
    std::string stitle, line;
    if (!std::getline(fstr, stitle)) {std::cerr << "Someting is wrong with file " << fname << "\n"; return 1;}
    TGraph gr(fname.Data());
    for (uint i = 0; i < static_cast<uint>(gr.GetN()); i++) gr.GetY()[i] /= gr.GetX()[i];
    for (uint i = 0; i < static_cast<uint>(gr.GetN()); i++) gr.GetY()[i] /= gr.GetY()[gr.GetN()-1];
    gr.SetName(TString::Format("g%d", i+1));
    gr.Write();
//    // std::vector<double>
//    while (std::getline(fstr, line)) {
//      std::cout << "line: \"" << line << "\"\n";
//      std::istringstream iss(line);
//      double x; double y;
//      iss >> x; iss >> y;
//      std::cout << "x = " << x << ", y = " << y << "\n";
//    }
//    double step = (X[N-1]-X[0])/(N-1);
//    std::cout << TString::Format("h%d", i+1) <<line.substr(2).c_str() << N << X[1]-step/2 << X[N-2]+step/2 << "\n";
//    TH1D *h = new TH1D(TString::Format("h%d", i+1), line.c_str(), N, X[0]-step/2, X[N-1]+step/2);
//     for (uint j = 0; j < N; j++) {
//       h->SetBinContent(j+1, Y[j]);
//     }
//     h->Write();
  }
  fout.Close();
  return 0;
}
