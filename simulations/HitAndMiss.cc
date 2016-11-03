// Copyright [2016] Misha Mikhasenko

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"



int main(int ac, char *av[]) {

  const char *fname_weighted = (ac >= 2) ? av[1] : "/tmp/tree_test.root";
  const char *fname_generated = (ac >= 3) ? av[2] : "/tmp/tree_test_generated.root";
  std::cout << " Input: " << fname_generated << "\n";
  std::cout << "Output: " << fname_weighted << "\n";
  
  TFile *fold = TFile::Open(fname_weighted);
  if (!fold) return 1;
  TTree *told; fold->GetObject("events", told);
  if (!told) return 1;

  // find max weight
  double max_weight = 0;
  double Nentries = told->GetEntries();
  
  double weight_ascoli_simplified;
  told->SetBranchAddress("weight_ascoli_simplified", &weight_ascoli_simplified);
  for (uint i = 0; i < Nentries; i++) {
    told->GetEntry(i);
    if (weight_ascoli_simplified > max_weight)
      max_weight = weight_ascoli_simplified;
  }

  // Hit and Miss
  TFile *fnew = TFile::Open(fname_generated, "recreate");
  TTree *tnew = told->CloneTree(0);
  TRandom3 r3(12345);
  uint Ngen = 0;
  for (uint i = 0; i < Nentries; i++) {
    told->GetEntry(i);
    if (r3.Rndm()*max_weight < weight_ascoli_simplified) {
      tnew->Fill();
      Ngen++;
    }
  }

  tnew->AutoSave();
  fnew->Close();
  fold->Close();

  std::cout << "Events are generated(" << Ngen << "), efficiency is " << 1.0*Ngen/Nentries << "\n";

  return 0;
}
