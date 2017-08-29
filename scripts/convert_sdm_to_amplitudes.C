#include <iostream>

#include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"
// #include ".h"

int convert_sdm_to_amplitudes(const char *sdm, const char *phase_space,
                              const char *fout_name) {
  TFile *fsdm = TFile::Open(sdm); if (!fsdm) return -3;
  TFile *fphsp = TFile::Open(phase_space); if (!fphsp) return -2;
  TFile *fout = TFile::Open(fout_name, "RECREATE");

  TH1F FF("fake_phase", "Fake phase", 100, 0.5, 2.5);
  const uint Nwaves = 88;
  for (uint w=1; w <= Nwaves; w++) {
    fsdm->cd();
    TH1F *hsdm = (TH1F*)gDirectory->Get(TString::Format("h%d",w));
    if (!hsdm) { std::cout << "Can not find h" << w << "\n"; return -1; }
    TH1F *hphi = (TH1F*)gDirectory->Get(TString::Format("h4002%03d",w));
    if (w < 3 || w > 81) hphi = &FF;
    if (!hphi) { std::cout << "Can not find " << TString::Format("h4002%03d",w) << "\n"; return -1; }
    fphsp->cd();
    TH1F *hphsp = (TH1F*)gDirectory->Get(TString::Format("h%d",w));
    if (!hphsp) return -1;
    // set up the output
    fout->cd();
    TH1D *hr = new TH1D(TString::Format("hr%d", w),
                        TString::Format("Expansion coeff, real part, %s", hphsp->GetTitle()),
                        100, 0.5, 2.5);
    TH1D *hi = new TH1D(TString::Format("hi%d", w),
                        TString::Format("Expansion coeff, imag part, %s", hphsp->GetTitle()),
                        100, 0.5, 2.5);
    TH1D *hint = new TH1D(TString::Format("h%d", w),
                        TString::Format("Intensity %s", hphsp->GetTitle()),
                        100, 0.5, 2.5);
    for (uint c=0; c<100; c++) {
      double M = 0.51+0.02*c;
      double intensity = hsdm->GetBinContent(c+1);
      double phsp = hphsp->GetBinContent(c+1);
      double phi = hphi->GetBinContent(c+1);
      if (phsp == 0 && intensity == 0) {std::cout << "Something is wrong!\n"; return -5;}
      // Dimas intensity = |A|^2 * M * phase space
      double amp_mod = sqrt(intensity / M /phsp * (4*M_PI)*(4*M_PI) * (8*M_PI) );
      hr->SetBinContent(c+1, amp_mod*cos(phi/180*M_PI));
      hi->SetBinContent(c+1, amp_mod*sin(phi/180*M_PI));
      hint->SetBinContent(c+1, intensity);  // the factor (4pi)^2(8pi) will be devided in later
    }
    fout->cd();
    hr->Write();
    hi->Write();
    hint->Write();
  }
  
  return 0;
}

