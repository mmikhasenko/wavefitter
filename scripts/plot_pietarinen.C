// [01.2017] Misha Mikhasenko

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <algorithm>

#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TMultiGraph.h"

#define MAX_N_PIETARINEN 10
#define MAX_N_CHANNEL 7

typedef std::complex<double> cd;

void plot_pietarinen(TTree *tin, int entry) {
  uint ind = 0, Nch = 0;
  double pars[2*MAX_N_PIETARINEN*MAX_N_CHANNEL];
  std::vector<std::pair<std::pair<int, std::string>,
                        std::pair<int, std::string> > > links[MAX_N_CHANNEL];

  TObjArray *bl = tin->GetListOfBranches();
    for (uint ch=0; ch < MAX_N_CHANNEL; ch++) {
      for (int i=0; i < MAX_N_PIETARINEN; i++) {
      // generate name
      std::string name_r = std::string(TString::Format("%dcr%d", i, ch).Data());
      std::string name_i = std::string(TString::Format("%dci%d", i, ch).Data());
      // check if name is present in ListOfBranches()
      if (!bl->FindObject(name_r.c_str())) continue;
      // create channel
      if (Nch <= ch) Nch++;
      links[ch].push_back(std::make_pair(std::make_pair(ind, name_r),
                                         std::make_pair(ind+1, name_i)));
      tin->SetBranchAddress(name_r.c_str(), &pars[ind++]);
      tin->SetBranchAddress(name_i.c_str(), &pars[ind++]);
    }
  }
  if (links[0].size() == 0) std::cout << "Warning<>: nothing is found, check branch names.\n";

  // slope and lhc
  double slope, lhc;
  tin->SetBranchAddress("slope", &slope);
  tin->SetBranchAddress("lhc", &lhc);

  tin->GetEntry(entry);
  for (uint ch = 0; ch < Nch; ch++) {
    std::cout << "Channel " << ch << ":\n";
    for (auto & p : links[ch]) {
      std::cout << p.first.second << " " << pars[p.first.first] << ", "
                << p.second.second << " " << pars[p.second.first] << "\n";
    }
  }
  std::cout << "Map adjustment: slope = " << slope << ", lhc = " << lhc << "\n";
  // define functions
  std::function<double(double)> omega = [slope, lhc](double s)->double{
    double r = (s-lhc)/slope;
    return (1.-sqrt(r))/(1.+sqrt(r));
  };
  std::vector<std::function<cd(double)> > piet(Nch);
  for (uint ch = 0; ch < Nch; ch++)
    piet[ch] = [&, ch](double s)->cd {
      cd value = 0.;
      for (uint ord = 0; ord < links[ch].size(); ord++) {
        value += cd(pars[links[ch][ord].first.first], pars[links[ch][ord].second.first]) *
          pow(omega(s), ord);
      }
      return value;
    };
  
  // plot
  uint Npoints = 100;
  std::pair<double, double> xrange(0.5, 2.5);
  std::vector<double> x(Npoints);
  for (uint i=0; i < Npoints; i++) x[i] = xrange.first + (xrange.second-xrange.first)/(Npoints-1)*i;

  std::vector<cd> y(Npoints);
  std::vector<double> yr(Npoints), yi(Npoints);  // real and imag
  std::vector<double> ya(Npoints), yp(Npoints);  // abs and phase
  // 0-function
  std::vector<TMultiGraph*> ms(Nch);
  std::vector<std::pair<TGraph*, TGraph*> > pg(Nch);
  for (uint ch = 0; ch < Nch; ch++) {
    std::transform(x.begin(), x.end(), y.begin(),  [&, ch](double m)->cd{return piet[ch](m*m);});
    std::transform(y.begin(), y.end(), yr.begin(), [&](cd v)->double{return real(v);});
    std::transform(y.begin(), y.end(), yi.begin(), [&](cd v)->double{return imag(v);});
    std::transform(y.begin(), y.end(), ya.begin(), [&](cd v)->double{return abs(v);});
    TGraph *gr = new TGraph(Npoints, x.data(), yr.data()); gr->SetLineColor(kBlack); gr->SetTitle("real"); gr->SetFillStyle(0);
    TGraph *gi = new TGraph(Npoints, x.data(), yi.data()); gr->SetLineColor(kRed);   gi->SetTitle("imag"); gi->SetFillStyle(0);
    ms[ch] = new TMultiGraph(); ms[ch]->SetTitle(TString::Format("Channel #%d", ch));
    ms[ch]->Add(gr); ms[ch]->Add(gi);
    pg[ch].first  = new TGraph(Npoints, x.data(), ya.data()); pg[ch].first->SetLineColor(kGreen+1); pg[ch].first->SetTitle("abs"); pg[ch].first->SetFillStyle(0);
    ms[ch]->Add(pg[ch].first);
    // abs phase
    std::transform(y.begin(), y.end(), yp.begin(), [&](cd v)->double{return arg(v);});
    pg[ch].second = new TGraph(Npoints, x.data(), yp.data()); pg[ch].second->SetLineColor(kRed); pg[ch].second->SetTitle("phase");
  }
  TCanvas *can = new TCanvas("can", "canva", 0, 0, 500*Nch, 500*2);
  can->Divide(Nch, 2);
  for (uint ch = 0; ch < Nch; ch++) {
    can->cd(ch+1);     ms[ch]->Draw("al");  can->cd(ch+1)->BuildLegend();
    can->cd(Nch+ch+1); pg[ch].second->Draw("al");
  }
}
