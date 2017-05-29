// Copyright [2017] Misha Mikhasenko
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TCanvas.h"
#include "THStack.h"

typedef std::complex<double> cd;

#define NAME_PATT "h2int%d"

void find_numbers_of_hists_which_much_the_pattern(const char *name_patt,
                                                  std::vector<uint> from,
                                                  const char* match_pattern,
                                                  std::vector<uint> *result);

TH1D *make_hsum(std::vector<TH1D*> vec);
TH1D *integrate_ysq(TH2D *h);

TCanvas *plot_analytical_deck_components(const char *fin_name) {

        TFile *fin = new TFile(fin_name);

        std::vector< std::pair<std::string, std::vector<uint>> > waves;
        // all
        { std::vector<uint> v(213);
          for (uint i = 0; i < v.size(); i++) v[i] = i+1;
          waves.push_back(std::make_pair("The sum, J < 5;M_{3#pi}", v));}
        // positive reflectivity
        std::vector<uint> vpos;
        find_numbers_of_hists_which_much_the_pattern(NAME_PATT, waves[0].second, "+(S", &vpos);
        std::cout << "Found " << vpos.size() << " hists with positive refl.\n";
        // 1++
        { std::vector<uint> v = {};
          find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "1++", &v);
          waves.push_back(std::make_pair("1^{++}", v));}
        // 0-+
        { std::vector<uint> v = {};
          find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "0-+", &v);
          waves.push_back(std::make_pair("0^{-+}", v));}
        // 2-+
        { std::vector<uint> v = {};
          find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "2-+", &v);
          waves.push_back(std::make_pair("2^{-+}", v));}
        // 1-+
        { std::vector<uint> v = {};
          find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "1-+", &v);
          waves.push_back(std::make_pair("1^{-+}", v));}
        // 2++
        { std::vector<uint> v = {};
          find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "2++", &v);
          waves.push_back(std::make_pair("2^{++}", v));}
        // // neg reflectivity
        // std::vector<uint> vneg;
        // find_numbers_of_hists_which_much_the_pattern(NAME_PATT, waves[0].second, "-(S", &vneg);
        // std::cout << "Found " << vneg.size() << " hists with negative refl.\n";
        // waves.push_back(std::make_pair("#epsilon = (-)", vneg));
        // // 3++
        // { std::vector<uint> v = {};
        //   find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "3++", &v);
        //   waves.push_back(std::make_pair("3^{++}", v));}
        // // 4-+
        // { std::vector<uint> v = {};
        //   find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "4-+", &v);
        //   waves.push_back(std::make_pair("4^{-+}", v));}
        // // 3-+
        // { std::vector<uint> v = {};
        //   find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "3-+", &v);
        //   waves.push_back(std::make_pair("3^{-+}", v));}
        // // 4++
        // { std::vector<uint> v = {};
        //   find_numbers_of_hists_which_much_the_pattern(NAME_PATT, vpos, "4++", &v);
        //   waves.push_back(std::make_pair("4^{++}", v));}


        TCanvas *c1 = new TCanvas("c1", "title", 1000, 1000);

        THStack *hs = new THStack("hs", "Intensities for the coherent sum of projections;M_{3#pi}");
        for (auto & wi : waves) {
                std::vector<TH1D*> comp;
                for (uint & w : wi.second) {
                        TH2D *h2; gDirectory->GetObject(TString::Format(NAME_PATT, w), h2);
                        if (!h2) {
                                std::cerr << "Error: hist " << w << " is not found!\n";
                                return 0;
                        }
                        TH1D *h1 = integrate_ysq(h2);
                        comp.push_back(h1);
                }
                TH1D *hsum = make_hsum(comp);
                for (auto & hi : comp) delete hi;
                hsum->SetTitle(wi.first.c_str());
                hsum->SetStats(kFALSE);
                hs->Add(hsum);
        }
        hs->Draw("pfc nostack");

        c1->BuildLegend();
        return c1;
}

TH1D *make_hsum(std::vector<TH1D*> vec) {
        // for (auto & h : vec) std::cout << h->GetTitle() << " " << h->GetBinContent(50) << "\n";
        const uint Nbins = 100;
        TH1D *hsum = new TH1D("hsum", "title", Nbins, 0.5, 2.5);
        for (uint b = 0; b < Nbins; b++) {
                hsum->SetBinContent(b+1, 0.0);
                for (auto & h : vec) {
                        hsum->SetBinContent(b+1, hsum->GetBinContent(b+1)+h->GetBinContent(b+1));
                }
        }
        return hsum;
}

TH1D *integrate_ysq(TH2D *h) {
        double widthX = h->GetXaxis()->GetBinWidth(1);
        double widthY = h->GetYaxis()->GetBinWidth(1);
        uint NbX = h->GetXaxis()->GetNbins();
        uint NbY = h->GetYaxis()->GetNbins();
        TH1D *hr = new TH1D(TString::Format("%s_int_ysq", h->GetName()),
                            TString::Format("%s_int_ysq", h->GetTitle()),
                            NbX,
                            h->GetXaxis()->GetBinLowEdge(1),
                            h->GetXaxis()->GetBinLowEdge(NbX)+widthX);
        for (uint i = 1; i <= NbX; i++) {
                double v = 0;
                for (uint j = 1; j <= NbY; j++) {
                        v += h->GetBinContent(i, j) *
                             2 * h->GetYaxis()->GetBinCenter(j) * // jacobian
                             widthY;
                }
                hr->SetBinContent(i, v);
        }
        return hr;
}

void find_numbers_of_hists_which_much_the_pattern(const char *name_patt,
                                                  std::vector<uint> from, const char* match_pattern,
                                                  std::vector<uint> *result) {
        for (uint & i : from) {
                TH2D *h; gDirectory->GetObject(TString::Format(name_patt, i), h);
                if (!h) {
                        std::cerr << "Error: hist with the name " << TString::Format(name_patt, i) << " not found\n";
                        return;
                }
                std::string title(h->GetTitle());
                if (title.find(match_pattern) != std::string::npos) result->push_back(i);
        }
}
