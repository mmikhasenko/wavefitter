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

TH1D *make_hsum(std::vector<TH1D*> vec);

TCanvas *plot_analytical_deck_components(const char *fin_name) {

        std::vector< std::pair<std::string, std::vector<uint>> > waves;
        { std::vector<uint> v = { 1, 15}; waves.push_back(std::make_pair("Incoherent sum;M_{3#pi}", v)); }
        { std::vector<uint> v = { 4,  6}; waves.push_back(std::make_pair("1^{?+}", v)); }
        { std::vector<uint> v = { 1,  3}; waves.push_back(std::make_pair("0^{?+}", v)); }
        { std::vector<uint> v = { 7,  9}; waves.push_back(std::make_pair("2^{?+}", v)); }
        { std::vector<uint> v = {10, 12}; waves.push_back(std::make_pair("3^{?+}", v)); }
        { std::vector<uint> v = {13, 15}; waves.push_back(std::make_pair("4^{?+}", v)); }

        TCanvas *c1 = new TCanvas("c1", "title", 1000, 1000);

        TFile *fin = new TFile(fin_name);
        THStack *hs = new THStack("hs", "Intensities for the coherent sum of projections;M_{3#pi}");
        for (auto & wi : waves) {
                std::vector<TH1D*> comp;
                for (uint w = wi.second[0]; w <= wi.second[1]; w++) {
                        TH2D *h2; gDirectory->GetObject(TString::Format("h2int%d", w), h2);
                        if (!h2) {
                                std::cerr << "Error: hist " << w << " is not found!\n";
                                return 0;
                        } else {
                                std::cout << "Success! " << h2->GetTitle() << "\n";
                        }
                        TH1D *h1 = h2->ProjectionX();
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
        for (auto & h : vec) std::cout << h->GetTitle() << " " << h->GetBinContent(50) << "\n";
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
