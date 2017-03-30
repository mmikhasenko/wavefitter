#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "TCanvas.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TList.h"
#include "TGraph.h"

TMultiGraph * applyShift(TMultiGraph *, int v);
double toPi(double v);

#define C180(v) (v/M_PI*180.)

int plotSpinDensity_v2(const char *fin_name, uint ind) {

  // uint ind = 1;
  std::string tranges[] = {"0.100", "0.113", "0.127", "0.144",
                           "0.164", "0.189", "0.220", "0.262",
                           "0.326", "0.449", "0.724", "1.000"};

  const uint Nsize = 4;
  std::string titles[Nsize] = {"2^{-+}0^{+} f_{2} #pi S",
                          "2^{-+}0^{+} f_{2} #pi D",
                          "2^{-+}0^{+} #rho #pi F",
                          "2^{-+}0^{+} (#pi#pi)_{S} #pi D"};

  int shift[Nsize][Nsize] = {{0, 0, 0, 0},
                     {0, 0, 0, 0},
                     {0, 0, 0, 0},
                     {0, 0, 0, 0}};
  
  TCanvas *c1 = new TCanvas("canva", "title", 0, 0, 1500, 1000);
  const double gmargin = 0.08;
  const double pmargin = 0.0;
  const double Xsize = (1. - 2*gmargin) / Nsize;
  const double Ysize = (1. - 2*gmargin) / Nsize;  

  TFile fin(fin_name);
  TList *ll = ((TCanvas*)fin.Get(TString::Format("page%d",ind)))->GetListOfPrimitives();
  TVirtualPad *p = (TVirtualPad*)ll->First();
  for (uint i = 0; i < Nsize; i++) {
    for (uint j = 0; j < Nsize; j++) {
      if (i != 0 || j != 0) p = (TVirtualPad*)ll->After(p);      
      if (i > j) continue;
      // create pad
      double xl = gmargin + Xsize*j;
      double xr = gmargin + Xsize*(j+1);
      double yl = 1. - (gmargin + Ysize*i);
      double yr = 1. - (gmargin + Ysize*(i+1));
      TPad *np = new TPad(TString::Format("pad_%02d",Nsize*i+j+1), "noTitle", xl, yl, xr, yr, 0); // i+j+1
      np->SetRightMargin(pmargin);
      np->SetTopMargin(pmargin);
      np->SetLeftMargin(pmargin);
      np->SetBottomMargin(pmargin);
      np->SetFillStyle(4000);

      c1->cd(); np->Draw();
      // get plot
      TList *pl = (TList*)p->GetListOfPrimitives();
      TObject *obj = pl->After(pl->First());
      std::cout << "title-("<<i<<","<<j<<"): " << obj->GetTitle() << "\n";
      TMultiGraph *m = (TMultiGraph*)obj;
      std::string wave_title(m->GetHistogram()->GetTitle());
      m->SetTitle(""); m->GetHistogram()->SetTitle("");
      m->GetHistogram()->SetFillStyle(4000);
      m->GetHistogram()->GetXaxis()->SetTitle("");
      np->cd(); m->Draw("al");
      if(i < j) applyShift(m, shift[i][j])->GetYaxis()->SetRangeUser(-3.3, 3.3);
      if (i == j) {
        TGaxis *ax = new TGaxis(xl, yr, xr, yr,
                                m->GetHistogram()->GetXaxis()->GetBinLowEdge(1),
                                m->GetHistogram()->GetXaxis()->GetBinUpEdge(m->GetHistogram()->GetXaxis()->GetNbins()),
                                m->GetHistogram()->GetXaxis()->GetNdivisions(), "L+");
        ax->SetTitle(""); // ax->SetTitleSize(0.); ax->SetTitleOffset(0.04);
        ax->SetLabelOffset(0.005); ax->SetLabelSize(0.015);
        c1->cd(); ax->Draw();
        // Y axis
        TGaxis *ay = new TGaxis(xl, yr, xl, yl,
                                m->GetHistogram()->GetYaxis()->GetBinLowEdge(1),
                                m->GetHistogram()->GetYaxis()->GetBinUpEdge(m->GetHistogram()->GetYaxis()->GetNbins()),
                                10, "L-");
        ay->SetTitle("");  // ay->SetTitleSize(0.03); ay->SetTitleOffset(+0.9);
        ay->SetLabelOffset(0.027); ay->SetLabelSize(0.015);
        c1->cd(); ay->Draw();
        TLatex lwave_title;
        lwave_title.SetTextSize(0.022);
        lwave_title.SetTextFont(132);
        lwave_title.SetTextAlign(22);
        lwave_title.SetTextAngle(90);
        // lwave_title.SetTextAngle(90);
        lwave_title.DrawLatex(xl-0.04, yr+Xsize/2., titles[i].c_str());
      }
      if (i == 0) {
        TGaxis *ax = new TGaxis(xl, yl, xr, yl,
                                m->GetHistogram()->GetXaxis()->GetBinLowEdge(1),
                                m->GetHistogram()->GetXaxis()->GetBinUpEdge(m->GetHistogram()->GetXaxis()->GetNbins()),
                                m->GetHistogram()->GetXaxis()->GetNdivisions(), "R+");
        ax->SetTitle(""); // ax->SetTitleSize(0.); ax->SetTitleOffset(0.04);
        ax->SetLabelOffset(-0.017); ax->SetLabelSize(0.015);
        c1->cd(); ax->Draw();
        TLatex lwave_title;
        lwave_title.SetTextSize(0.022);
        lwave_title.SetTextFont(132);
        lwave_title.SetTextAlign(22);
        lwave_title.DrawLatex(xl+Xsize/2., yl+0.027, titles[j].c_str());
      }
      if (j == Nsize-1 && i != Nsize-1) {
        TGaxis *ax = new TGaxis(xr, yr, xr, yl,
                                -3.3,//m->GetYaxis()->GetBinLowEdge(1),
                                3.3, //m->GetYaxis()->GetBinUpEdge(m->GetHistogram()->GetYaxis()->GetNbins()),
                                710,//m->GetYaxis()->GetNdivisions(),
                                "L-");
        ax->SetTitle(""); // ax->SetTitleSize(0.); ax->SetTitleOffset(0.04);
        ax->SetLabelOffset(-0.005); ax->SetLabelSize(0.015);
        c1->cd(); ax->Draw();
      }
    }
  }
  TLatex pi3;
  pi3.SetTextSize(0.03);
  pi3.SetTextFont(132);
  pi3.DrawLatex(1-gmargin-0.1, 1-gmargin+0.05,"m_{3#pi} (GeV/c^{2})");

  TLatex phases;
  phases.SetTextSize(0.03);
  phases.SetTextFont(132);
  phases.SetTextAngle(90);
  phases.DrawLatex(1-gmargin+0.04, 1-gmargin-0.24,"Relative phases (rad)");

  TLatex intensities;
  intensities.SetTextSize(0.03);
  intensities.SetTextFont(132);
  intensities.SetTextAngle(90);
  intensities.SetTextAlign(22);
  intensities.DrawLatex(gmargin-0.04, 0.6, "Intensity / (20 MeV)");

  TLatex t_interval;
  t_interval.SetTextSize(0.025);
  t_interval.SetTextFont(132);
  t_interval.DrawLatex(gmargin-0.05, gmargin + 0.15,
                       TString::Format("%s < t' < %s (GeV/c)^{2}",
                                       tranges[ind].c_str(), tranges[ind+1].c_str()));

  TLatex title_reaction;
  title_reaction.SetTextSize(0.025);
  title_reaction.SetTextFont(132);
  title_reaction.DrawLatex(gmargin-0.05, gmargin+0.1,
                       TString::Format("#pi^{-}p #rightarrow #pi^{-}#pi^{+}#pi^{-}p (COMPASS 2008)"));
  TLatex title_text;
  title_text.SetTextSize(0.025);
  title_text.SetTextFont(132);
  title_text.DrawLatex(gmargin-0.05, gmargin+0.07,
                       TString::Format("Mass-independent fit"));
  
  TLatex model;
  title_text.SetTextSize(0.025);
  title_text.SetTextFont(132);
  title_text.SetTextColor(kRed);
  title_text.DrawLatex(gmargin-0.05, gmargin+0.04,
                       TString::Format("Unitarized model"));

  TLatex components;
  components.SetTextSize(0.025);
  components.SetTextFont(132);
  components.DrawLatex(gmargin-0.05, gmargin+0.01,
                       TString::Format("K-matrix components:"));
  TLatex comp1;
  comp1.SetTextSize(0.02);
  comp1.SetTextFont(132);
  comp1.SetTextColor(38);
  comp1.DrawLatex(gmargin-0.04, gmargin-0.02,
                       TString::Format("Pole-1"));
  TLatex comp2;
  comp2.SetTextSize(0.02);
  comp2.SetTextFont(132);
  comp2.SetTextColor(39);
  comp2.DrawLatex(gmargin+0.00, gmargin-0.02,
                       TString::Format("Pole-2"));
  TLatex comp3;
  comp3.SetTextSize(0.02);
  comp3.SetTextFont(132);
  comp3.SetTextColor(40);
  comp3.DrawLatex(gmargin+0.04, gmargin-0.02,
                       TString::Format("Pole-3"));
  TLatex comp4;
  comp4.SetTextSize(0.02);
  comp4.SetTextFont(132);
  comp4.SetTextColor(41);
  comp4.DrawLatex(gmargin+0.08, gmargin-0.02,
                       TString::Format("Pole-4"));

//  TLatex deck;
//  deck.SetTextSize(0.025);
//  deck.SetTextFont(132);
//  deck.SetTextColor(8);
//  deck.DrawLatex(gmargin-0.05, gmargin-0.02,
//                       TString::Format("Nonunitarized long-range process"));

//   for (uint i = 0; i < Nsize; i++) {
//     for (uint j = 0; j < Nsize; j++) {
//       if (shift[i][j] != 0.) {
//         TLatex *lshift = new TLatex();
//         lshift->SetTextSize(0.04);
//         lshift->SetTextFont(132);
//         lshift->DrawLatex(gmargin + (j)*Xsize + Xsize * ( 0.6 ),
//                           gmargin + (Nsize-i-1)*Ysize+Ysize * ( 0.6 ),
//                           TString::Format("%f^{#circ}", shift[i][j]));
//       }
//     }
//   }  
  
  c1->SaveAs(TString::Format("/tmp/SDM%02d.pdf", ind));
  return 0;
}


TMultiGraph * applyShift(TMultiGraph *m, int shift) {
  // std::cout << "shift " << shift << std::endl;
  TList *ll = m->GetListOfGraphs();
  TGraph *g1 = (TGraph*)(ll->First());
  TGraph *g2 = (TGraph*)(ll->After(g1));
  TGraph *g3 = (TGraph*)(ll->After(g2));
  double rshift = (shift != 0) ? M_PI/shift : 0;
  for (int i=0; i < g1->GetN(); i++) g1->GetY()[i] = toPi(g1->GetY()[i] + rshift);
  for (int i=0; i < g2->GetN(); i++) g2->GetY()[i] = toPi(g2->GetY()[i] + rshift);
  for (int i=0; i < g3->GetN(); i++) g3->GetY()[i] = toPi(g3->GetY()[i] + rshift);
  
//   // ***********
  if (shift == 0) return m;
  TPad *np = new TPad("ttt", "noTitle", 0.0, 0.0, 1.0, 1.0, kRed); // i+j+1
  np->SetFillStyle(4000);
  np->Draw(); np->cd();
  TLatex tx;
  tx.SetTextSize(0.15);
  tx.SetTextFont(132);
  tx.SetTextAlign(32);
  tx.DrawLatex(0.97, 0.9,
               (shift > 0) ?
               ((shift == 1) ? "+#pi" : TString::Format("+#pi/%d", shift)) :
               ((shift == -1) ? "-#pi" : TString::Format("-#pi/%d", shift)));  
  return m;
}

double toPi(double v) {
  double ov = v;
  if (v >= M_PI) ov = v - 2.*M_PI;
  if (v < -M_PI) ov = v + 2.*M_PI;
  return ov;
}
