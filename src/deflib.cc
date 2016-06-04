#include <deflib.h>
#include <TGraphErrors.h>

//TMultiGraph *split_data(const TH1D* h, double low, double up) {
//  TGraphErrors *gray = new TGraphErrors(0); gray->SetLineColor(kGray); gray->SetMarkerColor(kGray);
//  TGraphErrors *dark = new TGraphErrors(0);
//
//  int Nbins = h->GetXaxis()->GetNbins();
//  double dm = h->GetBinCenter(2) - h->GetBinCenter(1);
//  for(int i=0;i<Nbins;i++) {
//    double mass = h->GetBinCenter(i+1);
//    if(mass<low || mass >up) {
//      if(i!=0&&h->GetBinContent(i+1)==0&&h->GetBinError(i+1)==0) continue;
//      gray->Set(gray->GetN()+1);
//      gray->SetPoint     (gray->GetN()-1,mass,h->GetBinContent(i+1));
//      gray->SetPointError(gray->GetN()-1,dm  ,h->GetBinError  (i+1));
//    } else {
//      dark->Set(dark->GetN()+1);
//      dark->SetPoint     (dark->GetN()-1,mass,h->GetBinContent(i+1));
//      dark->SetPointError(dark->GetN()-1,dm  ,h->GetBinError  (i+1));
//    }
//  }
//  TMultiGraph *m = new TMultiGraph();
//  m->Add(gray,"pz");
//  m->Add(dark,"pz");
//  m->SetTitle(TString::Format("%s;%s;%s",h->GetTitle(),h->GetXaxis()->GetTitle(),h->GetYaxis()->GetTitle()));  
//  return m;
//}
//
//double normPhase(double a) {
//  if(a>M_PI)  return normPhase(a-2*M_PI);
//  if(a<-M_PI) return normPhase(a+2*M_PI);
//  return a;
//}

