// "1-(1++)0+ f0(980) pi P"
// "1-(4++)1+ rho pi G"
// "1-(2++)1+ rho pi D"

int extraction_of_graphs_from_Stefan(const char *outname,
                                     const char *hname,
                                     const char *what_name = "intensity",
                                     const char *tslice = "0.100000GeV_0.112853GeV") {
  gDirectory->cd("TBins");   gDirectory->ls();
  gDirectory->cd(tslice);    gDirectory->ls();
  gDirectory->cd(what_name); gDirectory->ls();
  THStack *hs = (THStack*)gDirectory->Get(hname);
  for (uint ng = 0; ng < hs->GetNhists(); ng++) {
    std::ofstream ofstr(TString::Format("%s_%i", outname, ng).Data());
    TH1D *h = (TH1D*)hs->GetHists()->At(ng);
    for(int i=1;i <= h->GetXaxis()->GetNbins();i++)
      ofstr << h->GetBinCenter(i) << " " << h->GetBinContent(i) << "\n";
    ofstr.close();
  }
  return 0;
}
