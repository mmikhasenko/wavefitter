int extraction_of_graphs_from_Mathias(TCanvas *c, int ip, const char *outname) {
  TPad *p = (TPad*)c->GetListOfPrimitives()->At(ip);
  TMultiGraph *m = (TMultiGraph*)p->GetListOfPrimitives()->At(1);
  for (uint ng = 1; ng < m->GetListOfGraphs()->GetSize(); ng++) {
    std::ofstream ofstr(TString::Format("%s_%i", outname, ng).Data());
    TGraph *gr = (TGraph*)m->GetListOfGraphs()->At(ng);
    for(int i=0;i<gr->GetN();i++)
      ofstr << gr->GetX()[i] << " " << gr->GetY()[i] << "\n";
    ofstr.close();
  }
  return 0;
}
