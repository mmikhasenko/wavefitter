// Copyright of libconfig++
// My modification of standard example

// [] Use to compile
// [] g++ -o read_model_settings \
//       /localhome/mikhasenko/Tools/libconfig-1.5/lib/libconfig++.so \
//       -I/localhome/mikhasenko/Tools/libconfig-1.5/include read_model_settings.cc

#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <vector>

#include <libconfig.h++>

#include "MIsobar.h"
#include "MIsobarChannel.h"
#include "MParKeeper.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MRelationHolder.h"

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "mstructures.hh"

int main(int argc, char *argv[]) {
  const char *output_file = "updated.cfg";
  libconfig::Config cfg;

  // Read the file. If there is an error, report it and exit.
  try {
    if (argc > 1) cfg.readFile(argv[1]);
    else
      return EXIT_FAILURE;
  }
  catch(const libconfig::FileIOException &fioex) {
    std::cerr << "I/O error while reading file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const libconfig::ParseException &pex) {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    return(EXIT_FAILURE);
  }

  // Get the store name.
  try {
    std::string name = cfg.lookup("name");
    std::cout << "Store name: " << name << std::endl << std::endl;
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "No 'name' setting in configuration file." << std::endl;
  }

  const libconfig::Setting& root = cfg.getRoot();

  std::vector<DP> whole_data;
  // Get data fiels
  try {
    const libconfig::Setting &data = root["data"];
    std::string path;
    data.lookupValue("path", path);
    std::cout << path << "\n";

    const libconfig::Setting &dfiles = root["data"]["points"];
    const uint count = dfiles.getLength();

    whole_data.resize(count);
    for (uint i = 0; i < count; ++i) {
      const libconfig::Setting &dfile = dfiles[i];

      // Only output the record if all of the expected fields are present.
      std::string title, type, file_name;

      if (  !(dfile.lookupValue("type", type)
              && dfile.lookupValue("file_name", file_name)
              && dfile.lookupValue("title", title)))
        continue;

      const libconfig::Setting &trust_range = dfile["trust_range"];
      double tv1 = trust_range[0];
      double tv2 = trust_range[1];
      const libconfig::Setting &fit_range = dfile["fit_range"];
      double fv1 = fit_range[0];
      double fv2 = fit_range[1];

      std::cout << "READ: " << std::left << file_name << "  "
                << title << "  ("
                << fv1 << ", " << fv2 << "), ("
                << tv1 << ", " << tv2 << ")"
                << std::endl;

      TGraphErrors *g = (type == "txt") ? new TGraphErrors((path+file_name).c_str()) : 0;
      if (!g) {std::cerr << "Warning <> " << file_name << "not found!" << std::endl; continue; }
      std::string gname = path;
      gname.erase(std::remove(gname.begin(), gname.end(), '/'), gname.end());
      // convert to DP fortat
      if (g->GetN() == 0) {std::cerr << "Error<> Nothing is read from data file!" << "\n"; return EXIT_FAILURE;}
      for (uint j = 0; j < g->GetN(); j++) {
        if (g->GetX()[j] < tv1 || g->GetX()[j] > tv2) continue;
        whole_data[i].data.push_back(data_point{
            g->GetX()[j],
              g->GetY()[j],
              g->GetEX()[j]
              });
//        std::cout << g->GetX()[j] << " "
//                  << g->GetY()[j] << " "
//                  << g->GetEX()[j] << "\n";
      }
      whole_data[i].name = gname.c_str();
      whole_data[i].title = title.c_str();
      whole_data[i].lrange = fv1;
      whole_data[i].rrange = fv2;

      delete g;
    }
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"data\" secton!" << std::endl;
    return EXIT_FAILURE;
  }

  const uint Nhist = whole_data.size();
  TCanvas c1("c1");
  c1.DivideSquare(Nhist);
  for (uint i = 0; i < Nhist; i++) {
    c1.cd(i+1); draw(whole_data[i])->Draw("ap");
  }
  c1.SaveAs("/tmp/1.pdf");
  std::cout << "\n";

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// m o d e l   a d j u s t m e n t ////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  // channels
  std::vector<MChannel*> iset;
  // standard isobars
  MIsobar rho_iso(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
  MIsobar  f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.);
  // model content is vector to fill K-matrix
  std::vector<std::pair<std::string, std::string> > model_content;

  try {
    const libconfig::Setting &modelT = root["modelT"];

    const libconfig::Setting &channels = modelT["channels"];
    uint count = channels.getLength();

    for (uint i = 0; i < count; ++i) {
      const libconfig::Setting &channel = channels[i];

      std::string type = "?";
      int L = 0;
      double R = 3.0;
      const libconfig::Setting &spart = channel["particles"];
      std::string particles[2] = {spart[0], spart[1]};

      if (  !(channel.lookupValue("type", type) &&
              channel.lookupValue("L", L) &&
              channel.lookupValue("size", R)))
        continue;

      std::cout << "READ: " << std::left << type << ": "
                << particles[0] << " " << particles[1] << " "
                << (std::vector<char>{'S', 'P', 'D', 'F', 'G', 'H'})[L] << "-wave"
                << std::endl;

      if (type != "quasi-two-body") {
        std::cerr << "Error: type of model channel can be only quasi-two-body";
        return EXIT_FAILURE;
      }
      if (particles[0] == "rho") {
        MIsobarChannel *mCh = new MIsobarChannel(rho_iso, PI_MASS, L);
        mCh->makeLookupTable(mCh->sth(), 10., 100);
        iset.push_back(mCh);
      } else if (particles[0] == "f2") {
        MIsobarChannel *mCh = new MIsobarChannel(f2_iso, PI_MASS, L);
        mCh->makeLookupTable(mCh->sth(), 10., 100);
        iset.push_back(mCh);
      } else {
        std::cerr << "Error: isobar is not rho/f2. Only them are available.";
        return EXIT_FAILURE;
      }
    }
    for (auto && it : iset) it->makeDisperseLookupTable(0.01, 10., 100);

    // content of the model
    const libconfig::Setting &content = modelT["content"];
    count = content.getLength();
    for (uint i = 0; i < count; ++i) {
      const libconfig::Setting &kcomp = content[i];

      std::string type = "?",
        mass = "?",
        couplings = "?";

      if (  !(kcomp.lookupValue("type", type) &&
              kcomp.lookupValue("mass", mass) &&
              kcomp.lookupValue("couplings", couplings)))
        continue;

      std::cout << "READ: " << type << " ("
                << mass << ", " << couplings << "*)"
                << std::endl;
      if (type != "pole") {
        std::cerr << "Error: model blok is not a pole. Only poles are available.";
        return EXIT_FAILURE;
      }
      model_content.push_back(std::make_pair(mass, couplings));
    }  // icount
    // create the model
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"model\" secton!" << std::endl;
    return EXIT_FAILURE;
  }
  
  MmatrixK km(iset, model_content.size());
  MProductionPhysics pr(iset);
  pr.addScattering([&](double s)->b::matrix<cd>{return km.getValue(s);});
  pr.addShortRange();
  MParKeeper::gI()->printAll();

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// p a r a m e t e r s ///////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  try {
    const libconfig::Setting &list_of_parameters = root["parameters"];

    const libconfig::Setting &starting_values = list_of_parameters["start_value"];
    const uint count = starting_values.getLength();

    if (count != MParKeeper::gI()->nPars()) {
      std::cerr << "Warning: number of the specified parameters does not much to number of expected parameters!\n";
      // return EXIT_FAILURE;
    }
    std::vector<std::string> pars_names(count);
    std::vector<double> pars_values(count);
    std::vector<std::pair<double, double> > pars_ranges(count);
    for (uint i = 0; i < count; i++) {
      const libconfig::Setting &iPar = starting_values[i];
      const char* name = iPar[0];
      pars_names[i] = std::string(name);
      pars_values[i] = iPar[1];
      pars_ranges[i].first = iPar[2];
      pars_ranges[i].second = iPar[3];
      std::cout << "READ: " << name << ": " << pars_values[i] << " in ("
                << pars_ranges[i].first << ", " << pars_ranges[i].second << ")\n";
    }
    // create pool
    MParKeeper::gI()->makePool(pars_names);
    MParKeeper::gI()->pset(pars_values.data());
    for (uint j=0; j < MParKeeper::gI()->pnPars(); j++)
      MParKeeper::gI()->psetRange(j,
                                  pars_ranges[j].first,
                                  pars_ranges[j].second);
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in the \"parameters\" section" << std::endl;
    return EXIT_FAILURE;
  }

  const uint pnPars = MParKeeper::gI()->pnPars();

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// r e l a t i o n ///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  try {
    const libconfig::Setting &adjustment = root["adjustment"];

    const libconfig::Setting &list_of_relations = adjustment["relation"];
    const uint count = list_of_relations.getLength();

    for (uint i = 0; i < count; ++i) {
      const libconfig::Setting &iRel = list_of_relations[i];

      uint jData = iRel[0];
      if (jData >= whole_data.size()) std::cerr << "Error: jData<0 || jData >= Nhist." << std::endl;
      std::string type = iRel[1];
      if (type == "I@") {
        uint iCh = iRel[2];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh](double e)->double{
            auto v = pr.getValue(e*e);
            return norm(v(iCh))*iset[iCh]->rho(e*e);
          });        
      } else if (type == "Phi@") {
        uint iCh1 = iRel[2][0];
        uint iCh2 = iRel[2][1];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh1, iCh2](double e)->double{
            auto v = pr.getValue(e*e);
            return arg(v(iCh1)*conj(v(iCh2)));
          });
      } else if (type == "SinPhi@") {
        uint iCh1 = iRel[2][0];
        uint iCh2 = iRel[2][1];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh1, iCh2](double e)->double{
            auto v = pr.getValue(e*e);
            return sin(arg(v(iCh1)*conj(v(iCh2))));
          });
      } else if (type == "CosPhi@") {
        uint iCh1 = iRel[2][0];
        uint iCh2 = iRel[2][1];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh1, iCh2](double e)->double{
            auto v = pr.getValue(e*e);
            return cos(arg(v(iCh1)*conj(v(iCh2))));
          });
      }
    }
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"relation\" secton" << std::endl;
    return EXIT_FAILURE;
  }

  const uint Nrels = MRelationHolder::gI()->Nrels();
  MRelationHolder::gI()->Print();
  TCanvas *canva = new TCanvas("canva","title");
  canva->DivideSquare(Nrels);

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// f i t  s e t t i n g s /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  try {
    const libconfig::Setting &fit_settings = root["fit_settings"];

    const libconfig::Setting &strategy = fit_settings["strategy"];

    const uint nAttempts = fit_settings["nAttempts"];

    /*********************************** Fit itself *****************************************/
    /****************************************************************************************/
    TFile fout("/tmp/test.fit.result.root", "RECREATE");
    TTree tout("tout", "Results of fit");
    // set branches
    tout.Branch("can", "TCanvas", &canva);
    double chi2 = 0; tout.Branch("chi2", &chi2);
    double pars_mirrow[MParKeeper::gI()->nPars()];
    for (int i=0; i < MParKeeper::gI()->nPars(); i++)
      tout.Branch(MParKeeper::gI()->getName(i).c_str(), &pars_mirrow[i]);
    // to copy to array from where it is copied to tree
    for (uint e = 0; e < nAttempts; e++) {
      std::cout << "---------- Attempt " << e << " -----------" << std::endl;
      std::cout << "------------------------------------------" << std::endl;

      /**************************************MINIMIZE******************************************/
      // Build minimizer
      ROOT::Math::Minimizer* min =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

      // set tolerance , etc...
      min->SetMaxFunctionCalls(100000);
      min->SetTolerance(0.001);
      min->SetStrategy(1);
      min->SetPrintLevel(3);
      min->Options().Print();

      // Create funciton wrapper for minmizer a IMultiGenFunction type
      ROOT::Math::Functor functor([&](const double *pars)->double {
          MParKeeper::gI()->pset(pars);
          km.RecalculateNextTime();
          pr.RecalculateNextTime();
          return MRelationHolder::gI()->CalculateChi2();
        }, pnPars);
      min->SetFunction(functor);
      // specify parameters
      MParKeeper::gI()->randomizePool();
      MParKeeper::gI()->printAll();
      for (uint i=0; i < pnPars; i++) min->SetVariable(i,
                                                      MParKeeper::gI()->pgetName(i),
                                                      MParKeeper::gI()->pget(i),
                                                      0.1);
      /*ooooooooooooooooooooooooooooooooooooooo Fit itself ooooooooooooooooooooooooooooooooooo*/
      // step fit
      const uint count = strategy.getLength();
      for (uint i = 0; i < count; i++) {
        const libconfig::Setting &fit_step = strategy[i];
        const libconfig::Setting &relations = fit_step["relations_to_fit"];
        const libconfig::Setting &pars = fit_step["pars_to_vary"];
        // adjust which relation to fit
        const uint Nrelations = relations.getLength();
        MRelationHolder::gI()->passiveAll();
        for (uint r=0; r < Nrelations; r++) MRelationHolder::gI()->activateRelation(relations[r]);
        MRelationHolder::gI()->Print();
        // adjust which parameters to vary
        const uint nPars_to_vary = pars.getLength();
        // loop over all parameters for make fix them
        for (uint r=0; r < pnPars; r++) min->FixVariable(r);
        for (uint r=0; r < nPars_to_vary; r++) {
          const std::string &pname = pars[r];
          min->ReleaseVariable(MParKeeper::gI()->pgetIndex(pname));
        }
        // minimize
        min->Minimize();
      }
      // Plot all
      for (uint i=0; i < Nrels; i++) {
        const DP & data = MRelationHolder::gI()->GetRelation(i).data;
        std::function<double(double)> func = MRelationHolder::gI()->GetRelation(i).func;
        // draw
        canva->cd(i+1)->Clear();
        draw(data)->Draw("ap");
        SET1(draw(func,
                  (data.data.begin())->x, (--data.data.end())->x, 200),
             SetLineColor(kOrange) )->Draw("l");  // same
        if (MRelationHolder::gI()->relationStatus(i))
          SET1(draw(func,
                    data.lrange, data.rrange, 200),
               SetLineColor(kRed) )->Draw("l");
      }
      MParKeeper::gI()->printAll();
      canva->SaveAs(TString::Format("/tmp/e%d.pdf", e));
      delete min;

      // Fill result to tree
      // to copy to array from where it is copied to tree
      memcpy(pars_mirrow,
             MParKeeper::gI()->get().data(),
             sizeof(double)*MParKeeper::gI()->nPars());
      chi2 = MRelationHolder::gI()->CalculateChi2();
      tout.Fill();
    }

    tout.Write();
    fout.Close();
    delete canva;
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"fit\" secton" << std::endl;
    return EXIT_FAILURE;
  }

  // well done

  return(EXIT_SUCCESS);
}
