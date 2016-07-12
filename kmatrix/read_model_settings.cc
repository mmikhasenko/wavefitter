// Copyright of libconfig++
// My modification of standard example

// [] Use to compile
// [] g++ -o read_model_settings \
//       /localhome/mikhasenko/Tools/libconfig-1.5/lib/libconfig++.so \
//       -I/localhome/mikhasenko/Tools/libconfig-1.5/include read_model_settings.cc

#include <iostream>
#include <fstream>
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
#include "MDeck.h"

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TMultiGraph.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "mstructures.h"

int main(int argc, char *argv[]) {
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

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////// d a t a   s p e c i f i c a t i o n ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  std::vector<DP> whole_data;
  // Get data fiels
  try {
    const libconfig::Setting &data = root["data"];

    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cout << "/////////////// Data specification: //////////////////\n";

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
      for (uint j = 0; j < (uint)g->GetN(); j++) {
        if (g->GetX()[j] < tv1 || g->GetX()[j] > tv2) continue;
        whole_data[i].data.push_back(data_point{
            g->GetX()[j],
              g->GetY()[j],
              g->GetEX()[j]
              });
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

    std::cout << "\n\n";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cout << "/////////////// Model adjustment: ////////////////////\n";

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

  const uint Jsector = 2;
  MmatrixK km(iset, model_content.size());


  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////// Production model: /////////////////////////////////////////////////////
  /////////////////////// short range, long range, unitarisation ////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  std::vector<MProductionPhysics*> vpr;  // (iset);

  try {
    const libconfig::Setting &modelsA = root["modelA"];
    const uint Nmodels = modelsA.getLength();

    std::cout << "\n\n";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cout << "/////////////// Production model: ////////////////////\n";
    std::cout << "--> " << Nmodels << " production models will be constructed.\n";

    for (uint imodelA = 0; imodelA < Nmodels; imodelA++) { 
      const libconfig::Setting &modelA = modelsA[imodelA];

      vpr.push_back(new MProductionPhysics(iset));
      // set a reference and use as it was done
      MProductionPhysics &pr = *(vpr[imodelA]);

      if ( modelA.exists("scattering") )
        pr.addScattering([&](double s)->b::matrix<cd>{return km.getValue(s);});

      // short_range
      if (modelA.exists("short_range")) {
        const libconfig::Setting &short_range = modelA["short_range"];
        std::string type;
        double rhc, slope;
        std::vector<std::pair<int,std::string> > powers;
        if ( short_range.lookupValue("type", type) &&
             short_range.lookupValue("rhc", rhc) &&
             short_range.lookupValue("slope", slope) &&
             short_range.exists("powers")) {
          std::cout << "READ: " << type << ", s0 = "
                    << rhc << ", #alpha = " << slope << ")"
                    << std::endl;
          const libconfig::Setting &orders = short_range["powers"];
          const uint count = orders.getLength();
          std::cout << "Powers of expantion: ";
          for (uint i = 0; i < count; i++) {
            int power(orders[i][0]);
            std::string name = orders[i][1];
            std::cout << name << " " << power;
            powers.push_back(std::make_pair(power, name));
          }
          std::cout << "\n";
        }
        if (powers.size() != 1) {std::cout << "Warning: pietarinen expantion is not inplemented. A constant will be used"; }
        // temperary fix
        if (powers.size() > 0) pr.addShortRange(powers[0].second); /* to be improved */
        else pr.addShortRange(); /* to be improved */
      }

      std::vector<MDeck*> vdeck;
      if (modelA.exists("long_range")) {
        const libconfig::Setting &long_range = modelA["long_range"];
        std::string type;
        uint Sp; int M;
        double R;
        double tP;
        if ( long_range.lookupValue("type", type) &&
             long_range.lookupValue("damping_R", R) &&
             long_range.lookupValue("pomeron_virtuality", tP) &&
             long_range.lookupValue("pomeron_S", Sp) &&
             long_range.lookupValue("pomeron_M", M)) {
          if (type != "Deck") {
            std::cerr << "Error: long range interaction is not 'Deck', but no other options are available!";
            return 0;
          }
          std::cout << "Deck projections will be constructed with parameters: J = " << Jsector
                    << ", Sp = " << Sp << ", M = " << M << ", R = " << R << "\n";
          // check if there is lookup lable already
          std::vector<std::string> long_range_lookup_path;
          if (long_range.exists("long_range_lookup")) {
            const libconfig::Setting &long_range_lookup = long_range["long_range_lookup"];
            if (long_range_lookup.getLength() == static_cast<int>(iset.size())) {
              long_range_lookup_path.resize(iset.size());
              for (uint i=0; i < iset.size(); i++) {
                const std::string &name = long_range_lookup[i];
                long_range_lookup_path[i] = std::string(name);
              }
              std::cout << "--> Files for long range lookup tables are found!\n";
            } else {
              std::cout << "Warning: long_range.getLength() != static_cast<int>(iset.size())\n"; 
            }
          } else {
            std::cout << "--> Files for long range lookup tables are not found!\n";
          }
          // create deck and add functions
          vdeck.resize(iset.size());
          std::vector<std::function<cd(double)> > getB(iset.size());
          TCanvas c3("c3"); c3.DivideSquare(iset.size());
          for (uint i=0; i < iset.size(); i++) {
            MIsobarChannel *ich = dynamic_cast<MIsobarChannel*>(iset[i]);
            const MIsobar &iso = ich->getIsobar();
            vdeck[i] = new MDeck(POW2(PI_MASS), tP, iso.GetM(), POW2(PI_MASS), POW2(PI_MASS),
                                 Jsector, iset[i]->GetL(), Sp, -M,
                                 iso.GetL(), 0, R);
            vdeck[i]->Print();

            // check if it is possible to load it from somewhere
            if (long_range_lookup_path.size()) {
              TGraph *igr = new TGraph(long_range_lookup_path[i].c_str());
              if (igr->GetN()) {
                std::vector<std::pair<double, double> > ltable;
                for (int j = 0; j < igr->GetN(); j++) ltable.push_back(std::make_pair(igr->GetX()[j], igr->GetY()[j]));
                vdeck[i]->setLookupTable(ltable);
                std::cout << "-> Lookup table is filled\n";
              } else {
                std::cout << "Warning: igr.GetN() = 0\n";
                vdeck[i]->makeLookupTable(iso, ich->getBachelorMass(), ich->sth(), POW2(2.5), 30);
                const std::vector<std::pair<double, double> > &ltable = vdeck[i]->getLookupTable();
                // save to file
                std::ofstream myfile(long_range_lookup_path[i]);
                if (myfile.is_open()) {
                  for (auto && it : ltable) myfile << it.first << " "  << it.second << "\n";
                  myfile.close();
                }
                std::cout << "lookup table is saved at " << long_range_lookup_path[i] << "\n";
              }
            } else {
              vdeck[i]->makeLookupTable(iso, ich->getBachelorMass(), ich->sth(), POW2(2.5), 30);
            }

            //
            MDeck *iD = vdeck[i];
            getB[i] = [&, iD](double s)->cd { return iD->getPrecalculated(s); };
            c3.cd(i+1);
            draw([&, i](double s)->double{return vdeck[i]->getPrecalculated(s);}, 1.0, POW2(4.2))->Draw("al");
          }
          c3.SaveAs("/tmp/deck.pdf");
          // finally add Long Range
          if (long_range.exists("par_name")) {
            std::string par_name = long_range["par_name"];
            pr.addLongRange(getB, par_name);
          } else pr.addLongRange(getB);
        }
      }
      if (modelA.exists("unitarisation")) pr.unitarize();
    }
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"production model\" secton!" << std::endl;
    return EXIT_FAILURE;
  }
  MParKeeper::gI()->printAll();

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// p a r a m e t e r s ///////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  try {
    const libconfig::Setting &list_of_parameters = root["parameters"];

    std::cout << "\n\n";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cout << "/////////////// Parameters business: /////////////////\n";

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

    std::cout << "\n\n";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cout << "/////////////// Relations business: //////////////////\n";

    const libconfig::Setting &list_of_relations = adjustment["relation"];
    const uint count = list_of_relations.getLength();

    for (uint i = 0; i < count; ++i) {
      const libconfig::Setting &iRel = list_of_relations[i];

      uint jData = iRel[0];
      if (jData >= whole_data.size()) std::cerr << "Error: jData<0 || jData >= Nhist." << std::endl;
      std::string type = iRel[1];
      if (type == "I@") {
        uint iModel = iRel[2][0];
        uint iCh = iRel[2][1];
        MProductionPhysics *pr = vpr[iModel];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return norm(v(iCh))*iset[iCh]->rho(e*e);
          });
      } else if (type == "Re@") {
        uint iModel = iRel[2][0];
        uint iCh = iRel[2][1];
        MProductionPhysics *pr = vpr[iModel];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return real(v(iCh))*iset[iCh]->rho(e*e);
          });
      } else if (type == "Im@") {
        uint iModel = iRel[2][0];
        uint iCh = iRel[2][1];
        MProductionPhysics *pr = vpr[iModel];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return imag(v(iCh))*iset[iCh]->rho(e*e);
          });
      } else if (type == "Phi@") {
        uint iModel0 = iRel[2][0][0]; uint iCh0 = iRel[2][0][1];
        uint iModel1 = iRel[2][1][0]; uint iCh1 = iRel[2][1][1];
        if (iModel0 != iModel1) { std::cerr << "Error: iModel0!=iModel1!\n"; return 1;}
        MProductionPhysics *pr = vpr[iModel1];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh0, iCh1, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return arg(v(iCh0)*conj(v(iCh1)));
          });
      } else if (type == "SinPhi@") {
        uint iModel0 = iRel[2][0][0]; uint iCh0 = iRel[2][0][1];
        uint iModel1 = iRel[2][1][0]; uint iCh1 = iRel[2][1][1];
        if (iModel0 != iModel1) { std::cerr << "Error: iModel0!=iModel1!\n"; return 1;}
        MProductionPhysics *pr = vpr[iModel1];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh0, iCh1, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return sin(arg(v(iCh0)*conj(v(iCh1))));
          });
      } else if (type == "CosPhi@") {
        uint iModel0 = iRel[2][0][0]; uint iCh0 = iRel[2][0][1];
        uint iModel1 = iRel[2][1][0]; uint iCh1 = iRel[2][1][1];
        if (iModel0 != iModel1) { std::cerr << "Error: iModel0!=iModel1!\n"; return 1;}
        MProductionPhysics *pr = vpr[iModel1];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh0, iCh1, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return cos(arg(v(iCh0)*conj(v(iCh1))));
          });
      } else if (type == "ReInterf@") {
        uint iModel0 = iRel[2][0][0]; uint iCh0 = iRel[2][0][1];
        uint iModel1 = iRel[2][1][0]; uint iCh1 = iRel[2][1][1];
        if (iModel0 != iModel1) { std::cerr << "Error: iModel0!=iModel1!\n"; return 1;}
        MProductionPhysics *pr = vpr[iModel1];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh0, iCh1, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return real(v(iCh0)*conj(v(iCh1)))*sqrt(iset[iCh0]->rho(e*e)*iset[iCh1]->rho(e*e));
          });
      } else if (type == "ImInterf@") {
        uint iModel0 = iRel[2][0][0]; uint iCh0 = iRel[2][0][1];
        uint iModel1 = iRel[2][1][0]; uint iCh1 = iRel[2][1][1];
        if (iModel0 != iModel1) { std::cerr << "Error: iModel0!=iModel1!\n"; return 1;}
        MProductionPhysics *pr = vpr[iModel1];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh0, iCh1, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return imag(v(iCh0)*conj(v(iCh1)))*sqrt(iset[iCh0]->rho(e*e)*iset[iCh1]->rho(e*e));
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
  TCanvas *canva = new TCanvas("canva", "title");

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// p l o t  s e t t i n g s /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  try {
    if (root.exists("plot_settings")) {
      const libconfig::Setting &plot_settings = root["plot_settings"];

      std::cout << "\n\n";
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      std::cout << "/////////////// Plot settings: ///////////////////////\n";

      std::string type;
      if (!plot_settings.lookupValue("type", type)) type = "model";
      const libconfig::Setting &mapping = plot_settings["mapping"];
      const uint NrelsToPlot = mapping.getLength();
      canva->Clear(); canva->DivideSquare(NrelsToPlot);

      // if file is specified
      //  - open file
      //  - find tree thee
      //  - load perameters from particular entry of the tree
      if (plot_settings.exists("path_to_fit_results") && plot_settings.exists("entry_fit_result")) {
        std::string fres_name; plot_settings.lookupValue("path_to_fit_results", fres_name);
        uint entry =  plot_settings["entry_fit_result"];  // throw exception it something is wrong
        TFile *fres = TFile::Open(fres_name.c_str());
        if (fres) {
          std::cout << "File with results successfully opened!\n";
          TTree *tres; gDirectory->GetObject("tout", tres);
          if (tres) {
            const uint Npars = MParKeeper::gI()->nPars();
            double pars[Npars];
            for (uint i = 0; i < Npars; i++) {
              const std::string & name = MParKeeper::gI()->getName(i);
              // check if it is at list of branches
              tres->SetBranchAddress(name.c_str(), &pars[i]);
            }
            // set values from tree and set to keeper
            tres->GetEntry(entry);
            for (uint i = 0; i < Npars; i++) MParKeeper::gI()->set(i, pars[i]);
          } else {
            std::cerr << "Tree with results not found by name 'tout'!\n";
            return 1;
          }
        } else {
          std::cerr << "File with results specified but not found!\n";
          return 1;
        }
      }

      std::vector<TMultiGraph*> mgr(NrelsToPlot);
      // add data in order by mapping
      for (uint i=0; i < NrelsToPlot; i++) {
        if (mapping[i].getLength() == 0) continue;
        const uint iPad = mapping[i][0];
        mgr[i] = new TMultiGraph();
        // draw data
        const DP & data = MRelationHolder::gI()->GetRelation(iPad).data;
        if ( type.find("data") != std::string::npos ) mgr[i]->Add(draw(data), "p");
      }

      // add model model
      if ( type.find("model") != std::string::npos ) {
        if (plot_settings.exists("what_to_plot")) {
          const libconfig::Setting &what_to_plot = plot_settings["what_to_plot"];
          const uint Ncurves = what_to_plot.getLength();
          std::cout << "--> " << Ncurves << " model settings are found to be plotted.\n";
          for (uint j = 0; j < Ncurves; j++) {
            std::string title; uint color;
            if (what_to_plot[j].lookupValue("title", title) &&
                what_to_plot[j].lookupValue("color", color) &&
                what_to_plot[j].exists("set_to_zero")) {
              const libconfig::Setting &set_to_zero = what_to_plot[j]["set_to_zero"];
              const uint NparsToZero = set_to_zero.getLength();
              // set to zero
              std::vector<std::pair<uint, double> > parsBackUp;
              for (uint p = 0; p < NparsToZero; p++) {
                const char *pname = set_to_zero[p];
                uint ip = MParKeeper::gI()->getIndex(std::string(pname));
                double vp = MParKeeper::gI()->get(ip);
                parsBackUp.push_back(std::make_pair(ip, vp));
                MParKeeper::gI()->set(ip, 0.);
              }
              km.RecalculateNextTime();
              for (auto & pr : vpr) pr->RecalculateNextTime();
              // message
              std::cout << "--> Settings " << j << " <" << title << "> : \n";
              MParKeeper::gI()->printAll();

              // fill vector of historrams where the model will be plotted on
              std::vector<uint> vme;
              if (what_to_plot[j].exists("mapping_elements")) {
                const libconfig::Setting &mapping_elements = what_to_plot[j]["mapping_elements"];
                for (int me = 0; me < mapping_elements.getLength(); me++) vme.push_back(uint(mapping_elements[me]));
              } else {
                for (uint i=0; i < NrelsToPlot; i++) vme.push_back(i);
              }
              // add model curves to plot
              for (uint i : vme) {
                if (mapping[i].getLength() == 0) continue;
                const uint iPad = mapping[i][0];
                // get data and model function
                const DP & data = MRelationHolder::gI()->GetRelation(iPad).data;
                std::function<double(double)> func = MRelationHolder::gI()->GetRelation(iPad).func;
                // plot
                mgr[i]->Add(
                            SET3(draw(func,
                                      (data.data.begin())->x , (--data.data.end())->x,
                                      100),
                                 SetLineStyle(2),
                                 SetLineColor(color),
                                 SetTitle("") ), "l");
                mgr[i]->Add(
                            SET2(draw(func,
                                      data.lrange, data.rrange,
                                      100),
                                 SetLineColor(color),
                                 SetTitle(title.c_str()) ), "l");
              }
              // set back to nominal value
              for (uint p = 0; p < NparsToZero; p++)
                MParKeeper::gI()->set(parsBackUp[p].first, parsBackUp[p].second);
            }
          }
        }
      }

      // finally plot
      for (uint i=0; i < NrelsToPlot; i++) {
        if (mapping[i].getLength() == 0) continue;
        // setTitle
        const uint iPad = mapping[i][0];
        const DP & data = MRelationHolder::gI()->GetRelation(iPad).data;
        mgr[i]->SetTitle(data.title.c_str());
        // draw
        canva->cd(i+1); mgr[i]->Draw("a");
      }

      // save to pdf
      std::string fplot_name = "/tmp/default_plot.read_model_settings.pdf";
      if (!plot_settings.lookupValue("fplot_name", fplot_name))
        std::cerr << "Warning: fplot_name is not specified. A default name wil be used.\n";
      canva->SaveAs(fplot_name.c_str());
    }
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"plot_settings\" secton" << std::endl;
    return EXIT_FAILURE;
  }


  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// f i t  s e t t i n g s /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  canva->Clear(); canva->DivideSquare(Nrels);
  try {
    if (root.exists("fit_settings")) {
      const libconfig::Setting &fit_settings = root["fit_settings"];

      std::cout << "\n\n";
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      std::cout << "/////////////// Fit settings: ////////////////////////\n";

      const libconfig::Setting &strategy = fit_settings["strategy"];

      const uint nAttempts = fit_settings["nAttempts"];

      const std::string &dout_name = fit_settings["dout_name"];

      /*********************************** Fit itself *****************************************/
      /****************************************************************************************/
      TFile *fout = new TFile(TString::Format("%s/fit.results.root", dout_name.c_str()), "RECREATE");
      if (!fout) { std::cerr << "no fout acceptable!\n"; return 1; }
      TTree tout("tout", "Results of fit");
      // set branches
      tout.Branch("can", "TCanvas", &canva);
      double chi2 = 0; tout.Branch("chi2", &chi2);
      uint iStep = 0; tout.Branch("fit_step", &iStep);
      double pars_mirrow[MParKeeper::gI()->nPars()];
      for (uint i=0; i < MParKeeper::gI()->nPars(); i++)
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
            for (auto & pr : vpr) pr->RecalculateNextTime();
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
        for (iStep = 0; iStep < count; iStep++) {
          const libconfig::Setting &fit_step = strategy[iStep];
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

          // Plot all
          km.RecalculateNextTime();
          for (auto & pr : vpr) pr->RecalculateNextTime();
          for (uint i=0; i < Nrels; i++) {
            const DP & data = MRelationHolder::gI()->GetRelation(i).data;
            std::function<double(double)> func = MRelationHolder::gI()->GetRelation(i).func;
            // draw
            canva->cd(i+1)->Clear();
            draw(data)->Draw("ap");
            SET1(draw(func,
                      (data.data.begin())->x, (--data.data.end())->x, 200),
                 SetLineColor(kOrange) )->Draw("l");
            if (MRelationHolder::gI()->relationStatus(i))
              SET1(draw(func,
                        data.lrange, data.rrange, 200),
                   SetLineColor(kRed))->Draw("l");
          }
          MParKeeper::gI()->printAll();
          canva->SaveAs(TString::Format("%s/att%03d.step%d.pdf", dout_name.c_str(), e, iStep));

          // Fill result to tree
          // to copy to array from where it is copied to tree
          memcpy(pars_mirrow,
                 MParKeeper::gI()->get().data(),
                 sizeof(double)*MParKeeper::gI()->nPars());
          chi2 = MRelationHolder::gI()->CalculateChi2();
          tout.Fill();
        }
        delete min;
      }

      tout.Write();
      fout->Close();
    }
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"fit\" secton" << std::endl;
    return EXIT_FAILURE;
  }

  delete canva;

  // well done

  return(EXIT_SUCCESS);
}
