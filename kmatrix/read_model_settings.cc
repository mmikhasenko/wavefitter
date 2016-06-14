// Copyright of libconfig++
// My modification of standard example

// [] Use to compile
// g++ -o read_model_settings \
//    /localhome/mikhasenko/Tools/libconfig-1.5/lib/libconfig++.so \
//    -I/localhome/mikhasenko/Tools/libconfig-1.5/include read_model_settings.cc

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

  std::vector<TGraphErrors*> vdata;
  // Get data fiels
  try {
    const libconfig::Setting &data = root["data"];
    std::string path;
    data.lookupValue("path", path);
    std::cout << path << "\n";

    const libconfig::Setting &dfiles = root["data"]["points"];
    int count = dfiles.getLength();

    for (int i = 0; i < count; ++i) {
      const libconfig::Setting &dfile = dfiles[i];

      // Only output the record if all of the expected fields are present.
      std::string title, type, file_name;
      int index;

      if (  !(dfile.lookupValue("type", type)
           && dfile.lookupValue("file_name", file_name)
           && dfile.lookupValue("title", title)
           && dfile.lookupValue("index", index)))
        continue;

      const libconfig::Setting &trust_range = dfile["trust_range"];
      const libconfig::Setting &fit_range = dfile["fit_range"];
      double v1 = fit_range[0];
      double v2 = fit_range[1];

      std::cout << "READ: " << std::left << file_name << "  "
                << title << "  "
                << index << "  ("
                << v1 << ", " << v2 << ")"
                << std::endl;

      TGraphErrors *g = (type == "txt") ? new TGraphErrors((path+file_name).c_str()) : 0;
      if (!g) {std::cerr << "Warning <> " << file_name << "not found!" << std::endl; continue; }
      std::string gname = path;
      gname.erase(std::remove(gname.begin(), gname.end(), '/'), gname.end());
      g->SetName(gname.c_str());
      g->SetTitle(title.c_str());
      vdata.push_back(g);
    }
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException" << std::endl;
    // Ignore.
  }

  // Prepare data and fuctions
  const uint Nhist = vdata.size();
  DP whole_data[Nhist];
  for (uint i = 0; i < Nhist; i++) {
    for (uint j = 0; j < vdata[i]->GetN(); j++) {
      whole_data[i].data.push_back(data_point{vdata[i]->GetX()[j],
            vdata[i]->GetY()[j],
            vdata[i]->GetEX()[j],
            });
      vdata[i]->GetEY()[j] = vdata[i]->GetEX()[j];
      vdata[i]->GetEX()[j] = 0;
    }
    whole_data[i].name = vdata[i]->GetName();
    whole_data[i].title = vdata[i]->GetTitle();
    whole_data[i].lrange = 0.;
    whole_data[i].rrange = 3.;
  }

  TCanvas c1("c1");
  c1.Divide(Nhist, 2);
  for (uint i = 0; i < Nhist; i++) {
    c1.cd(i+1); vdata[i]->Draw("alp");
    c1.cd(Nhist+i+1); draw(whole_data[i])->Draw("alp");
  }
  c1.SaveAs("/tmp/1.pdf");
  std::cout << "\n";

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// m o d e l   a d j u s t m e n t ////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  // channels
  std::vector<MIsobarChannel*> iset;
  // standard isobars
  MIsobar rho_iso(RHO_MASS, RHO_WIDTH, PI_MASS, PI_MASS, 1, 5.);
  MIsobar  f2_iso(F2_MASS, F2_WIDTH,  PI_MASS, PI_MASS, 2, 5.);
  // model content is vector to fill K-matrix
  std::vector<std::pair<std::string, std::string> > model_content;

  try {
    const libconfig::Setting &modelT = root["modelT"];

    const libconfig::Setting &channels = root["modelT"]["channels"];
    int count = channels.getLength();

    for (int i = 0; i < count; ++i) {
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
      if (particles[0] == "rho") { iset.push_back(new MIsobarChannel(rho_iso, PI_MASS, L));
      } else if (particles[0] == "f2") { iset.push_back(new MIsobarChannel(f2_iso, PI_MASS, L));
      } else {
        std::cerr << "Error: isobar is not rho/f2. Only them are available.";
        return EXIT_FAILURE;
      }
    }

    // content of the model
    const libconfig::Setting &content = root["modelT"]["content"];
    count = content.getLength();
    for (int i = 0; i < count; ++i) {
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
    // MmatrixK k(Nchannels, Npoles);
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException" << std::endl;
    // Ignore.
  }

  MmatrixK km(iset, model_content.size());
  MProductionPhysics pr(iset);
  pr.addScattering([&](double s)->b::matrix<cd>{return km.getValue(s);});
  pr.addShortRange();

  MParKeeper::gI()->printAll();

  // Relation holder
  // map comes from user
  std::vector<uint> assignment = {1, 3, 2};
  for (int i=0; i < assignment.size(); i++)
    MRelationHolder::gI()->AddRelation(whole_data[i], [&](double s)->double{
        auto v = pr.getValue(s); return norm(v(assignment[i]));
      });


  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// f i t  s e t t i n g s /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  // TODO(mikhasenko): make config reading for the fit settings
  std::vector<std::string> pars_input = {"m1",
                                         "g1", "g2", "g3",
                                         "m2",
                                         "h1", "h2", "h3",
                                         "c1", "c2", "c3"};
  const uint nPars = pars_input.size();
  if (nPars != MParKeeper::gI()->nPars()) {
    std::cerr << "Error: number of the specified parameters does not much to number of recerved parameters!";
    return EXIT_FAILURE;
  }

  const uint nAttempt = 10;
  for (uint e = 0; e < nAttempt; e++) {
    std::cout << "---------- Attempt " << e << " -----------"<< std::endl;
    std::cout << "------------------------------------------"<< std::endl;

    /**************************************MINIMIZE******************************************/
    /****************************************************************************************/
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
        MParKeeper::gI()->setPool(pars);
        return MRelationHolder::gI()->CalculateChi2();
      }, nPars);
    min->SetFunction(functor);
    // specify parameters
    for (int i=0; i < nPars; i++) min->SetVariable(i, pars_input[i].c_str(), 0, 0.1);
    // minimize
    min->Minimize();

    MParKeeper::gI()->printAll();

    delete min;
  }

  return(EXIT_SUCCESS);
}

