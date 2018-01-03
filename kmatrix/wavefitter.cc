// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <climits>
#include "libconfig.h++"

#include "MIsobar.h"
#include "MIsobarPiPiS.h"
#include "MIsobarChannel.h"
#include "MTwoBodyChannel.h"
#include "MParKeeper.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MRelationHolder.h"
#include "MAscoli.h"

#include "MatrixInverse.h"

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TColor.h"
#include "TH2D.h"
#include "TText.h"
#include "TMultiGraph.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/GSLMinimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "mstructures.h"

#define E_BEAM_LAB 190

void save_config_to_file(const char *config) {
  std::ifstream fconf(config);
  if (!fconf.is_open()) {
    std::cerr << "Warning<save_config_to_file> : Can not open config file. "
              << "Check your program, the error is very strange!\n";
    return;
  }

  // read text from config file
  TString TSconf("\n");
  std::string line;
  while (!fconf.eof()) {
    std::getline(fconf, line);
    TSconf += line;
    TSconf += "\n";
  }
  fconf.close();

  // save text to file
  gDirectory->mkdir("config");
  gDirectory->cd("config");
  TText TTconf(0., 0., TSconf);
  TTconf.SetName(config);
  TTconf.Write();
  gDirectory->cd("../");
}

int main(int argc, char *argv[]) {
  const char * iconfig = argv[1];
  libconfig::Config cfg;
  // Read the file. If there is an error, report it and exit.
  try {
    if (argc > 1) cfg.readFile(iconfig);
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

  const libconfig::Setting& root = cfg.getRoot();

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////// d a t a   s p e c i f i c a t i o n ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  std::vector<DP> whole_data;
  // Get data fiels
  if (root.exists("data")) {
    const libconfig::Setting &data = root["data"];

    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cout << "/////////////// Data specification: //////////////////\n";

    std::string path;
    data.lookupValue("path", path);
    std::cout << path << "\n";

    if (!data.exists("points")) {
      std::cerr << "Error <data>: no points section found! It does not make sence "
                << " that data section exists while no data files are provided, does it?\n";
      return EXIT_FAILURE;
    }
    const libconfig::Setting &dfiles = data["points"];
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
  } else {
    std::cerr << "\nWarning <main>: \"data\" section has not been found!\n\n";
  }

  const uint Nhist = whole_data.size();
  TCanvas c1("c1");
  c1.DivideSquare(Nhist);
  for (uint i = 0; i < Nhist; i++) {
    c1.cd(i+1); draw(whole_data[i])->Draw("ap");
  }
  c1.SaveAs("/tmp/default.data.pdf");
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
  // MIsobar  pipiS_iso(0.5, 0.5,  PI_MASS, PI_MASS, 0);
  MIsobarPiPiS  pipiS_iso;
  // model content is vector to fill K-matrix

  // K-matrix, just a reference
  const uint Jsector = 2;
  MmatrixK *km = 0;

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
      double R = 5.0;

      if (  !(channel.lookupValue("type", type) &&
              channel.lookupValue("L", L) &&
              channel.lookupValue("size", R)))
        continue;

      std::cout << "READ: " << std::left << type << ": "
                << (std::vector<char>{'S', 'P', 'D', 'F', 'G', 'H'})[L] << "-wave\n";

      const libconfig::Setting &spart = channel["particles"];

      if (type == "quasi-two-body") {
        // for sure, it should be strings
        std::string particles[2] = {spart[0], spart[1]};
        std::cout << "READ: particle species (" << particles[0] << ", " << particles[1] <<  ")\n";

        if (particles[0] == "rho") {
          MIsobarChannel *mCh = new MIsobarChannel(rho_iso, PI_MASS, L, R);
          mCh->makeLookupTable(mCh->sth(), 10., 100);
          iset.push_back(mCh);
        } else if (particles[0] == "f2") {
          MIsobarChannel *mCh = new MIsobarChannel(f2_iso, PI_MASS, L, R);
          mCh->makeLookupTable(mCh->sth(), 10., 100);
          iset.push_back(mCh);
        } else if (particles[0] == "pipiS") {
          MIsobarChannel *mCh = new MIsobarChannel(pipiS_iso, PI_MASS, L, R);
          mCh->makeLookupTable(mCh->sth(), 10., 100);
          iset.push_back(mCh);
        } else {
          std::cerr << "Error: isobar is not rho/f2/pipiS. Only them are available.\n";
          return EXIT_FAILURE;
        }
      } else if (type == "two-body") {
        double fpartmass[2] = {PI_MASS, PI_MASS};
        if (spart.getLength() != 2) {
          std::cerr << "Error<particle description>: there should be exactly 2 particles specified!\n";
          return EXIT_FAILURE;
        }
        for (int ps = 0; ps < spart.getLength(); ps++) {
          if (spart[ps].getType() == libconfig::Setting::TypeFloat) {
            fpartmass[ps] = spart[ps];
          } else if (spart[ps].getType() == libconfig::Setting::TypeString) {
            std::string pname = spart[ps];
            if (pname == "pi") {           fpartmass[ps] = PI_MASS;
            } else if (pname == "rho") {   fpartmass[ps] = RHO_MASS;
            } else if (pname == "f2") {    fpartmass[ps] = F2_MASS;
            } else if (pname == "pipiS") { fpartmass[ps] = PIPIS_MASS;
            } else {  std::cerr << "Error<particle description>: particle must be mass or name (" << pname << ")\n"; return EXIT_FAILURE; }
          } else {
            std::cerr << "Error: isobar is not pi/rho/f2/pipiS. Only them are available.\n";
            return EXIT_FAILURE;
          }
        }
        std::cout << "READ: particle masses (" << fpartmass[0] << ", " << fpartmass[1] <<  ")\n";
        MTwoBodyChannel *mCh = new MTwoBodyChannel(fpartmass[0], fpartmass[1], L, R);
        iset.push_back(mCh);
      } else {
        std::cerr << "Error: type of model channel can be only two-body of quasi-two-body";
        return EXIT_FAILURE;
      }
    }
    for (auto && it : iset) it->makeDisperseLookupTable(0.01, 10., 100);

    // create k-matrix;
    km = new MmatrixK(iset, 0);

    // content of the model
    const libconfig::Setting &content = modelT["content"];
    count = content.getLength();
    for (uint i = 0; i < count; ++i) {
      const libconfig::Setting &kcomp = content[i];

      std::string type = "?",
        masssq = "?",
        couplings = "?";

      if (  !(kcomp.lookupValue("type", type) &&
              kcomp.lookupValue("masssq", masssq) &&
              kcomp.lookupValue("couplings", couplings)))
        continue;

      std::cout << "READ: " << type << " ("
                << masssq << ", " << couplings << "*)"
                << std::endl;

      if (type == "pole") {
        km->addPole(masssq, couplings);
      } else if (type == "pole-like-background") {
        km->addBackground(masssq, couplings);
      } else {
        std::cerr << "Error: model blok is not a pole or pole-like-background. Only poles and bgds are available.";
        return EXIT_FAILURE;
      }
    }  // icount
    // create the model
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"model\" section!" << std::endl;
    return EXIT_FAILURE;
  }


  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////// Production model: /////////////////////////////////////////////////////
  /////////////////////// short range, long range, unitarisation ////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  std::vector<MProductionPhysics*> vpr;  // (iset);
  // std::vector<std::vector<std::pair<double, double> > > long_range_lookup_values[iset.size()];
  std::vector<std::vector<std::vector<std::pair<double, double> > > > long_range_lookup_values(iset.size());
  // channels<     modelA<     points< m, value > > >
  try {
    const libconfig::Setting &modelsA = root["modelA"];
    const uint Nmodels = modelsA.getLength();
    for (uint i = 0; i < iset.size(); i++) long_range_lookup_values[i].resize(Nmodels);

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
        pr.addScattering([&](double s)->b::matrix<cd>{return km->getValue(s);});

      // short_range
      if (modelA.exists("short_range")) {
        const libconfig::Setting &short_range = modelA["short_range"];
        std::string type = short_range["type"];

        if (type == "pietarinen") {
          std::string rhc_name = short_range["rhc"]; uint irhc = MParKeeper::gI()->add(rhc_name);
          std::string slope_name = short_range["slope"]; uint islope = MParKeeper::gI()->add(slope_name);
          // loop over powers
          std::vector<std::string> powers;
          const libconfig::Setting &orders = short_range["powers"];
          const uint count = orders.getLength();
          std::cout << "READ: powers are ";
          for (uint i = 0; i < count; i++) {
            std::string name = orders[i];
            std::cout << name << ((i < count-1) ? ", " : "\n");
            powers.push_back(name);
          }
          pr.addShortRange(powers,
                           [&, irhc, islope](double s)->cd{
                             double rhc = MParKeeper::gI()->get(irhc);
                             double slope = MParKeeper::gI()->get(islope);
                             cd r = (s-rhc)/slope;
                             return (1.-sqrt(r)) / (1.+sqrt(r)); });

        } else if (type == "polinomial") {
          // loop over powers
          std::vector<std::string> powers;
          const libconfig::Setting &orders = short_range["powers"];
          const uint count = orders.getLength();
          std::cout << "READ: powers are ";
          for (uint i = 0; i < count; i++) {
            std::string name = orders[i];
            powers.push_back(name);
            std::cout << name << ((i < count-1) ? ", " : "\n");
          }
          pr.addShortRange(powers, [](double s)->cd{return s;});

        } else {
          std::cerr<< "Warning<addShortRange()>: started without arguments. "
                   << "Do you know what happens then?" << std::endl;
          pr.addShortRange();  /* to be improved */
        }
        if (short_range.exists("phase_lock")) {
          pr.lockShortRangePhase();
          std::cout << "PHASE IS LOCKED: Phase of short range production term will be locked. "
                    << "For every channel, imaginary part of the first term in the expansion "
                    << " determines now imaginary part of all futher terms.\n";
        }
      }

      if (modelA.exists("long_range")) {
        const libconfig::Setting &long_range = modelA["long_range"];
        std::string type;

        // check if title is specified
        if (!long_range.lookupValue("type", type)) {
          std::cerr << "Error<main,production>: long_range \"type\" is missing!\n";
          return EXIT_FAILURE;
        }
        // cases for different types
        std::vector<std::function<cd(double)> > getB(iset.size());
        if (type == "DeckJMS") {
          std::cerr << "Error<modelA setting> : DeckJMS is not supported anymore!\n";
          return 1.;
        } else if (type == "DeckAJ") {
          uint Sp; int M; double tP, R;
          if ( !( long_range.lookupValue("damping_R", R) &&
                  long_range.lookupValue("pomeron_virtuality", tP) &&
                  long_range.lookupValue("pomeron_S", Sp) &&
                  long_range.lookupValue("pomeron_M", M))
               ) {
            std::cerr << "Error<main,production>: some papapeters of Long range is missing!\n";
            return EXIT_FAILURE;
          }
          std::cout << "DeckAJ projections will be constructed with parameters: J = " << Jsector
                    << ", Sp = " << Sp << ", M = " << M << ", R = " << R << "\n";
          TCanvas c3("c3"); c3.DivideSquare(iset.size());
          for (uint i=0; i < iset.size(); i++) {
            MIsobarChannel *ich = dynamic_cast<MIsobarChannel*>(iset[i]);
            const MIsobar &iso = ich->getIsobar();
            double mR = iso.GetM();

            double from = POW2(3*PI_MASS);
            double to = POW2(2.5);
            uint Npoints = 100;
            long_range_lookup_values[i][imodelA].resize(Npoints);
            for (uint t = 0; t < Npoints; t++) {
              double wsq = from + (to-from)/(Npoints-1)*t;
              double w = sqrt(wsq);
              double m23 = 2*PI_MASS+(mR-2*PI_MASS)*(1.-exp(-1./mR*(w-3*PI_MASS)));
              double mAsq = POW2(PI_MASS);  // pion beam
              double mBsq = POW2(PROT_MASS);  // prothon target
              double mDsq = POW2(PROT_MASS);  // elastic recoil
              double m1sq = POW2(PI_MASS);  // bachelor is pion
              double value_deck_AJ =
                MAscoli::getProjectedReducedDeck(Jsector, M, ich->GetL(),
                                                 m23*m23, iso.GetL(), R,
                                                 wsq, tP,
                                                 POW2(PI_MASS),
                                                 2*PROT_MASS*E_BEAM_LAB,
                                                 mAsq, mBsq, mDsq, m1sq);
              long_range_lookup_values[i][imodelA][t] = std::make_pair(wsq, value_deck_AJ);
            }
            const std::vector<std::pair<double, double> > *ltable = &(long_range_lookup_values[i][imodelA]);
            getB[i] = [&, ltable, mR, Jsector](double s)->cd {
              auto it = --(ltable->end());
              // std::cout << it->first << "\n";
              if (s >= it->first) return it->second * it->first / s;
              return getvalue(s, ltable->data(), ltable->size());
            };
            c3.cd(i+1);
            draw([&, i](double s)->double{return real(getB[i](s));}, 1.0, POW2(4.2))->Draw("al");
          }

          c3.SaveAs("/tmp/deckAJ.pdf");
          // finally add Long Range
          if (long_range.exists("par_name")) {
            std::string par_name = long_range["par_name"];
            pr.addLongRange(getB, par_name);
          } else pr.addLongRange(getB);


          // check and save in case
          if (long_range.exists("export_lookup_tables")) {

            // check if paths are present
            if (!long_range.exists("long_range_lookup")) {
              std::cerr << "Error<main,production> : No paths \"long_range_lookup\" to save lookup are provided!\n";
              return EXIT_FAILURE;
            }
            const libconfig::Setting &long_range_lookup = long_range["long_range_lookup"];
            // check amount of lookup tables
            if (long_range_lookup.getLength() != static_cast<int>(iset.size())) {
              std::cerr << "Error<main,production>: long_range_lookup.getLength() != iset.size()\n";
              return EXIT_FAILURE;
            }
            for (uint i=0; i < iset.size(); i++) {
              // save to file
              const std::string &lut_path = long_range_lookup[i];
              std::ofstream myfile(lut_path);
              if (myfile.is_open()) {
                const std::vector<std::pair<double, double> > &ltable = long_range_lookup_values[i][imodelA];
                for (auto && it : ltable) myfile << it.first << " "  << it.second << "\n";
                myfile.close();
              }
              std::cout << "lookup table is saved at " << lut_path << "\n";
            }
          } // if export is present

        } else if (type == "txt") {

          // check of paths are present
          if (!long_range.exists("long_range_lookup")) {
            std::cerr << "Error<main,production> : long range type is \"txt\" but paths \"long_range_lookup\" are not provided!\n";
            return EXIT_FAILURE;
          }
          const libconfig::Setting &long_range_lookup = long_range["long_range_lookup"];

          // check amount of lookup tables
          if (long_range_lookup.getLength() != static_cast<int>(iset.size())) {
            std::cerr << "Error<main,production>: long_range_lookup.getLength() != iset.size()\n";
            return EXIT_FAILURE;
          }
          std::cout << "DeckAJ projections will be constructed from txt lookup tables\n";

          TCanvas c3("c3"); c3.DivideSquare(iset.size());
          for (uint i=0; i < iset.size(); i++) {
            const std::string &lut_path = long_range_lookup[i];
            TGraph *igr = new TGraph(lut_path.c_str());
            if (igr && igr->GetN() <= 0) {
              std::cerr << "Error<main,production>: long_range_lookup_table is empry or not appropriate format\n";
              return EXIT_FAILURE;
            }
            uint Npoints = igr->GetN();
            long_range_lookup_values[i][imodelA].resize(Npoints);
            for (int j = 0; j < igr->GetN(); j++) long_range_lookup_values[i][imodelA][j] = std::make_pair(igr->GetX()[j], igr->GetY()[j]);
            std::cout << "-> Lookup table is filled, Np = " << Npoints << "\n";
            delete igr;

            const std::vector<std::pair<double, double> > *ltable = &(long_range_lookup_values[i][imodelA]);
            getB[i] = [&, ltable](double s)->cd {
              auto it = --(ltable->end());
              // std::cout << it->first << "\n";
              if (s >= it->first) return it->second * it->first / s;
              return getvalue(s, ltable->data(), ltable->size());
            };
            c3.cd(i+1);
            draw([&, i](double s)->double{return real(getB[i](s));}, 1.0, POW2(4.2))->Draw("al");
          }

          c3.SaveAs("/tmp/deckTXT.pdf");

          // finally add Long Range
          std::string par_name = "B";
          if (long_range.exists("par_name")) long_range.lookupValue("par_name", par_name);
          if (long_range.exists("separate_couplings_per_channels")) {
            pr.addLongRangeSeparated(getB, par_name);
            std::cout << "READ: \"separate_couplings_per_channels\" will be used!\n";
          } else {
            pr.addLongRange(getB, par_name);
            std::cout << "DEFAULT: a single coupling for the background will be used!\n";
          }
        } else {
          std::cerr << "Error: long range interaction is not 'DeckJMS' neither 'DeckAscoli', but no other options are available!";
          return 0;
        }
        if (modelA.exists("unitarisation")) pr.unitarize();
      }
    }
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"production model\" section!" << std::endl;
    return EXIT_FAILURE;
  }
  MParKeeper::gI()->printAll();

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// r e l a t i o n ///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  if (root.exists("adjustment")) {
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
        MChannel *ciso = iset[iCh];
        std::function<double(double)> int_lambda_function;
        if (iRel.getLength() <= 3) {
          int_lambda_function = [&, iCh, ciso, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return norm(v(iCh))*ciso->rho(e*e) * e;
          };
        } else {
          std::cout << "READ: fourth parameter to scale intensity: ";
          if (iRel[3].getType() == libconfig::Setting::TypeFloat) {
            double cbFactor = iRel[3];
            std::cout << "cbFactor = " << cbFactor << "\n";
            int_lambda_function = [&, iCh, ciso, pr, cbFactor](double e)->double{
              auto v = pr->getValue(e*e);
              return cbFactor*norm(v(iCh))*ciso->rho(e*e) * e;
            };
          } else if (iRel[3].getType() == libconfig::Setting::TypeString) {
            const std::string cbFactor_str = iRel[3];
            uint cbFactor_index = MParKeeper::gI()->add(cbFactor_str);
            std::cout << "cbFactor named as \"" << cbFactor_str << "\"("<< cbFactor_index << ")\n";
            int_lambda_function = [&, iCh, ciso, pr, cbFactor_index](double e)->double{
              auto v = pr->getValue(e*e);
              double cbFactor = MParKeeper::gI()->get(cbFactor_index);
              return cbFactor*norm(v(iCh))*ciso->rho(e*e) * e;
            };
          } else {
            std::cerr << "Error<main,relations>: fourth argument have to be float of string!\n";
            return 1;
          }
        }
        MRelationHolder::gI()->AddRelation(whole_data[jData], int_lambda_function);
      } else if (type == "DECAY_dGdM@") {
        if (iRel.getLength() != 4 && iRel.getLength() != 5) {
           std::cerr << "Error<main,relations>: Expectation does not match a number of arguments in DECAY_dGdM!!\n";
        }
        if (iRel.getLength() == 5) {
          uint iModel = iRel[2][0];
          uint iCh = iRel[2][1];
          MProductionPhysics *pr = vpr[iModel];
          MChannel *ciso = iset[iCh];
          // get decay mass A -> SYSTEM + B
          double massA = iRel[3], massB = iRel[4];
          // std::cout << "Constract differential rate dGdM for the A(" << massA << ")->(SYSTEM)+B(" << massB << ")\n";
          std::function<double(double)> int_lambda_function;
          int_lambda_function = [&, iCh, ciso, pr, massA, massB](double e)->double{
            if (e+massB > massA) return 0.0;
            auto v = pr->getValue(e*e);
            return norm(v(iCh))*ciso->rho(e*e)*sqrt(LAMBDA(massA*massA,e*e,massB*massB)) * e;  // x E because of jacobian ds = 2e de
          };
          MRelationHolder::gI()->AddRelation(whole_data[jData], int_lambda_function);
        } else if (iRel.getLength() == 4) {
          uint iModel = 0;
          MProductionPhysics *pr = vpr[iModel];
          // get decay mass A -> SYSTEM + B
          double massA = iRel[2], massB = iRel[3];
          std::function<double(double)> int_lambda_function;
          int_lambda_function = [&, pr, massA, massB](double e)->double{
            if (e+massB > massA) return 0.0;
            double intensity_sum = 0;
            auto v = pr->getValue(e*e);
            for (uint iCh=0; iCh < v.size(); iCh++) intensity_sum += norm(v(iCh))*iset[iCh]->rho(e*e);
            return intensity_sum*sqrt(LAMBDA(massA*massA,e*e,massB*massB)) * e;  // x E because of jacobian ds = 2e de
          };
          MRelationHolder::gI()->AddRelation(whole_data[jData], int_lambda_function);
        }
      } else if (type == "Re@") {
        uint iModel = iRel[2][0];
        uint iCh = iRel[2][1];
        MProductionPhysics *pr = vpr[iModel];
        MChannel *ciso = iset[iCh];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh, ciso, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return real(v(iCh));
          });
      } else if (type == "Im@") {
        uint iModel = iRel[2][0];
        uint iCh = iRel[2][1];
        MProductionPhysics *pr = vpr[iModel];
        MChannel *ciso = iset[iCh];
        MRelationHolder::gI()->AddRelation(whole_data[jData], [&, iCh, ciso, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return imag(v(iCh));
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

        std::function<double(double)> int_lambda_function;
        if (iRel.getLength() <= 3) {
          int_lambda_function = [&, iCh0, iCh1, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return real(v(iCh0)*conj(v(iCh1)))*sqrt(iset[iCh0]->rho(e*e)*iset[iCh1]->rho(e*e));
          };
        } else {
          std::cout << "READ: fourth parameter to scale intensity: ";
          if (iRel[3].getType() == libconfig::Setting::TypeFloat) {
            double cbFactor = iRel[3];
            std::cout << "cbFactor = " << cbFactor << "\n";
            int_lambda_function = [&, iCh0, iCh1, pr, cbFactor](double e)->double{
              auto v = pr->getValue(e*e);
              return cbFactor * real(v(iCh0)*conj(v(iCh1)))*sqrt(iset[iCh0]->rho(e*e)*iset[iCh1]->rho(e*e));
            };
          } else {
            std::cerr << "Error<main,relations>: fourth argument have to be float!\n";
            return 1;
          }
        }

        MRelationHolder::gI()->AddRelation(whole_data[jData], int_lambda_function);

      } else if (type == "ImInterf@") {
        uint iModel0 = iRel[2][0][0]; uint iCh0 = iRel[2][0][1];
        uint iModel1 = iRel[2][1][0]; uint iCh1 = iRel[2][1][1];
        if (iModel0 != iModel1) { std::cerr << "Error: iModel0!=iModel1!\n"; return 1;}
        MProductionPhysics *pr = vpr[iModel1];
        std::function<double(double)> int_lambda_function;
        if (iRel.getLength() <= 3) {
          int_lambda_function = [&, iCh0, iCh1, pr](double e)->double{
            auto v = pr->getValue(e*e);
            return imag(v(iCh0)*conj(v(iCh1)))*sqrt(iset[iCh0]->rho(e*e)*iset[iCh1]->rho(e*e));
          };
        } else {
          std::cout << "READ: fourth parameter to scale intensity: ";
          if (iRel[3].getType() == libconfig::Setting::TypeFloat) {
            double cbFactor = iRel[3];
            std::cout << "cbFactor = " << cbFactor << "\n";
            int_lambda_function = [&, iCh0, iCh1, pr, cbFactor](double e)->double{
              auto v = pr->getValue(e*e);
              return cbFactor * imag(v(iCh0)*conj(v(iCh1)))*sqrt(iset[iCh0]->rho(e*e)*iset[iCh1]->rho(e*e));
            };
          } else {
            std::cerr << "Error<main,relations>: fourth argument have to be float!\n";
            return 1;
          }
        }
        MRelationHolder::gI()->AddRelation(whole_data[jData], int_lambda_function);
      }
    }
  } else {
    std::cerr << "\nWarning <main>: \"adjustment\" section has not been found! No relations are set.\n\n";
  }

  const uint Nrels = MRelationHolder::gI()->Nrels();
  MRelationHolder::gI()->Print();

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
  //////////////////////////// p l o t  s e t t i n g s /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  try {
    if (root.exists("plot_settings")) {
      const libconfig::Setting &plot_settings = root["plot_settings"];

      std::cout << "\n\n";
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      std::cout << "/////////////// Plot settings: ///////////////////////\n";

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
	      if (tres->GetListOfBranches()->FindObject(name.c_str())) {
		// std::cout << "FOUND brahch with name " << name << "\n";
		tres->SetBranchAddress(name.c_str(), &pars[i]);
	      } else {
		std::cout << "NOT FOUND brahch with name " << name << "\n";
	      }
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

      // draw data
      std::string type;
      if (!plot_settings.lookupValue("type", type)) type = "model";


      uint Npages = 0;
      if (plot_settings.exists("mapping")) {
        Npages = 1;
        std::cout << "A parameter \"mapping\" is found and one-page pdf will be produced!\n";
      } else if (plot_settings.exists("mapping_list")) {
        Npages = plot_settings["mapping_list"].getLength();
        std::cout << "A parameter \"mapping_list\" is found and " << Npages << "-pages pdf will be produced!";
      } else {
        std::cerr << "Error<main,plot>: mapping or mapping_list has to be specified!";
        return EXIT_FAILURE;
      }

      TCanvas *canva = new TCanvas("canva", "title");
      if (plot_settings.exists("canva_settings")) {
        const libconfig::Setting &canva_sets = plot_settings["canva_settings"];
        delete canva;
        double canva_width, canva_height;
        if ( !canva_sets.lookupValue("width", canva_width) ||
             !canva_sets.lookupValue("height", canva_height) ) {
          std::cerr << "Error<main,plot>:\"canva_settings\" is present but does not have \"width\" or \"heigth\"!\n";
          return EXIT_FAILURE;
        }
        canva = new TCanvas("canva", "title", 0., 0.,
                            canva_width, canva_height);
      }

      // name of pdf output
      std::string fplot_name = "/tmp/default_plot.read_model_settings.pdf";
      if (!plot_settings.lookupValue("fplot_name", fplot_name))
	std::cerr << "Warning: fplot_name is not specified. A default name wil be used.\n";

      // open file for canvas
      TFile fcanvas(TString::Format("%s.root", fplot_name.c_str()), "RECREATE");
      save_config_to_file(iconfig);

      // loop ove pages
      for (uint pg = 0; pg < Npages; pg++) {
        const libconfig::Setting &mapping = (plot_settings.exists("mapping")) ? plot_settings["mapping"] : plot_settings["mapping_list"][pg];
        const uint NrelsToPlot = mapping.getLength();
        canva->Clear();

        /**************  TEMPERARY FIX ******************/
	uint npX = 1, npY = 1;
	if(!plot_settings.lookupValue("nplotsX", npX) ||
	   !plot_settings.lookupValue("nplotsY", npY)) {
	  canva->DivideSquare(NrelsToPlot);
	} else {
	  std::cout << "READ: canva for plotting will be divided " << npX << " x " << npY << "\n";
	  canva->Divide(npX, npY);
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

        // add model
        if ( type.find("model") != std::string::npos ) {
          if (plot_settings.exists("what_to_plot")) {
            const libconfig::Setting &what_to_plot = plot_settings["what_to_plot"];
            const uint Ncurves = what_to_plot.getLength();
            std::cout << "--> " << Ncurves << " model settings are found to be plotted.\n";
            for (uint j = 0; j < Ncurves; j++) {
              std::string title = "Default"; if (what_to_plot[j].exists("title")) what_to_plot[j].lookupValue("title", title);
              uint color = 0;  if (what_to_plot[j].exists("color")) what_to_plot[j].lookupValue("color", color);
              std::vector<std::pair<uint, double> > parsBackUp;
              if (what_to_plot[j].exists("set_to_value")) {
                const libconfig::Setting &setv = what_to_plot[j]["set_to_value"];
                const uint Nv = setv.getLength();
                for (uint i = 0; i < Nv; i++) {
                  std::string pname;
                  double value = 0.0;
                  if (setv[i].getType() == libconfig::Setting::TypeString) {
                    const std::string pname0 = setv[i];
                    pname = pname0;
                  } else if (setv[i].getType() == libconfig::Setting::TypeList) {
                    if (setv[i].getLength() == 2) {
                      const std::string pname0 = setv[i][0]; pname = pname0;
                      value = setv[i][1];
                    } else if (setv[i].getLength() == 3) {
                      const std::string pname0 = setv[i][0]; pname = pname0;
                      const std::string rname = setv[i][2];
                      value = setv[i][1];
                      value *= MParKeeper::gI()->get(rname);
                    } else {
                      std::cerr << "Error<plot settings>: set_to_value arguments can be only in format\n"
                                << "    (\"name\",value) - should be clear"
                                << "    (\"name1\",value,\"name2\") - to set value*get(\"name2\") to par \"name1\"\n\n";
                    }
                  } else {
                    std::cerr << "Error<plot, settings>: set_to_value requires can be list of names or list of lists\n";
                    return EXIT_FAILURE;
                  }
                  // backup
                  uint ip = MParKeeper::gI()->getIndex(std::string(pname));
                  if (ip == MParKeeper::error_uint) {std::cerr << "Error<main,plot>: \"set_to_zero\" check parameter name!\n"; return EXIT_FAILURE;}
                  double vp = MParKeeper::gI()->get(ip);
                  parsBackUp.push_back(std::make_pair(ip, vp));
                  // set
                  MParKeeper::gI()->set(pname, value);
                  std::cout << "-------> set_to_value: \"" << pname << "\" is set to " << value << "\n";
                }

                // Important to recalculate later
                km->RecalculateNextTime();
                for (auto & pr : vpr) pr->RecalculateNextTime();
              }
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
              const uint Nv = parsBackUp.size();
              for (uint p = 0; p < Nv; p++)
                MParKeeper::gI()->set(parsBackUp[p].first, parsBackUp[p].second);
            }
          } else {
            std::cerr << "Error<main,plot>: \"what_to_plot\" have to be specified for the model!\n";
          }
        }

        // finally plot
        for (uint i=0; i < NrelsToPlot; i++) {
          if (mapping[i].getLength() == 0) continue;
          // setTitle
          const uint iPad = mapping[i][0];
          const DP & data = MRelationHolder::gI()->GetRelation(iPad).data;
          // draw
          TVirtualPad *pd = canva->cd(i+1);
          pd->SetRightMargin(0.);
          pd->SetTopMargin(0.);

          mgr[i]->SetTitle(data.title.c_str());
          // gStyle->SetTitleSize(0.0015);
          // gStyle->SetTitleOffset(-0.2);
          mgr[i]->Draw("a");
          mgr[i]->GetHistogram()->SetTitle(data.title.c_str());
          mgr[i]->GetHistogram()->SetTitleSize(0.0015);
          mgr[i]->GetHistogram()->SetTitleOffset(-0.2);

          mgr[i]->GetXaxis()->SetLabelSize(0.07);
          mgr[i]->GetYaxis()->SetLabelSize(0.05);
          mgr[i]->GetXaxis()->SetTitleSize(0.07);
          mgr[i]->GetXaxis()->SetTitleOffset(-0.35);
        }

        // save multipage pdf;
        if (Npages != 1 && pg == Npages-1) { canva->Print(TString::Format("%s)", fplot_name.c_str()), "pdf");
        } else if (Npages != 1 && pg == 0) { canva->Print(TString::Format("%s(", fplot_name.c_str()), "pdf");
        } else {
          std::cout << "-----> Saving page " << pg << "\n";
          canva->Print(TString::Format("%s", fplot_name.c_str()), "pdf");
	}  // if
	// save root pictures to file
	canva->SetName(TString::Format("page%d", pg));
	fcanvas.cd(); canva->Write();
      }  // Npages
      delete canva;
    }  // exists plot_settings
  }  // try
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"plot_settings\" section" << std::endl;
    return EXIT_FAILURE;
  }


  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// f i t  s e t t i n g s /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  try {
    if (root.exists("fit_settings")) {
      const libconfig::Setting &fit_settings = root["fit_settings"];

      std::cout << "\n\n";
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      std::cout << "/////////////// Fit settings: ////////////////////////\n";

      const libconfig::Setting &strategy = fit_settings["strategy"];

      const uint nAttempts = fit_settings["nAttempts"];

      int seed;
      uint pid = ::getpid();
      if (!fit_settings.lookupValue("seed", seed)) {
        std::srand(pid); const uint frand = std::rand();
        seed = frand+std::time(0);
      }
      std::srand(seed);

      const std::string &dout_name = fit_settings["dout_name"];


      TCanvas *canva = new TCanvas("canva", "title"); canva->DivideSquare(Nrels);
      /*********************************** Fit itself *****************************************/
      /****************************************************************************************/
      const int rand_file_id = std::rand()%1000;
      TString fout_name = TString::Format("%s/fit.results.root.pid%d.rand%03d", dout_name.c_str(), pid, rand_file_id);
      std::cout << "IMPORTANT: results of the fit are stored in the file " << fout_name.Data() << "\n";
      TFile *fout = new TFile(fout_name, "RECREATE");
      if (!fout) { std::cerr << "no fout acceptable!\n"; return EXIT_FAILURE; }
      save_config_to_file(iconfig);
      TTree tout("tout", "Results of fit");
      // set branches
      // tout.Branch("can", "TCanvas", &canva);
      double chi2 = 0; tout.Branch("chi2", &chi2);
      uint iStep = 0; tout.Branch("fit_step", &iStep);
      double status = 0; tout.Branch("status", &status);
      double edm = 0; tout.Branch("edm", &edm);
      uint eAtt; tout.Branch("eAtt", &eAtt);
      tout.Branch("pid", &pid);
      double pars_mirrow[MParKeeper::gI()->nPars()];
      for (uint i=0; i < MParKeeper::gI()->nPars(); i++)
        tout.Branch(MParKeeper::gI()->getName(i).c_str(), &pars_mirrow[i]);

      // length of chi2 relations arrays
      uint bins_in_relation[Nrels];
      tout.Branch("bins_in_relation", &bins_in_relation,
                  TString::Format("bins_in_relations[%d]/i", Nrels));
      // creare a branch for every relation
      double *chi2_relation[Nrels];
      for (uint i=0; i < Nrels; i++) {
        if (MRelationHolder::gI()->relationStatus(i)) {
          uint sizeR = MRelationHolder::gI()->ExtractChi2Array(i);
          chi2_relation[i] = new double[sizeR+1];
          tout.Branch(TString::Format("chi2_relation_%d", i), chi2_relation[i],
                      TString::Format("chi2_relation_%d[%d]/D", i, sizeR+1));
        }
      }
      // to copy to array from where it is copied to tree
      for (eAtt = 0; eAtt < nAttempts; eAtt++) {
        std::cout << "---------- Attempt " << eAtt << " -----------" << std::endl;
        std::cout << "------------------------------------------" << std::endl;

        /**************************************MINIMIZE******************************************/
        // Build minimizer
        ROOT::Math::Minimizer* min =
          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

        // set tolerance , etc...
        min->SetMaxFunctionCalls(UINT_MAX);
        min->SetMaxIterations(UINT_MAX);
        min->SetTolerance(0.001);
        min->SetStrategy(1);
        min->SetPrintLevel(2);
        // min->Options().Print();

        // Create funciton wrapper for minmizer a IMultiGenFunction type
        ROOT::Math::Functor functor([&](const double *pars)->double {
            MParKeeper::gI()->pset(pars);
            MParKeeper::gI()->satisfyConstraints();
            km->RecalculateNextTime();
            for (auto & pr : vpr) pr->RecalculateNextTime();
            return MRelationHolder::gI()->CalculateChi2();
          }, pnPars);
        min->SetFunction(functor);

        // set starting parameters
        // if possible then from the file
        if (fit_settings.exists("path_to_starting_value") &&
            fit_settings.exists("entry")) {
          const std::string spar_path = fit_settings["path_to_starting_value"];
          const uint entry = fit_settings["entry"];
          TFile *fres = TFile::Open(spar_path.c_str());
          if (!fres) {
            std::cerr << "File with results specified but not found!\n";
            return EXIT_FAILURE;
          }
          std::cout << "File with parameter values successfully opened!\n";
          TTree *tres; gDirectory->GetObject("tout", tres);
          if (!tres) {
            std::cerr << "Tree with starting values not found by name 'tout'!\n";
            return EXIT_FAILURE;
          }
          const uint Npars = MParKeeper::gI()->nPars();
          double pars[Npars];
          for (uint i = 0; i < Npars; i++) {
            const std::string & name = MParKeeper::gI()->getName(i);
            // check if it is at list of branches
	    if (tres->GetListOfBranches()->FindObject(name.c_str())) {
	      // std::cout << "FOUND brahch with name " << name << "\n";
	      tres->SetBranchAddress(name.c_str(), &pars[i]);
	    } else {
	      std::cout << "NOT FOUND brahch with name " << name << "\n";
	      pars[i] = MParKeeper::gI()->get(i);
	    }
	  }
          // set values from tree and set to keeper
          tres->GetEntry(entry);
          for (uint i = 0; i < Npars; i++) {
            const std::string & name = MParKeeper::gI()->getName(i);
	    if (tres->GetListOfBranches()->FindObject(name.c_str()))
              MParKeeper::gI()->set(i, pars[i]);
	  }
          fres->Close(); fout->cd();
        } else {
          // otherwise from random
          MParKeeper::gI()->randomizePool();
        }
        MParKeeper::gI()->printAll();
        for (uint i=0; i < pnPars; i++)
	  min->SetVariable(i,
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

          if (fit_step.exists("set_to_value")) {
            const libconfig::Setting &setv = fit_step["set_to_value"];
            const uint Nv = setv.getLength();
            for (uint i = 0; i < Nv; i++) {
              const std::string pname = setv[i][0];
              const double value = setv[i][1];
              MParKeeper::gI()->set(pname, value);
              min->SetVariableValue(min->VariableIndex(pname), value);
              std::cout << "-------> set_to_value: \"" << pname << "\" is set to " << value << "\n";
            }
          }
          if (fit_step.exists("randomize")) {
            const libconfig::Setting &setv = fit_step["randomize"];
            const uint Nv = setv.getLength();
            for (uint i = 0; i < Nv; i++) {
              const std::string pname = setv[i];
	      uint ind = MParKeeper::gI()->getIndex(pname);
              MParKeeper::gI()->setRandom(ind);
	      double value = MParKeeper::gI()->get(ind);
              min->SetVariableValue(min->VariableIndex(pname), value);
              std::cout << "-------> randomize: \"" << pname << "\" is set to " << value << "\n";
            }
          }
          MParKeeper::gI()->removeConstraints();
          if (fit_step.exists("constraints")) {
            const libconfig::Setting &setv = fit_step["constraints"];
            const uint Nv = setv.getLength();
            for (uint i = 0; i < Nv; i++) {
              const libconfig::Setting &constr = setv[i];
              if (constr.getLength() != 3) {
                std::cerr << "Error<main,fit,constraints>: only multiplicative constraints are implemented! "
                          <<" (\"par1\",factor,\"par2\") means \"par1\"=factor*\"par2\"\n";
                continue;
              }
              std::string i_pname = constr[0]; uint i_ind = MParKeeper::gI()->getIndex(i_pname);
              std::string j_pname = constr[2]; uint j_ind = MParKeeper::gI()->getIndex(j_pname);
              double factor = constr[1];
              MParKeeper::gI()->addConstraint(i_ind, j_ind, factor);
              std::cout << "-------> added constraint: \"" << i_pname << "\" is set to " << factor << " * " << j_pname << "\n";
            }
          }
          MParKeeper::gI()->printAll();
          // adjust which parameters to vary
          const uint nPars_to_vary = pars.getLength();
          // loop over all parameters for make fix them
          for (uint r=0; r < pnPars; r++) min->FixVariable(r);
          for (uint r=0; r < nPars_to_vary; r++) {
            const std::string &pname = pars[r];
            const int index = MParKeeper::gI()->pgetIndex(pname);
            if (index == MParKeeper::error_uint) {std::cerr << "Error<main,fit>: \"pars_to_vary\" check parameter name!\n"; return EXIT_FAILURE;}
            min->ReleaseVariable(index);
          }
          // minimize
          min->Minimize();
          status = min->Status();
	  edm = min->Edm();

	  MParKeeper::gI()->pset(min->X());
          MParKeeper::gI()->satisfyConstraints();  // dirty fix

          if (fit_settings.exists("save_preview")) {
            // Plot all
            km->RecalculateNextTime();
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
            std::string save_preview = fit_settings["save_preview"];
            if (save_preview != "no" && save_preview != "false")
              canva->SaveAs(TString::Format("%s/att%03d.step%d.pid%d.rand%03d.pdf", dout_name.c_str(), eAtt, iStep, pid, rand_file_id));
          } else {
            std::cerr << "Warning: save_preview option is not specified. \"no\" is set by default.";
          }
          // Fill result to tree
          // to copy to array from where it is copied to tree
          memcpy(pars_mirrow,
                 MParKeeper::gI()->get().data(),
                 sizeof(double)*MParKeeper::gI()->nPars());
          chi2 = MRelationHolder::gI()->CalculateChi2();
          // Fill relation blocks
           for (uint i=0; i < Nrels; i++) {
             if (MRelationHolder::gI()->relationStatus(i)) {
               uint sizeR = MRelationHolder::gI()->ExtractChi2Array(i, &chi2_relation[i][1]);
               chi2_relation[i][0] = 0.;
               for (uint p = 1; p <= sizeR; p++) chi2_relation[i][0] += chi2_relation[i][p];
               bins_in_relation[i] = sizeR;
             } else {
               uint sizeR = MRelationHolder::gI()->ExtractChi2Array(i);
               for (uint p = 0; p <= sizeR; p++) chi2_relation[i][p] = 0;
               bins_in_relation[i] = 0;
             }
           }
          tout.Fill();
        }
        delete min;
      }

      tout.Write();
      fout->Close();
      for (uint i=0; i < Nrels; i++) delete [] chi2_relation[i];
      delete canva;
    }
  }
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"fit\" section" << std::endl;
    return EXIT_FAILURE;
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  //////////////////// c o n t i n u a t i o n  s e t t i n g s /////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  try {
    if (root.exists("continuation_settings")) {
      const libconfig::Setting &continuation_settings = root["continuation_settings"];

      std::cout << "\n\n";
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      std::cout << "/////////////// Continuation settings: ///////////////////////\n";

      // if file is specified
      //  - open file
      //  - find tree thee
      //  - load perameters from particular entry of the tree
      if (continuation_settings.exists("path_to_fit_results") && continuation_settings.exists("entry_fit_result")) {
        std::string fres_name; continuation_settings.lookupValue("path_to_fit_results", fres_name);
        uint entry =  continuation_settings["entry_fit_result"];  // throw exception it something is wrong
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
            return EXIT_FAILURE;
          }
        } else {
          std::cerr << "File with results specified but not found!\n";
          return EXIT_FAILURE;
        }
      }

      if (continuation_settings.exists("set_to_value")) {
        const libconfig::Setting &setv = continuation_settings["set_to_value"];
        const uint Nv = setv.getLength();
        for (uint i = 0; i < Nv; i++) {
          std::string pname;
          double value = 0.0;
          if (setv[i].getType() == libconfig::Setting::TypeString) {
            const std::string pname0 = setv[i];
            pname = pname0;
            continue;
          } else if (setv[i].getType() == libconfig::Setting::TypeList) {
            if (setv[i].getLength() == 2) {
              const std::string pname0 = setv[i][0]; pname = pname0;
              value = setv[i][1];
              continue;
            } else if (setv[i].getLength() == 3) {
              const std::string pname0 = setv[i][0]; pname = pname0;
              const std::string rname = setv[i][2];
              value = setv[i][1];
              value *= MParKeeper::gI()->get(rname);
            } else {
              std::cerr << "Error<continue settings>: set_to_value arguments can be only in format\n"
                        << "    (\"name\",value) - should be clear"
                        << "    (\"name1\",value,\"name2\") - to set value*get(\"name2\") to par \"name1\"\n\n";
            }
          } else {
            std::cerr << "Error<continue settings>: set_to_value requires can be list of names or list of lists\n";
            return EXIT_FAILURE;
          }
          MParKeeper::gI()->set(pname, value);
          std::cout << "-------> set_to_value: \"" << pname << "\" is set to " << value << "\n";
        }
      }
      MParKeeper::gI()->printAll();



      TCanvas canva_sheets("sheets");
      const libconfig::Setting &ranges = continuation_settings["plot_range"];
      const libconfig::Setting &real_range = ranges["real_range"];
      const libconfig::Setting &imag_range = ranges["imag_range"];
      if (real_range.getLength() != 3) {std::cerr << "Error<continuation_settings>: in plot_range, real_range "
                                                  << "is expected in the form [Nbins, left_value, right_value]\n"; return EXIT_FAILURE; }
      if (imag_range.getLength() != 3) {std::cerr << "Error<continuation_settings>: in plot_range, imag_range "
                                                  << "is expected in the form [Nbins, left_value, right_value]\n"; return EXIT_FAILURE; }
      const uint Nbx = real_range[0]; double lrx = real_range[1]; double rrx = real_range[2];
      const uint Nby = imag_range[0]; double lry = imag_range[1]; double rry = imag_range[2];
      TH2D hreal("realTm1","Real part of det[T_{I,II}^{-1}K];M=Re[#sqrt{s}];#Gamma=2Im[#sqrt{s}]", Nbx, lrx, rrx, Nby, lry, rry);
      TH2D himag("imagTm1","Imag part of det[T_{I,II}^{-1}K];M=Re[#sqrt{s}];#Gamma=2Im[#sqrt{s}]", Nbx, lrx, rrx, Nby, lry, rry);
      TH2D habs ( "absTm1", "Ln@Abs part of det[T_{I,II}^{-1}K];M=Re[#sqrt{s}];#Gamma=2Im[#sqrt{s}]", Nbx, lrx, rrx, Nby, lry, rry);

      // decide what to plot
      const std::string what_to_plot = continuation_settings["what_to_plot"];
      uint sheet = 0;
      if ( what_to_plot.find("first") != std::string::npos ) { sheet = 1; std::cout << "---> First sheet will be plotted!\n";}
      if ( what_to_plot.find("second") != std::string::npos ) { sheet = 2; std::cout << "---> Second sheet will be plotted!\n";}
      if ( what_to_plot.find("first" ) != std::string::npos &&
           what_to_plot.find("second") != std::string::npos) { sheet = 12; std::cout << "---> First and Second sheets will be plotted together!\n";}

      // calculation loop
      km->RecalculateNextTime();
      km->Print();
      std::cout << "\nCalculations for " << Nbx*Nby << "points started:\n";
      for (uint ix = 0; ix < Nbx; ix++) {
        for (uint iy = 0; iy < Nbx; iy++) {
          cd s(habs.GetXaxis()->GetBinCenter(ix+1), // s = (M+iG/2)^2;
	       habs.GetYaxis()->GetBinCenter(iy+1)/2.); s = s*s;
          cd Tm1 = 0;
	  if (sheet == 1)  Tm1 = 1./det_fast(km->getFSdenominator(s));
	  if (sheet == 2)  Tm1 = 1./det_fast(km->getSSdenominator(s));
	  if (sheet == 12) Tm1 = (imag(s) > 0) ? det_fast(km->getFSdenominator(s)) : det_fast(km->getSSdenominator(s));
          // if (sheet == 1)  Tm1 = 1./km->getValue(s)(0,0);
          // if (sheet == 2)  Tm1 = 1./km->getSSvalue(s)(0,0);
          // if (sheet == 12) Tm1 = (imag(s) > 0) ? 1./km->getValue(s)(0,0) : 1./km->getSSvalue(s)(0,0);
          if((ix*Nbx+iy) % 50 == 0) std::cout << std::setprecision(3) << 100.*(ix*Nby+iy)/(Nbx*Nby) << "%: s = " << s << ", Tm1 = " << Tm1 << "\n";
          hreal.SetBinContent(ix+1, iy+1, real(Tm1));
          himag.SetBinContent(ix+1, iy+1, imag(Tm1));
          habs .SetBinContent(ix+1, iy+1, log10(abs(Tm1)));
        }
      }
      std::cout << "\n";
      // save to pdf
      std::string fplot_name = "/tmp/default_plot.read_model_settings.pdf";
      if (!continuation_settings.lookupValue("fplot_name", fplot_name))
        std::cerr << "Warning: fplot_name is not specified. A default name wil be used.\n";

      /* TEMPERARY FIX */
      // Adjustion style
      const Int_t NRGBs = 5;
      const Int_t NCont = 70;
      Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
      Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
      Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
      Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
      TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
      gStyle->SetNumberContours(NCont);

      // save multipage pdf;
      habs .SetStats(kFALSE); habs .Draw("col"); // habs.Draw("cont3 same");
      TH2D *habs_less_contours = static_cast<TH2D*>(habs.Clone("absTm1_clone_less_contours"));
      if (continuation_settings.exists("ncontours")) {
        uint ncont = continuation_settings["ncontours"];
        habs_less_contours->SetContour(ncont);
      }
      habs_less_contours->SetLineWidth(1);
      habs_less_contours->SetLineColor(kBlack);
      habs_less_contours->Draw("cont3 same");
      canva_sheets.Print(TString::Format("%s", fplot_name.c_str()), "pdf");
      delete habs_less_contours;
      // hreal.SetStats(kFALSE); hreal.Draw("colz"); hreal.Draw("cont3 same"); canva_sheets.Print(TString::Format("%s" , fplot_name.c_str()), "pdf");
      // himag.SetStats(kFALSE); himag.Draw("colz"); himag.Draw("cont3 same"); canva_sheets.Print(TString::Format("%s)", fplot_name.c_str()), "pdf");
      TFile fout(TString::Format("%s.root", fplot_name.c_str()), "RECREATE");
      save_config_to_file(iconfig);
      habs .Write();
      hreal.Write();
      himag.Write();
      fout.Close();
    }  // exists continuation_settings
  }  // try
  catch(const libconfig::SettingNotFoundException &nfex) {
    std::cerr << "Error <> libconfig::SettingNotFoundException in \"continuation_settings\" section" << std::endl;
    return EXIT_FAILURE;
  }

  // well done

  return(EXIT_SUCCESS);
}
