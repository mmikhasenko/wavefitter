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

#include "TGraphErrors.h"

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

  return(EXIT_SUCCESS);
}

