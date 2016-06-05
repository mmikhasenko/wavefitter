// Copyright [2016] Mikhail Mikhasenko

#include <iostream>
#include "MParKeeper.h"

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  int m1 = MParKeeper::gI()->add("m1", 1.1);
  int m2 = MParKeeper::gI()->add("m2", 2.1);
  int m3 = MParKeeper::gI()->add("m3", 1.1);
  int g1 = MParKeeper::gI()->add("g1", 8.8);
  int g2 = MParKeeper::gI()->add("g2", 4.3);

  cout << "m1 = " << MParKeeper::gI()->get(m1) << "\n"
       << "m2 = " << MParKeeper::gI()->get(m2) << "\n"
       << "m3 = " << MParKeeper::gI()->get(m3) << "\n"
       << "g1 = " << MParKeeper::gI()->get(g1) << "\n"
       << "g2 = " << MParKeeper::gI()->get(g2) << "\n";

  return 0;
}

