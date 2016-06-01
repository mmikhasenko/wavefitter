// Copyright [2016] Mikhail Mikhasenko

#ifndef __MSTRUCTURES__
#define __MSTRUCTURES__

#include <list>
#include <string>

struct data_points {
  std::list< std::array<double, 3> >;
  std::string name;
  std::string title;
};

#endif // __MSTRUCTURES__

