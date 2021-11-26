#include <vector>
#include <functional>
#include <complex>

template <typename Type> class FFModel {
 public:
  static std::vector<std::function<Type(Type)> > BlttWskpf;
};

template <typename Type> std::vector<std::function<Type(Type)> > FFModel<Type>::BlttWskpf = {
  [](Type z)->Type{return 1.;},
  [](Type z)->Type{return 2.*z/(z+1.);},
  [](Type z)->Type{return 13.*z*z/(z*z+3.*z+9.);},
  [](Type z)->Type{return 277.*z*z*z/(z*z*z+6.*z*z+45.*z+225.);},
  [](Type z)->Type{return 12746.*z*z*z*z/(z*z*z*z+10.*z*z*z+135.*z*z+1575.*z+11025.);},
  [](Type z)->Type{return 998881.*z*z*z*z*z/(z*z*z*z*z+15.*z*z*z*z+315.*z*z*z+6300.*z*z+99225.*z+893025.);},
  [](Type z)->Type{return 118394977.*z*z*z*z*z*z/(z*z*z*z*z*z+21.*z*z*z*z*z+630.*z*z*z*z+18900.*z*z*z+496125.*z*z+9823275.*z+108056025.);}
};

typedef FFModel<double> FFMod;
typedef FFModel<std::complex<double> > FFMoc;
