// Copyright [2016] Mikhail Mikhasenko

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH2D.h"
// #include "Math/SpecFuncMathMore.h"

// #include "MIsobar.h"
// #include "MAscoli.h"
// #include "MmatrixK.h"
// #include "MProductionPhysics.h"
// #include "MParKeeper.h"
#include "mstructures.h"
// #include "mintegrate.h"
#include "deflib.h"
#include "constants.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "mathUtils.hpp"

namespace b = boost::numeric::ublas;

cd ampAMP(double s) {

  /* digitize table */
  const double f[2] = {0.1968, -0.0154};  // AMP Table 1, M solution: f_1^1 and f_2^1
  std::vector<b::matrix<cd> > _a(2, b::matrix<cd>(2,2) );
  
  _a[0](0, 0) =  0.1131;  // AMP Table 1, M solution: f_2^2
  _a[0](0, 1) =  0.0150;  // AMP Table 1, M solution: f_1^3
  _a[0](1, 0) =  0.0150;  // AMP Table 1, M solution: f_1^3
  _a[0](1, 1) = -0.3216;  // AMP Table 1, M solution: f_2^3
  _a[1](0, 0) = f[0] * f[0];
  _a[1](0, 1) = f[0] * f[1];
  _a[1](1, 0) = f[1] * f[0];
  _a[1](1, 1) = f[1] * f[1];
  
  std::vector<b::matrix<cd> > _c(5, b::matrix<cd>(2,2) );
  _c[0](0, 0) =  0.0337;                // AMP Table 1, M solution: c_11^0
  _c[1](0, 0) = -0.3185;                // AMP Table 1, M solution: c_11^1
  _c[2](0, 0) = -0.0942;                // AMP Table 1, M solution: c_11^2
  _c[3](0, 0) = -0.5927;                // AMP Table 1, M solution: c_11^3
  _c[4](0, 0) =  0.1957;                // AMP Table 1, M solution: c_11^4
  _c[0](0, 1) = _c[0](1, 0) = -0.2826;  // AMP Table 1, M solution: c_12^0
  _c[1](0, 1) = _c[1](1, 0) =  0.0918;  // AMP Table 1, M solution: c_12^1
  _c[2](0, 1) = _c[2](1, 0) =  0.1669;  // AMP Table 1, M solution: c_12^2
  _c[3](0, 1) = _c[3](1, 0) = -0.2082;  // AMP Table 1, M solution: c_12^3
  _c[4](0, 1) = _c[4](1, 0) = -0.1386;  // AMP Table 1, M solution: c_12^4
  _c[0](1, 1) =  0.3010;                // AMP Table 1, M solution: c_22^0
  _c[1](1, 1) = -0.5140;                // AMP Table 1, M solution: c_22^1
  _c[2](1, 1) =  0.1176;                // AMP Table 1, M solution: c_22^2
  _c[3](1, 1) =  0.5204;                // AMP Table 1, M solution: c_22^3
  _c[4](1, 1) = -0.3977;                // AMP Table 1, M solution: c_22^4
  
  b::matrix<double> _sP(1, 2);
  _sP(0, 0) = -0.0074;  // AMP Table 1, M solution: s_0
  _sP(0, 1) =  0.9828;  // AMP Table 1, M solution: s_1
  

	// change parameters according to Kachaev's prescription
	_c[4](0, 0) = 0; // was 0.1957;
	_c[4](1, 1) = 0; // was -0.3977;

	_a[0](0, 1) = 0; // was 0.0150
	_a[0](1, 0) = 0; // was 0.0150

	// _a[1] are the f's from the AMP paper
	_a[1](0, 0) = 0;
	_a[1](0, 1) = 0;
	_a[1](1, 0) = 0;
	_a[1](1, 1) = 0;


  /* ************** */
  
  double mass = sqrt(s);
  // clarify 
  if (fabs(s - _sP(0, 1)) < 1e-6) s  = POW2(sqrt(s)+1e-6);

  double _piChargedMass = PI_MASS;
  double _piNeutralMass = PI0_MASS;
  double _kaonChargedMass = K_MASS;
  double _kaonNeutralMass = K0_MASS;
  double _kaonMeanMass    = (_kaonChargedMass + _kaonNeutralMass) / 2;
  
  const cd qPiPi   = sqrt( LAMBDA(cd(s,0.), POW2(_piChargedMass), POW2(_piChargedMass))/(4.*s) );
  const cd qPi0Pi0 = sqrt( LAMBDA(cd(s,0.), POW2(_piNeutralMass), POW2(_piNeutralMass))/(4.*s) );
  const cd qKK     = sqrt( LAMBDA(cd(s,0.), POW2(_kaonChargedMass), POW2(_kaonChargedMass))/(4.*s) );
  const cd qK0K0   = sqrt( LAMBDA(cd(s,0.), POW2(_kaonNeutralMass), POW2(_kaonNeutralMass))/(4.*s) );
  // cd qKmKm         = LAMBDA(cd(s,0.), POW2(_kaonMeanMass), POW2(_kaonMeanMass));

  b::matrix<cd> rho(2, 2);
  rho(0, 0) = ((2. * qPiPi) / mass + (2. * qPi0Pi0) / mass) / 2.;
  rho(1, 1) = ((2. * qKK)   / mass + (2. * qK0K0)   / mass) / 2.;
  rho(0, 1) = rho(1, 0) = 0;

  const double scale = (s / (4 * _kaonMeanMass * _kaonMeanMass)) - 1;

  b::matrix<cd> M(b::zero_matrix<cd>(2, 2));
  for (uint i = 0; i < _sP.size2(); ++i) {
    const cd fa = 1. / (s - _sP(0, i));
    M += fa * _a[i];
  }
  for (uint i = 0; i < _c.size(); ++i) {
    const cd sc = pow(scale, (int)i);
    M += sc *_c[i];
  }

  // modification: off-diagonal terms set to 0
  M(0, 1) = 0;
  M(1, 0) = 0;

  b::matrix<cd> _T(2, 2);
  rpwa::invertMatrix<cd>(M - cd(0, 1.) * rho, _T);
//  b::matrix<cd> din_inv = gjinverse(din, sing);
//  // std::cout << din_inv << std::endl;
//  if (sing) {std::cerr << "ERROR: SINGULAR" << std::endl;}  // exit(); }

  const cd amp = _T(0, 0);
  return amp;
}

int main () {
  
  TCanvas c1("c1");
  combine(
          SET2(
               draw([](double m)->double{return real(ampAMP(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Real part"),
               SetLineColor(kBlack)),
          SET2(
               draw([](double m)->double{return imag(ampAMP(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Imag part"),
               SetLineColor(kRed) ),
          SET2(
               draw([](double m)->double{return abs(ampAMP(m*m));}, 2*PI_MASS, 2.2),
               SetTitle("Abs part"),
               SetLineColor(kGreen) ) )->Draw("al");
  c1.SaveAs("/tmp/test_AMPK.pdf(");

  SET2(
       draw([](double m)->double{return 180./M_PI*arg(ampAMP(m*m));}, 2*PI_MASS, 2.2),
       SetTitle("Phase shift"),
       SetLineColor(kBlack) )->Draw("al");
  c1.SaveAs("/tmp/test_AMPK.pdf)");

  return 0;
}
