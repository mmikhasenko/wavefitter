// Copyright [14.07.2015] Misha Mikhasenko

#ifndef _MISOBARPIPIS_H_
#define _MISOBARPIPIS_H_

#include "MIsobar.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/triangular.hpp>

//////////////////////////////////////////////////////////////////////////////
/// Brief Au-Morgan-Pennington parameterization of pi pi s-wave
/// [K.L. Au et al, Phys. Rev. D35, 1633] M solution.
/// 
/// Kachaev has splitted pipiS and f0 waves
/// 
/// a small modification by setting the off-diagonal elements of the M-matrix
/// to zero was introduced by Munuch people (Florinan, Boris?)


class MIsobarPiPiS : public MIsobar {
 public:
  explicit MIsobarPiPiS();

 private:
  std::vector<double> _a;
  std::vector<double> _c;
  std::vector<double> _sP;
  
 public:
  virtual double U(double s12) const;
  virtual cd     U(cd s)       const;
  virtual cd     T(double s)   const;
  virtual cd     T(cd s)       const;
  virtual cd ToneVertex(double s) const;

};

#endif  // SRC_MISOBARPIPIS_H_
