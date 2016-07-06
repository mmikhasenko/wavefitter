// Copyright [2016] Mikhail Mikhasenko

#include <functional>

#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "MIsobar.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "mstructures.hh"
#include "MRelationHolder.h"



int main(int argc, char *argv[]) {
  DP dat[2];
  dat[0].lrange =  0.0;
  dat[0].rrange =  3.0;
  dat[1].lrange = -3.0;
  dat[1].rrange =  1.0;

  double par = 5.;
  std::function<double(double)> fun[2];
  fun[0] = [&](double s)->double{return exp(-par*s)+s;};
  fun[1] = [&](double s)->double{return exp(-s)+par*s;};
  const uint Np = 20;

  for (uint j = 0; j < 2;  j++) {
    for (uint i = 0; i < Np; i++) {
      double s = dat[j].lrange + (dat[j].rrange-dat[j].lrange)/(Np-1)*i;
      double val = fun[j](s);

      dat[j].data.push_back({s,
            val+(2*gRandom->Rndm()-1),
            2.*gRandom->Rndm()});
    }
  }

  TCanvas c1("c1"); c1.Divide(2, 1);
  c1.cd(1); draw(dat[0])->Draw("ap");
  c1.cd(2); draw(dat[1])->Draw("ap");
  c1.SaveAs("/tmp/a.test_Fit.pdf");

  ////////////////////////////////////////////////////////////////////////
  //////////////////////// f i t ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  MRelationHolder::gI()->AddRelation(dat[0], fun[0]);
  MRelationHolder::gI()->AddRelation(dat[1], fun[1]);

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
      par = pars[0];
      return MRelationHolder::gI()->CalculateChi2();
    }, 1);
  min->SetFunction(functor);
  // specify parameters
  min->SetVariable(0, "par0", -2, 0.01);
  // minimize
  min->Minimize();

  // draw finally
  for (int i=0; i < 2; i++) {
    c1.cd(i+1);
    combine(draw(dat[i]),
            SET2(draw(fun[i], dat[i].lrange, dat[i].rrange),
                 SetLineColor(kRed), SetLineWidth(2) ) )->Draw("a");
  }
  c1.SaveAs("/tmp/b.test_Fit.pdf");

  return 0;
}
