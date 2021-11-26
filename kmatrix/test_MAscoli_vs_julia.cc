// Copyright [2016] Mikhail Mikhasenko
#include <iostream>

#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLegend.h"
#include "TH2D.h"
#include "Math/SpecFuncMathMore.h"

#include "MIsobar.h"
#include "MAscoli.h"
#include "MmatrixK.h"
#include "MProductionPhysics.h"
#include "MParKeeper.h"
#include "mstructures.h"
#include "mintegrate.h"
#include "deflib.h"

int main(int argc, char *argv[]) {
        double mrho = 0.77;  double mrho2 = POW2(mrho );
        double mpi  = 0.14;  double mpi2  = POW2(mpi  );
        double mp   = 0.938; double mp2   = POW2(mp   );
        double Ebeam = 190;
        double t = -0.1;
        double stot = POW2(mpi) + POW2(mp) + 2*mp*Ebeam;
        
        std::cout << "MAscoli::sPionProton: " << MAscoli::sPionProton(0.5, M_PI/2, POW2(mrho),
                                    POW2(5), t, stot, POW2(mpi),
                                    POW2(mp), POW2(mp), POW2(mpi)) << "\n";
        std::cout << "MAscoli::upperPart: " << MAscoli::upperPart(0.5, POW2(mrho), 2, 0,
                                  5, POW2(5), t, POW2(mpi),
                                  POW2(mpi), POW2(mpi)) << "\n";
        std::cout << "MAscoli::getProjectedReducedDeck: " << MAscoli::getProjectedReducedDeck(5, 0, 5, mrho2,
                                                0, 5, POW2(5),
                                                t, mpi2, stot,
                                                mpi2, mp2, mp2,
                                                mpi2) << "\n";

        return 0;
}
