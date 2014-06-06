#include <iostream>
#include <hfsolver.h>
#include <electronsystem.h>
#include <turbomoleparser.h>
#include <nucleus.h>
#include <contracted.h>
#include <cmath>

using namespace std;

int main()
{

    arma::vec posAr(3);
    posAr.zeros();

    Nucleus Ar(posAr,18);

    TurboMoleParser parserAr("../../data/basis_sets/Ar_3-21g.txt");

    std::vector<Contracted> Arcontracted = parserAr.returnContracted(posAr);

    ElectronSystem system(9); // 10 electrons in argon

    for (Contracted contracted : Arcontracted) {
        system.addContracted(contracted);
    }

    system.addNucleus(Ar);

    HFSolver solver(system);
    solver.solve();
    solver.dumpDensity2D("/scratch/densityAr_321g.bin", 200, 0, 0, 1, 1);
}

