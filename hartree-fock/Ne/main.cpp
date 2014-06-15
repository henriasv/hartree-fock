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

    arma::vec posNe(3);
    posNe.zeros();

    Nucleus Ne(posNe,10);

    TurboMoleParser parserNe("../../data/basis_sets/Ne_4-21g.txt");

    std::vector<Contracted> Necontracted = parserNe.returnContracted(posNe);

    ElectronSystem system(5); // 10 electrons in argon

    for (Contracted contracted : Necontracted) {
        system.addContracted(contracted);
    }

    system.addNucleus(Ne);

    HFSolver solver(system);
    double energy = solver.solve();
    std::cout << "Energy " << energy << std::endl;
    //solver.dumpDensity2D("/scratch/densityNe_321g.bin", 200, 0, 0, 0.5, 0.5);
}

