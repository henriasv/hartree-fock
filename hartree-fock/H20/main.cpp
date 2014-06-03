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
    double bondDistance = 1.809;
    double bondAngle = 104.5*PI/180;

    arma::vec posH1(3), posH2(3), posO(3);
    posO.zeros();
    posH1.zeros(); posH1(0) = bondDistance;
    posH2.zeros(); posH2(0) = bondDistance*cos(bondAngle); posH2(1) = bondDistance*sin(bondAngle);
    std::cout << posH2 << std::endl;

    Nucleus H1(posH1,1);
    Nucleus H2(posH2,1);
    Nucleus O(posO,8);

    TurboMoleParser parserO("../../data/basis_sets/O_4-31g.txt");
    TurboMoleParser parserH("../../data/basis_sets/H_4-31g.txt");

    std::vector<Contracted> H1contracted = parserH.returnContracted(posH1);
    std::vector<Contracted> H2contracted = parserH.returnContracted(posH2);
    std::vector<Contracted> Ocontracted = parserO.returnContracted(posO);

    ElectronSystem system(5); // 2 H-electrons, 8 O-electrons. Divided by 2 because of spin symmetry

    for (Contracted contracted : H1contracted) {
        system.addContracted(contracted);
    }
    for (Contracted contracted : H2contracted) {
        system.addContracted(contracted);
    }
    for (Contracted contracted : Ocontracted) {
        system.addContracted(contracted);
    }

    system.addNucleus(H1);
    system.addNucleus(H2);
    system.addNucleus(O);

    HFSolver solver(system);
    solver.solve();
    solver.dumpDensity2D("/scratch/densityH2O.bin", 800, 0,0, 3, 3);
}

