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
    double bondDistance = 3.062;
    double bondAngle = 144*PI/180;

    arma::vec posO1(3), posO2(3), posSi(3);
    posSi.zeros();
    posO1.zeros(); posO1(0) = bondDistance;
    posO2.zeros(); posO2(0) = bondDistance*cos(bondAngle); posO2(1) = bondDistance*sin(bondAngle);
    std::cout << posO2 << std::endl;

    Nucleus O1(posO1,1);
    Nucleus O2(posO2,1);
    Nucleus Si(posSi,8);

    TurboMoleParser parserSi("../../data/basis_sets/Si_3-21g.txt");
    TurboMoleParser parserO("../../data/basis_sets/O_3-21g.txt");

    std::vector<Contracted> O1contracted = parserO.returnContracted(posO1);
    std::vector<Contracted> O2contracted = parserO.returnContracted(posO2);
    std::vector<Contracted> Sicontracted = parserSi.returnContracted(posSi);

    ElectronSystem system(15); // 14 Si-electrons, 2*8 = 16 O-electrons

    for (Contracted contracted : O1contracted) {
        system.addContracted(contracted);
    }
    for (Contracted contracted : O2contracted) {
        system.addContracted(contracted);
    }
    for (Contracted contracted : Sicontracted) {
        system.addContracted(contracted);
    }

    system.addNucleus(O1);
    system.addNucleus(O2);
    system.addNucleus(Si);

    HFSolver solver(system);
    solver.solve();
    solver.dumpDensity2D("/scratch/densitySiO2.bin", 400, 0,0, 5, 5);
}

