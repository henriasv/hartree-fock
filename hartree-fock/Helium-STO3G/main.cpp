#include <iostream>
#include <electronsystem.h>
#include <armadillo>
#include <primitive.h>
#include <hfsolver.h>

using namespace std;

int main()
{



    // Helium STO 3G basis Coefficients, exponents
    /*
    $basis
    *
    he   STO-3G
    *
        3  s
          6.36242139             0.15432897
          1.15892300             0.53532814
          0.31364979             0.44463454
    *
    $end

    */

    // Building STO-3G orbitals for helium
    // 1s orbital
    arma::vec nucleusPosition(3);
    nucleusPosition(0) = 0.0;
    nucleusPosition(1) = 0.0;
    nucleusPosition(2) = 0.0;

    int i, j, k;
    i=j=k=0; // 1s orbital
    double weight = 0.15432897;
    double a = 6.36242139;
    Primitive prim1(weight, i, j, k, a, nucleusPosition);

    weight = 0.53532814;
    a = 1.15892300;
    Primitive prim2(weight, i, j, k, a, nucleusPosition);

    weight = 0.44463454;
    a = 0.31364979;
    Primitive prim3(weight, i, j, k, a, nucleusPosition);

    Contracted orb1s_Helium;
    orb1s_Helium.addPrimitive(prim1);
    orb1s_Helium.addPrimitive(prim2);
    orb1s_Helium.addPrimitive(prim3);

    int charge = 2;
    Nucleus nucleus(nucleusPosition, charge);


    ElectronSystem electronSystem(1);
    electronSystem.addContracted(orb1s_Helium);
    electronSystem.addNucleus(nucleus);

    HFSolver solver(electronSystem);

    std::cout << solver.overlapMatrix() << std::endl;
    std::cout << solver.uncoupledMatrix() << std::endl;
    std::cout << solver.coupledMatrix() << std::endl;

    solver.solve();

    cout << "Hello World!" << endl;
    return 0;
}

