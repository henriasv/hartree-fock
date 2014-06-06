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
    double weights[] = {0.15432897, 0.53532814, 0.44463454};
    double a[] = {6.36242139, 1.15892300, 0.31364979};

    Primitive prim1(weights[0], i, j, k, a[0], nucleusPosition);
    Primitive prim2(weights[1], i, j, k, a[1], nucleusPosition);
    Primitive prim3(weights[2], i, j, k, a[2], nucleusPosition);

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

    solver.solve();
    solver.dumpDensity2D("/scratch/density.bin", 100, 0, 0, 1, 1);

    cout << "Hello World!" << endl;
    return 0;
}

