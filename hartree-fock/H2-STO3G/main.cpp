#include <iostream>
#include <primitive.h>
#include <armadillo>
#include <contracted.h>
#include <electronsystem.h>
#include <hfsolver.h>
#include <nucleus.h>


using namespace std;

int main()
{

    /*
$basis
*
h   STO-3G
*
    3  s
      3.42525091             0.15432897
      0.62391373             0.53532814
      0.16885540             0.44463454
*
$end
     **/

    double weights[] = {0.15432897,0.53532814,0.44463454};
    double as [] = {3.42525091, 0.62391373, 0.16885540};

    arma::vec nuc1(3), nuc2(3);
    nuc1(0) = 0; nuc1(1) = 0; nuc1(2) = 0;
    nuc2(0) = 1.4; nuc1(1) = 0; nuc1(2) = 0;

    Primitive prim1a(weights[0], 0, 0, 0, as[0], nuc1);
    Primitive prim2a(weights[1], 0, 0, 0, as[1], nuc1);
    Primitive prim3a(weights[2], 0, 0, 0, as[2], nuc1);
    Contracted basisA;
    basisA.addPrimitive(prim1a);
    basisA.addPrimitive(prim2a);
    basisA.addPrimitive(prim3a);

    Primitive prim1b(weights[0], 0, 0, 0, as[0], nuc2);
    Primitive prim2b(weights[1], 0, 0, 0, as[1], nuc2);
    Primitive prim3b(weights[2], 0, 0, 0, as[2], nuc2);
    Contracted basisB;
    basisB.addPrimitive(prim1b);
    basisB.addPrimitive(prim2b);
    basisB.addPrimitive(prim3b);

    Nucleus nucA(nuc1, 1);
    Nucleus nucB(nuc2, 1);

    ElectronSystem system(1);
    system.addContracted(basisA);
    system.addContracted(basisB);
    system.addNucleus(nucA);
    system.addNucleus(nucB);

    HFSolver solver(system);
    solver.solve();

    solver.dumpDensity2D("/scratch/densityH2.bin", 200, 0.7, 0, 2, 2);

    cout << "Hello World!" << endl;
    return 0;
}

