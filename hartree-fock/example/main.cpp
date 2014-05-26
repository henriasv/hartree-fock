#include <iostream>
#include <electronicsystem.h>
#include <berylliumhf.h>
#include <hartreefocksolver.h>
#include <unittest++/UnitTest++.h>
#include <primitive.h>
#include <integrator.h>
using namespace std;

int main()
{

    // Choose system
    int n = 2;
    int numElectrons = n;
    int nuclearCharge = n;
    int numBasisFunctions = 2;
    ElectronicSystem* system = new AtomicOrbitals(numElectrons, nuclearCharge, numBasisFunctions);
    // Initialize solver with system
    HartreeFockSolver solver = HartreeFockSolver(system);
    // Solve
    solver.solve(10);
    // Print resultsx

    return 0;
}


