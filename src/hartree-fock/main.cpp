#include <iostream>
#include "electronicsystem.h"
#include "berylliumhf.h"
#include "hartreefocksolver.h"

using namespace std;

int main()
{
    // Choose system
    int n = 2;
    int numElectrons = n;
    int nuclearCharge = n;
    int numBasisFunctions = 4;
    ElectronicSystem* system = new AtomicOrbitals(numElectrons, nuclearCharge, numBasisFunctions);
    // Initialize solver with system
    HartreeFockSolver solver = HartreeFockSolver(system);
    // Solve
    solver.solve(10);
    // Print results

    return 0;
}

