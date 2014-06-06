#include <iostream>
#include <electronsystem.h>
#include <armadillo>
#include <primitive.h>
#include <hfsolver.h>

using namespace std;

int main()
{



    // Beryllium STO 3G basis Coefficients, exponents
    /*
    $basis
    *
    be   STO-3G
    *
        3  s
         30.1678710              0.15432897
          5.4951153              0.53532814
          1.4871927              0.44463454
        3  s
          1.3148331             -0.09996723
          0.3055389              0.39951283
          0.0993707              0.70011547
        3  p
          1.3148331              0.15591627
          0.3055389              0.60768372
          0.0993707              0.39195739
    *
    $end

    */

    // Building STO-3G orbitals for helium
    // 1s orbital
    arma::vec nucleusPosition1(3);
    nucleusPosition1(0) = -4.63/2;
    nucleusPosition1(1) = 0.0;
    nucleusPosition1(2) = 0.0;

    arma::vec nucleusPosition2(3);
    nucleusPosition2(0) = 4.63/2;
    nucleusPosition2(1) = 0.0;
    nucleusPosition2(2) = 0.0;

    std::vector<arma::vec> nuclei;
    nuclei.push_back(nucleusPosition1);
    nuclei.push_back(nucleusPosition2);

    ElectronSystem electronSystem(4);
    for (int nucNum = 0; nucNum <2; nucNum++) {
        arma::vec nucleusPosition(3);
        nucleusPosition = nuclei.at(nucNum);
    // 1s orbital
    int i, j, k;
    i=j=k=0;
    double weights1s[] = {0.15432897, 0.53532814, 0.44463454};
    double a1s[] = {30.1678710, 5.4951153, 1.4871927};

    Primitive prim11s(weights1s[0], i, j, k, a1s[0], nucleusPosition);
    Primitive prim21s(weights1s[1], i, j, k, a1s[1], nucleusPosition);
    Primitive prim31s(weights1s[2], i, j, k, a1s[2], nucleusPosition);

    Contracted orb1s_Beryllium;
    orb1s_Beryllium.addPrimitive(prim11s);
    orb1s_Beryllium.addPrimitive(prim21s);
    orb1s_Beryllium.addPrimitive(prim31s);

    // 2s orbital
    double weights2s[] = {-0.09996723, 0.39951283, 0.70011547};
    double a2s[] = {1.3148331, 0.3055389, 0.0993707};
    Primitive prim12s(weights2s[0], i, j, k, a2s[0], nucleusPosition);
    Primitive prim22s(weights2s[1], i, j, k, a2s[1], nucleusPosition);
    Primitive prim32s(weights2s[2], i, j, k, a2s[2], nucleusPosition);

    Contracted orb2s_Beryllium;
    orb2s_Beryllium.addPrimitive(prim12s);
    orb2s_Beryllium.addPrimitive(prim22s);
    orb2s_Beryllium.addPrimitive(prim32s);

    // 2p orbitals
    double weights2p[] = {0.15591627, 0.60768372, 0.39195739};
    double a2p[] = {1.3148331,0.3055389,0.0993707 };

    Primitive prim12p(weights2p[0], 1, 0, 0, a2p[0], nucleusPosition);
    Primitive prim22p(weights2p[1], 1, 0, 0, a2p[1], nucleusPosition);
    Primitive prim32p(weights2p[2], 1, 0, 0, a2p[2], nucleusPosition);

    Contracted orb2px_Beryllium;
    orb2px_Beryllium.addPrimitive(prim12p);
    orb2px_Beryllium.addPrimitive(prim22p);
    orb2px_Beryllium.addPrimitive(prim32p);

    Primitive prim12py(weights2p[0], 0, 1, 0, a2p[0], nucleusPosition);
    Primitive prim22py(weights2p[1], 0, 1, 0, a2p[1], nucleusPosition);
    Primitive prim32py(weights2p[2], 0, 1, 0, a2p[2], nucleusPosition);

    Contracted orb2py_Beryllium;
    orb2py_Beryllium.addPrimitive(prim12py);
    orb2py_Beryllium.addPrimitive(prim22py);
    orb2py_Beryllium.addPrimitive(prim32py);

    Contracted orb2pz_Beryllium;
    orb2pz_Beryllium.addPrimitive((Primitive(weights2p[0], 0, 0, 1, a2p[0], nucleusPosition)));
    orb2pz_Beryllium.addPrimitive((Primitive(weights2p[1], 0, 0, 1, a2p[1], nucleusPosition)));
    orb2pz_Beryllium.addPrimitive((Primitive(weights2p[2], 0, 0, 1, a2p[2], nucleusPosition)));


    int charge = 4;
    Nucleus nucleus(nucleusPosition, charge);


    electronSystem.addContracted(orb1s_Beryllium);
    electronSystem.addContracted(orb2s_Beryllium);
    electronSystem.addContracted(orb2px_Beryllium);

    electronSystem.addContracted(orb2py_Beryllium);
    electronSystem.addContracted(orb2pz_Beryllium);

    electronSystem.addNucleus(nucleus);

    }


    HFSolver solver(electronSystem);

    solver.solve();
    solver.dumpDensity2D("/scratch/densityBe2_sto3g.bin", 200,0, 0, 5, 5);

    cout << "Hello World!" << endl;
    return 0;
}

