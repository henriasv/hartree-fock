#include <iostream>
#include <hfsolver.h>
#include <electronsystem.h>
#include <turbomoleparser.h>
#include <nucleus.h>
#include <contracted.h>
#include <cmath>
#include <armadillo>

using namespace std;

int main()
{
    //double bondDistance = 1.809;
    //double bondAngle = 104.5*PI/180;

    int N_theta = 10;
    int N_r = 10;
    double theta_min = 95*PI/180;
    double theta_max = 120*PI/180;
    double r_min = 1.6;
    double r_max = 2.1;

    arma::vec thetas = arma::linspace(theta_min, theta_max, N_theta);
    arma::vec rs = arma::linspace(r_min, r_max, N_r);
    arma::mat energies = arma::zeros(N_r, N_theta);

    int counter = 0;
    int countmax = N_r*N_theta;

    ofstream outfile("/scratch/hfdata/H2O_energies_431g.dat");
    for (int rCnt = 0; rCnt < N_r; rCnt ++) {
        for (int thetaCnt = 0; thetaCnt < N_theta; thetaCnt ++) {

            double bondAngle = thetas(thetaCnt);
            double bondDistance = rs(rCnt);

            arma::vec posH1(3), posH2(3), posO(3);
            posO.zeros();
            posH1.zeros();
            posH1(0) = bondDistance*cos((PI-bondAngle)/2);
            posH1(1) = bondDistance*sin((PI-bondAngle)/2);
            posH2.zeros();
            posH2(0) = -bondDistance*cos((PI-bondAngle)/2);
            posH2(1) = bondDistance*sin((PI-bondAngle)/2);

            Nucleus H1(posH1,1);
            Nucleus H2(posH2,1);
            Nucleus O(posO,8);

            TurboMoleParser parserO("../../data/basis_sets/O_3-21g.txt");
            TurboMoleParser parserH("../../data/basis_sets/H_3-21g.txt");

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
            double energy = solver.solve();
            //solver.dumpDensity2D("/scratch/densityH2O_sym_1045deg431.bin", 200, 0,0, 3, 3);
            energies(rCnt, thetaCnt) = energy;
            std::cout << energy << " " << solver.iterationsUsed() << std::endl;
            counter ++;
            std::cout << counter << "/" << countmax << std::endl;
            outfile << energy << " " << bondDistance << " "<< bondAngle<< std::endl;
        }
    }
}

