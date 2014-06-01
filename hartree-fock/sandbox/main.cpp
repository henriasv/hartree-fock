#include <turbomoleparser.h>
#include <armadillo>
#include <hfsolver.h>
#include <electronsystem.h>
#include <contracted.h>

using namespace std;

int main()
{
    char infile[] = "/Users/henriksveinsson/git_repos/hartree-fock/hartree-fock/data/basis_sets/Be_gto6g.txt";
    TurboMoleParser parser;// = new TurboMoleParser;
    std::cout << infile << std::endl;
    parser.setFile(infile);

    arma::vec nuc1pos(3);
    nuc1pos.zeros();
    std::vector<Contracted> contracted1 = parser.returnContracted(nuc1pos);
    Nucleus nuc1(nuc1pos, 4);

    arma::vec nuc2pos(3);
    nuc2pos.zeros();
    nuc2pos(0) = 4.63;
    std::vector<Contracted> contracted2 = parser.returnContracted(nuc2pos);
    Nucleus nuc2(nuc2pos, 4);

    ElectronSystem system(4);
    for (Contracted contracted : contracted1) {
        system.addContracted(contracted);
    }
    for (Contracted contracted : contracted2) {
        system.addContracted(contracted);
    }
    system.addNucleus(nuc1);
    system.addNucleus(nuc2);

    HFSolver solver(system);
    solver.solve();

    return 0;
}

