#include "hartreefocksolver.h"
#include <iostream>


HartreeFockSolver::HartreeFockSolver(ElectronicSystem *system) : electronicSystem(system), convergenceCriterion(1e-8), m_energy(0)
{
    std::cout << "Adding system " << system->str().c_str() << " to solver " << std::endl;
    setupQ();
    setuph();
    setupS();
    resetC();
    setupP();
    setupF();

    printMatrices();

    std::cout << "Printed matrices" << std::endl;
}

double HartreeFockSolver::solve(int maxIterations)
{
    std::cout << "In solver" << std::endl;

    // Set up matrices

    // advance
    for (int i = 0; i<maxIterations; i++) {
        advance();
        //printMatrices();
    }
    return m_energy;
}

void HartreeFockSolver::advance()
{
    int n = electronicSystem->numBasisFunctions;
    int numK = electronicSystem->numK;
    setupF();

    arma::mat Ctmp;
    arma::vec fockEnergies;
    eig_sym(fockEnergies, Ctmp, fockMatrix);
    //Ctmp.print("Ctmp");
    coefficientMatrix = Ctmp.submat(0, 0, n-1, numK-1);
    setupP();
    double energy = calcEnergy();
    m_energy = energy;
    std::cout << "Energy" << std::endl << energy << std::endl;
    /*std::cout << "Fock energies" << std::endl << fockEnergies << std::endl;
    std::cout << "Coefficient Matrix" << std::endl << coefficientMatrix << std::endl;
    densityMatrix.print("Density Matrix)");
    std::cout << "Norm of columns in Coefficient matrix" << std::endl;
    for (int i = 0; i< coefficientMatrix.n_cols; i++) {
       std::cout << arma::norm(coefficientMatrix.col(i),2) << std::endl;
    }
    std::cout << fockMatrix << std::endl;
*/
}

void HartreeFockSolver::setuph()
{
    int n = electronicSystem->numBasisFunctions;
    uncoupledMatrix = arma::zeros(n, n);
    for (int p = 0; p<n; p++) {
        for (int q = 0; q<n; q++) {
            uncoupledMatrix(p, q) = electronicSystem->uncoupledIntegral(p, q);
        }
    }
}

void HartreeFockSolver::setupS()
{
    int n = electronicSystem->numBasisFunctions;
    arma::mat &S = overlapMatrix;
    S = arma::eye(n,n);
}


void HartreeFockSolver::setupF()
{
    int n = electronicSystem->numBasisFunctions;
    fockMatrix = arma::zeros(n, n);
    arma::mat& F = fockMatrix;
    arma::mat& P = densityMatrix;
    arma::mat& h = uncoupledMatrix;
    arma::field<arma::mat>& Q = coupledMatrix;

    for (int p = 0; p<n; p++) {
        for (int q = 0; q<n; q++) {
            F(p, q) = h(p, q);
            for(int r = 0; r<n; r++) {
                for(int s = 0; s<n; s++){
                    F(p, q) += P(s, r) * (Q(p, r)(q, s)-0.5*Q(p, r)(s, q));
                }
            }
        }
    }
}

double HartreeFockSolver::calcEnergy()
{
    int n = electronicSystem->numBasisFunctions;
    double energy = 0;
    arma::mat& P = densityMatrix;
    arma::mat& h = uncoupledMatrix;
    arma::field<arma::mat>& Q = coupledMatrix;

    for (int p = 0; p<n; p++) {
        for (int q = 0; q<n; q++) {
            energy += P(p, q) * h(p, q);
            for (int r = 0; r<n; r++) {
                for (int s = 0; s<n; s++) {
                    energy += 0.5*P(p,q)*P(s,r)*(Q(p,r)(q,s)-0.5*Q(p,r)(s,q));
                }
            }
        }
    }
    return energy;
}

void HartreeFockSolver::setupP()
{
    arma::mat &P = densityMatrix;
    arma::mat &C = coefficientMatrix;
    P = 2 * C * C.t();
}

void HartreeFockSolver::setupQ()
{
    int n = electronicSystem->numBasisFunctions;
    coupledMatrix.set_size(n, n);
    for (int p = 0; p<n; p++)
        for (int q = 0; q<n; q++)
            coupledMatrix(p, q) = arma::zeros(n, n);

    for (int p = 0; p<n; p++)
        for (int q = 0; q<n; q++)
            for (int r = 0; r<n; r++)
                for (int s = 0; s<n; s++)
                    coupledMatrix(p, r)(q, s) = electronicSystem->coupledIntegral(p, q, r, s);

}

void HartreeFockSolver::resetC()
{
    int n = electronicSystem->numBasisFunctions;
    int K = electronicSystem->numK;
    arma::mat &C = coefficientMatrix;
    C = arma::zeros(n, K);
}

void HartreeFockSolver::printMatrices()
{
    std::cout << "Fock matrix" << std::endl << fockMatrix << std::endl
              << "Coefficient Matrix" << std::endl << coefficientMatrix << std::endl
              << "Density Matrix" << std::endl << densityMatrix << std::endl
              << "Overlap Matrix" << std::endl << overlapMatrix << std::endl
              << "Uncoupled Matrix" << std::endl << uncoupledMatrix << std::endl;
}
