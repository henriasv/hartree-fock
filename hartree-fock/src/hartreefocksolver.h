#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include "electronicsystem.h"
#include <armadillo>

class HartreeFockSolver
{
public:
    HartreeFockSolver(ElectronicSystem* system);
    double solve(int);

protected:
    void advance();
    void setuph();
    void setupS();
    void setupP();
    void setupQ();
    void resetC();
    void setupF();
    double calcEnergy();


    void printMatrices();

    ElectronicSystem* electronicSystem;
    double convergenceCriterion;
    arma::mat densityMatrix;
    arma::mat fockMatrix;
    arma::mat overlapMatrix;
    arma::mat coefficientMatrix;
    arma::mat uncoupledMatrix;
    arma::mat h;
    arma::field<arma::mat> coupledMatrix;
    double m_energy;
};

#endif // HARTREEFOCKSOLVER_H
