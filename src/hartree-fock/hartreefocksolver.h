#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include "electronicsystem.h"
#include <armadillo>

class HartreeFockSolver
{
public:
    HartreeFockSolver(ElectronicSystem* system);
    void solve(int);

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
};

#endif // HARTREEFOCKSOLVER_H
