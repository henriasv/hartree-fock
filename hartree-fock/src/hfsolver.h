#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <electronsystem.h>

class HFSolver
{
public:
    HFSolver(ElectronSystem & system);
    void solve();
    void advance();
    void setupUncoupledMatrix();
    void setupOverlapMatrix();
    void setupDensityMatrix();
    void setupCoupledMatrix();
    void setupCoefficientMatrix();
    void normalizeCoefficientMatrix();
    void resetCoefficientMatrix();
    void setupFockMatrix();
    double calcEnergy();

    arma::mat overlapMatrix();
    arma::mat uncoupledMatrix();
    arma::field<arma::mat> coupledMatrix();

    ElectronSystem* m_electronSystem; // To be private

private:

    double m_convergenceCriterion;
    int m_maxIterations;
    arma::mat m_densityMatrix;
    arma::mat m_fockMatrix;
    arma::mat m_overlapMatrix;
    arma::mat m_coefficientMatrix;
    arma::mat m_uncoupledMatrix;
    arma::field<arma::mat> m_coupledMatrix;
};

#endif // HFSOLVER_H
