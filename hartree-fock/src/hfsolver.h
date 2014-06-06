#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <electronsystem.h>
#include <armadillo>

class HFSolver
{
public:
    HFSolver(ElectronSystem & system);
    double solve();
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
    int iterationsUsed();
    void dumpDensity2D(char filename[], int resolution,double x_mid, double y_mid, double x_max, double y_max);

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
    arma::vec m_fockEnergies;
    int m_iterationsUsed;
    double m_energy;
};

#endif // HFSOLVER_H
