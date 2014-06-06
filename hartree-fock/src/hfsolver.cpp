#include "hfsolver.h"

HFSolver::HFSolver(ElectronSystem & system) :
    m_electronSystem(&system),
    m_convergenceCriterion(1e-8),
    m_maxIterations(10000),
    m_iterationsUsed(0),
    m_energy(0)
{
    //std::cout << " Nuclear_enery_terms " << m_electronSystem->nuclearEnergyTerms() << std::endl;
    setupOverlapMatrix();
    setupUncoupledMatrix();
    setupCoupledMatrix();
    setupCoefficientMatrix();
    setupDensityMatrix();
    setupFockMatrix();
}

double HFSolver::solve()
{
    for (int i = 0; i<m_maxIterations; i++) {
        arma::vec prevFockEnergies = m_fockEnergies;
        advance();
        if (i>0) {
            double convergenceMeasure = arma::sum(arma::abs(m_fockEnergies-prevFockEnergies));
            m_iterationsUsed = i;
            if (convergenceMeasure<m_convergenceCriterion) {
                break;
            }
        }
    }
    return m_energy;
}

void HFSolver::advance()
{
    int np = m_electronSystem->numParticles();
    int no = m_electronSystem->numBasisFunctions();
    setupFockMatrix();
    arma::mat Ctmp;

    arma::vec s;
    arma::mat U;

    arma::eig_sym(s, U, m_overlapMatrix, "std");
    arma::mat V = U*arma::diagmat(1.0/arma::sqrt(s));

    arma::mat Fprime = V.t() * m_fockMatrix *V;
    arma::mat Cprime;
    arma::eig_sym(m_fockEnergies, Cprime, Fprime, "std");

    m_coefficientMatrix = V*Cprime.submat(0, 0, no-1, np-1);
    normalizeCoefficientMatrix();
    setupDensityMatrix();
    m_energy = calcEnergy();
}

void HFSolver::setupUncoupledMatrix()
{
    int nOrbitals = m_electronSystem->numBasisFunctions();
    m_uncoupledMatrix = arma::zeros(nOrbitals, nOrbitals);
    for (int p = 0; p<nOrbitals; p++) {
        for (int q = 0; q<nOrbitals; q++) {
            m_uncoupledMatrix(p, q) = m_electronSystem->uncoupledIntegral(p,q);
        }
    }
}

void HFSolver::setupOverlapMatrix()
{
    int nOrbitals = m_electronSystem->numBasisFunctions();
    m_overlapMatrix = arma::zeros(nOrbitals, nOrbitals);
    for (int p = 0; p<nOrbitals; p++) {
        for (int q = 0; q<nOrbitals; q++) {
            m_overlapMatrix(p, q) = m_electronSystem->overlapIntegral(p,q);
        }
    }
}

void HFSolver::setupDensityMatrix()
{
    m_densityMatrix = 2*m_coefficientMatrix*m_coefficientMatrix.t();
}

void HFSolver::setupCoupledMatrix()
{
    int nOrbitals = m_electronSystem->numBasisFunctions();
    m_coupledMatrix.set_size(nOrbitals, nOrbitals);
    for (int i = 0; i<nOrbitals; i++) {
        for (int j = 0; j<nOrbitals; j++){
            m_coupledMatrix(i, j) = arma::zeros(nOrbitals, nOrbitals);
        }
    }

    for (int p = 0; p<nOrbitals; p++) {
        for (int q = 0; q<nOrbitals; q++) {
            for (int r = 0; r<nOrbitals; r++) {
                for (int s = 0; s<nOrbitals; s++) {
                    m_coupledMatrix(p, r)(q, s) = m_electronSystem->coupledIntegral(p, q, r, s);
                }
            }
        }
    }
}

void HFSolver::setupCoefficientMatrix()
{
    int n = m_electronSystem->numBasisFunctions();
    m_coefficientMatrix = arma::zeros(n, n);
}

void HFSolver::normalizeCoefficientMatrix()
{
    int np = m_electronSystem->numParticles();
    int no = m_electronSystem->numBasisFunctions();
    for (int k = 0; k<np; k++) {
        double factor = 0;
        for (int p = 0; p<no; p++) {
            for (int q = 0; q<no; q++) {
                factor += m_coefficientMatrix(p,k)*m_overlapMatrix(p,q)*m_coefficientMatrix(q,k);
            }
        }
        m_coefficientMatrix.col(k) = m_coefficientMatrix.col(k)/sqrt(factor);
    }
}

void HFSolver::setupFockMatrix()
{
    int n = m_electronSystem->numBasisFunctions();
    m_fockMatrix = arma::zeros(n, n);
    arma::mat& F = m_fockMatrix;
    arma::mat& P = m_densityMatrix;
    arma::mat& h = m_uncoupledMatrix;
    arma::field<arma::mat>& Q = m_coupledMatrix;

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

double HFSolver::calcEnergy()
{
    int n = m_electronSystem->numBasisFunctions();
    double energy = 0;
    arma::mat& P = m_densityMatrix;
    arma::mat& h = m_uncoupledMatrix;
    arma::field<arma::mat>& Q = m_coupledMatrix;

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
    energy += m_electronSystem->nuclearEnergyTerms();
    return energy;
}

int HFSolver::iterationsUsed()
{
    return m_iterationsUsed;
}

/**
 * @brief HFSolver::dumpDensity2D Dump density from xy-plane
 * @param filename
 * @param resolution
 */
void HFSolver::dumpDensity2D(char filename[], int resolution, double x_mid, double y_mid, double x_max, double y_max)
{
    arma::mat densities(resolution, resolution);
    densities.zeros();
    double dx = 2*x_max/resolution;
    double dy = 2*y_max/resolution;
    for (int i = -resolution/2; i<resolution/2; i++) {
        for (int j = -resolution/2; j<resolution/2; j++) {
            densities(i+resolution/2, j+resolution/2) = m_electronSystem->particleDensity(m_coefficientMatrix, i*dx+x_mid, j*dy+y_mid, 0);
        }
    }
    densities.save(filename, arma::raw_binary);
}


arma::mat HFSolver::overlapMatrix()
{
    return m_overlapMatrix;
}

arma::mat HFSolver::uncoupledMatrix()
{
    return m_uncoupledMatrix;
}

arma::field<arma::mat> HFSolver::coupledMatrix()
{
    return m_coupledMatrix;
}



