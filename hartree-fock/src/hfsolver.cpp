#include "hfsolver.h"

HFSolver::HFSolver(ElectronSystem system) :
    m_electronSystem(&system),
    m_convergenceCriterion(1e-8),
    m_maxIterations(10)
{
    setupOverlapMatrix();
    setupUncoupledMatrix();
    setupCoupledMatrix();
    setupCoefficientMatrix();
    setupDensityMatrix();
    setupFockMatrix();
}

void HFSolver::solve()
{
    std::cout << "In solver" << std::endl;

    // advance
    for (int i = 0; i<m_maxIterations; i++) {
        advance();
    }
}

void HFSolver::advance()
{
    int np = m_electronSystem->numParticles();
    int no = m_electronSystem->numBasisFunctions();

    arma::mat Ctmp;
    arma::vec fockEnergies;

    arma::vec s;
    arma::mat U;

    arma::eig_sym(s, U, m_overlapMatrix);
    arma::mat V = U*arma::diagmat(1.0/arma::sqrt(s));

    arma::mat Fprime = V.t() * m_fockMatrix *V;
    arma::mat Cprime;
    arma::eig_sym(fockEnergies, Cprime, Fprime);

    m_coefficientMatrix = V*Cprime.submat(0, 0, no-1, np-1);
    normalizeCoefficientMatrix();
    setupDensityMatrix();
    double energy = calcEnergy();
    std::cout << "Energy " << energy << std::endl;

    /*
    HartreeFockSolver::advance();
      ElectronSystem* f = electronSystem();
      uint no = f->nBasisFunctions();
      uint nk = f->nParticles() / 2;
      setupFockMatrix();

      vec s;
      mat U;
      eig_sym(s, U, overlapMatrix());

      mat V = U*diagmat(1.0/sqrt(s));

      const mat &F = m_fockMatrix;
      mat Fprime = V.t() * F * V;

      mat Cprime;
      eig_sym(m_fockEnergies, Cprime, Fprime);


      mat &C = m_coefficientMatrix;
      C = V*Cprime.submat(0, 0, no - 1, nk - 1);
      normalizeCoefficientMatrix(nk, C);

      setupDensityMatrix();

      double energy = 0;
      */
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
                factor += m_coefficientMatrix(p, k)*m_overlapMatrix(p,q)*m_coefficientMatrix(q,k);
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
    return energy;
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



