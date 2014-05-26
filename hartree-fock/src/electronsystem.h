#ifndef ELECTRONSYSTEM_H
#define ELECTRONSYSTEM_H

#include <armadillo>
#include <contracted.h>
#include <nucleus.h>
#include <vector>
#include <integrator.h>
#include <cmath>

/**
 * @brief The ElectronSystem class
 * "Interface" between gaussian primitives and hartree fock solver. pq-indexes matrixelements over contracted basis functions are requested by the solver. This class delivers them by combining integrals over gaussian primitives.
 */
class ElectronSystem
{
public:
    ElectronSystem(int num_particles);
    double overlapIntegral(int p, int q);
    double coupledIntegral(int p, int q, int r, int s);
    double uncoupledIntegral(int p, int q);
    double nuclearEnergyTerms();
    int nContractedBasisFunctions(); // max value of p, q, r, s
    void addContracted(Contracted in_contracted);
    void addNucleus(Nucleus in_nucleus);
    int numBasisFunctions();
    int numNuclei();
    int numParticles();

private:
    std::vector<Contracted> m_contracted;
    std::vector<Nucleus> m_nuclei;
    Integrator m_integrator;
    int m_numBasisFunctions;
    int m_numNuclei;
    int m_numParticles;
};

#endif // ELECTRONSYSTEM_H
