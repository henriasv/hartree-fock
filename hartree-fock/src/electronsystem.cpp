#include "electronsystem.h"


ElectronSystem::ElectronSystem(int numParticles) :
    m_numParticles(numParticles),
    m_numBasisFunctions(0),
    m_numNuclei(0)
{
}

double ElectronSystem::overlapIntegral(int p, int q)
{
    double result = 0;
    Contracted& pBF = m_contracted.at(p);
    Contracted& qBF = m_contracted.at(q);
    for (Primitive pP : pBF.primitives()) {
        for (Primitive qP : qBF.primitives()) {
            result += pP.weight()*qP.weight()*m_integrator.overlapIntegral(pP, qP);
        }
    }
    return result;
}

double ElectronSystem::coupledIntegral(int p, int q, int r, int s)
{
    double result = 0;
    Contracted& pBF = m_contracted.at(p);
    Contracted& qBF = m_contracted.at(q);
    Contracted& rBF = m_contracted.at(r);
    Contracted& sBF = m_contracted.at(s);

    for (Primitive pP : pBF.primitives()) {
        for (Primitive qP : qBF.primitives()) {
            for (Primitive rP : rBF.primitives()) {
                for (Primitive sP : sBF.primitives()) {
                    result += pP.weight()*qP.weight()*rP.weight()*sP.weight()*m_integrator.electronElectronIntegral(pP, qP, rP, sP);
                }
            }
        }
    }
    return result;
}

double ElectronSystem::uncoupledIntegral(int p, int q)
{
    double result = 0;
    Contracted& pBF = m_contracted.at(p);
    Contracted& qBF = m_contracted.at(q);
    for (Primitive pP : pBF.primitives()) {
        for (Primitive qP : qBF.primitives()) {
            result += pP.weight()*qP.weight()*m_integrator.kineticIntegral(pP, qP);
            for (Nucleus nuc : m_nuclei) {
                result -= nuc.charge()*pP.weight()*qP.weight() * m_integrator.nuclearElectronIntegral(pP, qP, nuc.position());
            }
        }
    }
    return result;
}

double ElectronSystem::nuclearEnergyTerms()
{
    double result = 0;
    for(int i = 0; i < m_numNuclei; i++) {
        Nucleus & nuc1 = m_nuclei.at(i);
        for(int j = i + 1; j < m_numNuclei; j++) {
            Nucleus & nuc2 = m_nuclei.at(j);
            result += nuc1.charge()*nuc2.charge() / sqrt(arma::dot(nuc1.position() - nuc2.position() , nuc1.position() - nuc2.position()));
        }
    }
    return result;
}




void ElectronSystem::addContracted(Contracted in_contracted)
{
    m_contracted.push_back(in_contracted);
    m_numBasisFunctions ++;
}

void ElectronSystem::addNucleus(Nucleus in_nucleus)
{
    m_nuclei.push_back(in_nucleus);
    m_numNuclei++;
}

int ElectronSystem::numBasisFunctions()
{
    return m_numBasisFunctions;
}

int ElectronSystem::numNuclei()
{
    return m_numNuclei;
}

int ElectronSystem::numParticles()
{
    return m_numParticles;
}
