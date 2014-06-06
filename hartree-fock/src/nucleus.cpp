#include "nucleus.h"
static int num = 0;
Nucleus::Nucleus(arma::vec in_nucleusPosition, int in_nuclearCharge):
    m_nucleusPosition(in_nucleusPosition),
    m_nuclearCharge(in_nuclearCharge)
{
}

const arma::vec & Nucleus::position() const
{
    //std::cout << "Trying to return " << std::endl;
    //std::cout << m_nucleusPosition  << std::endl;
    return m_nucleusPosition;
}

int Nucleus::charge()
{
    return m_nuclearCharge;
}
