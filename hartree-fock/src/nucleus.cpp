#include "nucleus.h"

Nucleus::Nucleus(arma::vec in_nucleusPosition, int in_nuclearCharge)
{
    m_nucleusPosition = in_nucleusPosition;
    m_nuclearCharge = in_nuclearCharge;
}

arma::vec Nucleus::position()
{
    return m_nucleusPosition;
}

int Nucleus::charge()
{
    return m_nuclearCharge;
}
