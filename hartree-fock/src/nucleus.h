#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <armadillo>

class Nucleus
{
public:
    Nucleus(arma::vec in_nucleusPosition, int in_nuclearCharge);
    arma::vec position();
    int charge();

private:
    arma::vec m_nucleusPosition;
    int m_nuclearCharge;
};

#endif // NUCLEUS_H
