#ifndef ATOMICORBITALS_H
#define ATOMICORBITALS_H

#include "electronicsystem.h"
#include </usr/include/armadillo>
#include <string>

class AtomicOrbitals : public ElectronicSystem
{
public:
    AtomicOrbitals(int, int, int);
    double overlapIntegral(int p, int q);
    double coupledIntegral(int p, int q, int r, int s);
    double uncoupledIntegral(int p, int q);
    virtual std::string str();

protected:
    arma::field<arma::mat> Q;
    arma::field<arma::mat> createQ(int Z);
    
};

#endif // ATOMICORBITALS_H
