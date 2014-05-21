#ifndef ELECTRONICSYSTEM_H
#define ELECTRONICSYSTEM_H

#include <string>

class ElectronicSystem
        /** Effectively defining a hamiltonian and some actions on that hamiltonian.
         *
         **/
{
public:
    ElectronicSystem();
    virtual double overlapIntegral(int p, int q) = 0;
    virtual double coupledIntegral(int p, int q, int r, int s) = 0;
    virtual double uncoupledIntegral(int p, int q) = 0;
    virtual std::string str() = 0;
    int numBasisFunctions;
    int numParticles;
    int numK;
    int Z;
};

#endif // ELECTRONICSYSTEM_H
