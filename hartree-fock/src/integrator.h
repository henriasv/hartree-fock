#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#define PI 3.14159265358979323846264338
#include <primitive.h>
#include <armadillo>

class Integrator
{
public:
    Integrator();
    double kineticIntegral();
    double overlapIntegral(Primitive &A, Primitive &B);
    double coupledIntegral(Primitive &A, Primitive &B, Primitive &C, Primitive &D);
    double unCoupledIntegral(Primitive &A, Primitive &B);
    arma::field<arma::cube> E = arma::field<arma::cube>(3);

private:
    void setupE(const Primitive& A, const Primitive& B);
    bool checkIndexCombinationForE(int iA, int iB, int t);


};

#endif // INTEGRATOR_H
