#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#define PI 3.14159265358979323846264338
#include <primitive.h>
#include <armadillo>
#include <boysfunction.h>

class Integrator
{
public:
    Integrator();
    double kineticIntegral(Primitive &A, Primitive &B);
    double overlapIntegral(Primitive &A, Primitive &B);
    double overlapIntegral_dim(int dim, int iA, int iB, Primitive &A, Primitive &B);
    double coupledIntegral(Primitive &A, Primitive &B, Primitive &C, Primitive &D);
    double unCoupledIntegral(Primitive &A, Primitive &B);


private:
    void setupE(const Primitive& A, const Primitive& B);
    void setupE(double a, double b, const arma::vec& A_pos, const arma::vec& B_pos, int iA, int iB, int jA, int jB, int kA, int kB);
    bool checkIndexCombinationForE(int iA, int iB, int t);

    void setupR(double a, double b, const arma::vec& A_pos, const arma::vec& B_pos, const arma::vec nuc_pos, int iA, int iB, int jA, int jB, int kA, int kB);

    arma::field<arma::cube> E = arma::field<arma::cube>(3);
    arma::field<arma::mat> R;
    BoysFunction boys = BoysFunction(10);

};

#endif // INTEGRATOR_H
