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
    double electronElectronIntegral(Primitive &A, Primitive &B, Primitive &C, Primitive &D);
    double nuclearElectronIntegral(Primitive &A, Primitive &B, const arma::vec nuc_pos);
    void setupHermiteIntegrals(double a, double b, const arma::vec& A_pos, const arma::vec& B_pos, const arma::vec& nuc_pos, int iA, int iB, int jA, int jB, int kA, int kB); //Should be private
    void setupHermiteIntegrals(const Primitive& A, const Primitive& B, const Primitive& C, const Primitive& D);


private:
    void setupE(const Primitive& A, const Primitive& B);
    void setupE(double a, double b, const arma::vec& A_pos, const arma::vec& B_pos, int iA, int iB, int jA, int jB, int kA, int kB);
    void setupE2(const Primitive& C, const Primitive& D);
    void setupE2(double a, double b, const arma::vec& C_pos, const arma::vec& D_pos, int iC, int iD, int jC, int jD, int kC, int kD);
    bool checkIndexCombinationForE(int iA, int iB, int t);



    arma::field<arma::cube> E = arma::field<arma::cube>(3);
    arma::field<arma::cube> E2 = arma::field<arma::cube>(3);
    arma::field<arma::cube> R;
    BoysFunction boys;

};

#endif // INTEGRATOR_H
