#include "integrator.h"
#include <cmath>

Integrator::Integrator()
{
}

double Integrator::overlapIntegral(Primitive &A, Primitive &B)
{
    Integrator::setupE(A,B);
    double p = A.exponent()+B.exponent();

    int i = A.xExponent();
    int j = B.xExponent();
    int k = A.yExponent();
    int l = B.yExponent();
    int m = A.zExponent();
    int n = B.zExponent();

    double Ex = E(0)(i, j, 0);
    double Ey = E(1)(k, l, 0);
    double Ez = E(2)(m, n, 0);

    double S = Ex*Ey*Ez*pow((PI/p),1.5);

    return S;
}

void Integrator::setupE(const Primitive &A, const Primitive &B)
{
    //E = arma::field<arma::cube>(3);
    // Primitive A brings a, i, k, m, and its position
    // Primitive B brings b, j, l, n, and its position
    // This function vill set up E to be able to return:
    // E[0]^{ij} E[1]^{kl} E[2]^{mn} to calculate an overlap integral S_{AB} = E[0]^{ij} E[1]^{kl} E[2]^{mn}
    // The cubes with more E-values than necessary are just made because this is an efficient way to find the right Es.
    double a = A.exponent();
    double b = B.exponent();
    double mu = a*b/(a+b);
    double p = a+b;

    arma::vec X_AB = A.nucleusPosition()-B.nucleusPosition();
    arma::vec X_P = (a*A.nucleusPosition()+b*B.nucleusPosition())/p;
    arma::vec X_PA = X_P-A.nucleusPosition();
    arma::vec X_PB = X_P-B.nucleusPosition();

    int maxiA[3];
    int maxiB[3];
    int max_t[3];

    maxiA[0] = A.xExponent()+1;
    maxiA[1] = A.yExponent()+1;
    maxiA[2] = A.zExponent()+1;
    maxiB[0] = B.xExponent()+1;
    maxiB[1] = B.yExponent()+1;
    maxiB[2] = B.zExponent()+1;



    for (int dim = 0; dim<3; dim++) {
        max_t[dim] = maxiA[dim]+maxiB[dim]+1;
        E(dim) = arma::cube(maxiA[dim], maxiB[dim], max_t[dim]);
        E(dim)(0, 0, 0) = exp(-mu*X_AB(dim)*X_AB(dim));

        // First build for iA = 0
        for (int iB = 0; iB<maxiB[dim]; iB++) {
            for (int t = 0; t<max_t[dim]; t++) {
                int iA = 0;

                // This is really important, since E_0^00 has been explicitly set, and will be set to zero if recurrence relation is induced.
                if (iA == 0 && iB == 0 && t == 0) {
                    continue;
                }
                double E_t_prev = 0;
                if (checkIndexCombinationForE(iA, iB-1, t-1)) {
                    E_t_prev = E(dim)(iA, iB-1, t-1);
                }
                double E_t_present = 0;
                if (checkIndexCombinationForE(iA, iB-1, t)) {
                    E_t_present = E(dim)(iA, iB-1, t);
                }
                double E_t_next = 0;
                if (checkIndexCombinationForE(iA, iB-1, t+1)) {
                    E_t_next = E(dim)(iA, iB-1, t+1);
                }
                E(dim)(iA, iB, t) = 1/(2*p)*(E_t_prev)+X_PB(dim)*E_t_present + (t+1)*E_t_next;
            }
        }

        // Then build for iA > 0
        for (int iA = 1; iA<maxiA[dim]; iA++){
            for (int iB = 0; iB<maxiB[dim]; iB++){
                for (int t = 0; t<max_t[dim]; t ++) {
                    double E_t_prev = 0;
                    if (checkIndexCombinationForE(iA-1, iB, t-1)) {
                        E_t_prev = E(dim)(iA-1, iB, t-1);
                    }
                    double E_t_present = 0;
                    if(checkIndexCombinationForE(iA-1, iB, t)){
                        E_t_present = E(dim)(iA-1, iB, t);
                    }
                    double E_t_next = 0;
                    if (checkIndexCombinationForE(iA-1, iB, t+1)) {
                        E_t_next = E(dim)(iA-1, iB, t+1);
                    }
                    E(dim)(iA, iB, t) = 1/(2*p)*E_t_prev+X_PA(dim)*E_t_present + (t+1)*E_t_next;
                }
            }
        }
    }
}

bool Integrator::checkIndexCombinationForE(int iA, int iB, int t)
{
    if ( t<0 || t>(iA+iB) || iA<0 || iB<0)
    {
        return false;
    }
    else
    {
        return true;
    }
}
