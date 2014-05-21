#include "integrator.h"
#include <cmath>

Integrator::Integrator()
{
}

double Integrator::kineticIntegral(Primitive &A, Primitive &B)
{
    int i = A.xExponent();
    int j = B.xExponent();
    int k = A.yExponent();
    int l = B.yExponent();
    int m = A.zExponent();
    int n = B.zExponent();

    double a = A.exponent();
    double b = B.exponent();

    double T_x = 4*b*b*overlapIntegral_dim(0, i, j+2, A, B)-2*b*(2*j+1)*overlapIntegral_dim(0, i, j, A, B);
    if (j >=2)
        T_x += j*(j-1)*overlapIntegral_dim(0, i, j-2, A, B);
    double T_y = 4*b*b*overlapIntegral_dim(1, k, l+2, A, B)-2*b*(2*l+1)*overlapIntegral_dim(1, k, l, A, B);
    if (l>=2)
        T_y += l*(l-1)*overlapIntegral_dim(1, k, l-2, A, B);
    double T_z = 4*b*b*overlapIntegral_dim(2, m, n+2, A, B)-2*b*(2*n+1)*overlapIntegral_dim(2, m, n, A, B);
    if (n>=2)
        T_z += n*(n-1)*overlapIntegral_dim(2, m, n-2, A, B);

    double S_x = overlapIntegral_dim(0, i, j, A, B);
    double S_y = overlapIntegral_dim(1, k, l, A, B);
    double S_z = overlapIntegral_dim(2, m, n, A, B);

    double T_ab = -0.5*(T_x*S_y*S_z + S_x*T_y*S_z + S_x*S_y*T_z);
    return T_ab;

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

double Integrator::overlapIntegral_dim(int dim, int iA, int iB, Primitive &A, Primitive &B)
{
    double a = A.exponent();
    double b = B.exponent();
    double p = a+b;
    arma::vec A_pos = A.nucleusPosition();
    arma::vec B_pos = B.nucleusPosition();
    int i, j, k, l, m, n;
    i=j=k=l=m=n=0;
    switch (dim)
    {
    case 0:
        i = iA;
        j = iB;
        break;
    case 1:
        k = iA;
        l = iB;
        break;
    case 2:
        m = iA;
        n = iB;
        break;
    }
    setupE(a, b, A_pos, B_pos, i, j, k, l, m, n);
    return E(dim)(iA, iB, 0)*pow(PI/p, 0.5);
}

void Integrator::setupE(const Primitive &A, const Primitive &B)
{
    //E = arma::field<arma::cube>(3);
    // Primitive A brings a, i, k, m, and its position
    // Primitive B brings b, j, l, n, and its position
    // This function vill set up E to be able to return:
    // E[0]^{ij} E[1]^{kl} E[2]^{mn} to calculate an overlap integral S_{AB} = E[0]^{ij} E[1]^{kl} E[2]^{mn}
    // The cubes with more E-values than necessary are just made because this is an efficient way to find the right Es.

    int i = A.xExponent();
    int k = A.yExponent();
    int m = A.zExponent();
    int j = B.xExponent();
    int l = B.yExponent();
    int n = B.zExponent();
    setupE(A.exponent(), B.exponent(), A.nucleusPosition(), B.nucleusPosition(), i, j, k, l, m, n);
}

void Integrator::setupE(double a, double b, arma::vec A_pos, arma::vec B_pos, int i, int j, int k, int l, int m, int n)
{
    int maxiA[3];
    int maxiB[3];
    int max_t[3];

    maxiA[0] = i+1;
    maxiA[1] = k+1;
    maxiA[2] = m+1;
    maxiB[0] = j+1;
    maxiB[1] = l+1;
    maxiB[2] = n+1;

    double mu = a*b/(a+b);
    double p = a+b;

    arma::vec X_AB = A_pos-B_pos;
    arma::vec X_P = (a*A_pos+b*B_pos)/p;
    arma::vec X_PA = X_P-A_pos;
    arma::vec X_PB = X_P-B_pos;

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
