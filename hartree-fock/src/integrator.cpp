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

double Integrator::electronElectronIntegral(Primitive &A, Primitive &B, Primitive &C, Primitive &D)
{


    int i = A.xExponent();
    int j = B.xExponent();
    int k = A.yExponent();
    int l = B.yExponent();
    int m = A.zExponent();
    int n = B.zExponent();

    int i2 = C.xExponent();
    int j2 = D.xExponent();
    int k2 = C.yExponent();
    int l2 = D.yExponent();
    int m2 = C.zExponent();
    int n2 = D.zExponent();

    double a = A.exponent();
    double b = B.exponent();
    double c = C.exponent();
    double d = D.exponent();

    double p = a+b;
    double q = c+d;

    int tMax = i+j; int uMax = k+l; int vMax = m+n;
    int tauMax = i2+j2; int nuMax = k2+l2; int phiMax = m2+n2;

    setupE(A, B);
    setupE2(C, D);
    setupHermiteIntegrals(A, B, C, D);

    double integral = 0;
    for (int t = 0; t<tMax+1; t++) {
        for (int u = 0; u<uMax+1; u++) {
            for (int v = 0; v<vMax+1; v++) {
                for (int tau = 0; tau<tauMax+1; tau++) {
                    for (int nu = 0; nu<nuMax+1; nu++) {
                        for (int phi = 0; phi<phiMax+1; phi++){
                            integral += E(0)(i,j,t)*E2(0)(i2,j2,tau) * E(1)(k,l,u)*E2(1)(k2,l2,nu) * E(2)(m,n,v)*E2(2)(m2,n2,phi) * pow(-1, tau+nu+phi) * R(0)(t+tau, u+nu, v+phi);
                        }
                    }
                }
            }
        }
    }
    integral *= 2*pow(PI, 5.0/2)/(p*q*sqrt(p+q));
    return integral;
}

double Integrator::nuclearElectronIntegral(Primitive &A, Primitive &B, const arma::vec nuc_pos)
{
    int i = A.xExponent();
    int j = B.xExponent();
    int k = A.yExponent();
    int l = B.yExponent();
    int m = A.zExponent();
    int n = B.zExponent();

    double a = A.exponent();
    double b = B.exponent();
    double p = a+b;

    int tMax = i+j;
    int uMax = k+l;
    int vMax = m+n;

    setupHermiteIntegrals(a, b, A.nucleusPosition(), B.nucleusPosition(), nuc_pos, i, j, k, l, m, n);
    setupE(A, B);

    double integral = 0;
    for (int t = 0; t<tMax+1; t++) {
        for (int u = 0; u<uMax+1; u++) {
            for (int v = 0; v<vMax+1; v++) {
                integral += E(0)(i, j, t)*E(1)(k, l, u)*E(2)(m, n, v)*R(0)(t, u, v);
            }
        }
    }
    integral *= 2*PI/p;
    return integral;
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

void Integrator::setupE(double a, double b, const arma::vec& A_pos, const arma::vec& B_pos, int i, int j, int k, int l, int m, int n)
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

void Integrator::setupE2(const Primitive &C, const Primitive &D)
{
    int i = C.xExponent();
    int k = C.yExponent();
    int m = C.zExponent();
    int j = D.xExponent();
    int l = D.yExponent();
    int n = D.zExponent();
    setupE2(C.exponent(), D.exponent(), C.nucleusPosition(), D.nucleusPosition(), i, j, k, l, m, n);
}

void Integrator::setupE2(double a, double b, const arma::vec& C_pos, const arma::vec& D_pos, int iC, int iD, int jC, int jD, int kC, int kD)
{
    // This function is a slight hack, E is temporarly stored, E2 is created in E, and then put to E2. Then E is restored. This is just because I didnt think there would be necessary with two E's when i first implemented E.
    arma::field<arma::cube> E_tmp=E;
    setupE(a, b, C_pos, D_pos, iC, iD, jC, jD, kC, kD);
    E2 = E;
    E = E_tmp;
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

void Integrator::setupHermiteIntegrals(double a, double b, const arma::vec &A_pos, const arma::vec &B_pos, const arma::vec& nuc_pos, int iA, int iB, int jA, int jB, int kA, int kB)
{
    int tMax = iA+iB;
    int uMax = jA+jB;
    int vMax = kA+kB;
    int nMax = tMax+uMax+vMax+1;

    int subMax = std::max(tMax, std::max(uMax, vMax));

    double p = a+b;

    arma::vec P = (a*A_pos + b*B_pos)/p;
    arma::vec PC = P-nuc_pos;

    boys.set(p*arma::dot(PC, PC), nMax);

    R.set_size(nMax+1);
    for (int i = 0; i<nMax+1; i++) {
        R(i) = arma::zeros(subMax+1, subMax+1, subMax+1);
    }

    // Setup R_000N
    for (int n = 0; n<nMax+1; n++)
    {
        R(n)(0, 0, 0) = pow(-2*p, n)*boys.returnValue(n);
    }

    for (int tuvSum = 1; tuvSum<nMax; tuvSum ++) {
        for (int n = 0; n<nMax-tuvSum; n++) {
            for (int t = 0; t<subMax+1; t++) {
                for (int u = 0; u<subMax+1; u++) {
                    for (int v = 0; v<subMax+1; v++) {
                        // Check if this combination is available now:
                        if (t+u+v != tuvSum || t+u+v == 0) {
                            continue;
                        }
                        // The largest element of t, u, v is the one that can absorb subtraction! There are several ways to the elements of R, but the one here should work.
                        int largestElement = std::max(t,std::max(u, v));
                        if (largestElement == t) {
                            R(n)(t, u, v) = PC(0)*R(n+1)(t-1, u, v);
                            if (t>=2)
                                R(n)(t,u,v) += (t-1)*R(n+1)(t-2, u,v);
                        }
                        else if (largestElement == u) {
                            R(n)(t, u, v) = PC(1)*R(n+1)(t, u-1, v);
                            if(u>=2)
                                R(n)(t,u,v) += (u-1)*R(n+1)(t, u-2, v);
                        }
                        else {
                            R(n)(t, u, v) = PC(2)*R(n+1)(t, u, v-1);
                            if (v>=2)
                                R(n)(t,u,v) += (v-1)*R(n+1)(t, u, v-2);

                        }
                    }
                }
            }
        }
    }
}

void Integrator::setupHermiteIntegrals(const Primitive &A, const Primitive &B, const Primitive &C, const Primitive &D)
{
    double a = A.exponent();
    double b = B.exponent();
    double c = C.exponent();
    double d = D.exponent();

    double p = a+b;
    double q = c+d;
    double alpha = p*q/(p+q);

    arma::vec P = (a*A.nucleusPosition()+b*B.nucleusPosition())/p;
    arma::vec Q = (c*C.nucleusPosition()+d*D.nucleusPosition())/q;
    arma::vec PC = P-Q; // This is really important P-Q, not Q-P!

    int tMax = A.xExponent()+B.xExponent()+C.xExponent()+D.xExponent();
    int uMax = A.yExponent()+B.yExponent()+C.yExponent()+D.yExponent();
    int vMax = A.zExponent()+B.zExponent()+C.zExponent()+D.zExponent();
    int nMax = tMax+uMax+vMax+1;
    int subMax = std::max(tMax, std::max(uMax, vMax));

    boys.set(alpha*arma::dot(PC, PC), nMax);

    R.set_size(nMax+1);
    for (int i = 0; i<nMax+1; i++) {
        R(i) = arma::zeros(subMax+1, subMax+1, subMax+1);
    }

    // Setup R_000N
    for (int n = 0; n<nMax+1; n++)
    {
        R(n)(0, 0, 0) = pow(-2*(alpha), n)*boys.returnValue(n);
    }

    for (int tuvSum = 1; tuvSum<nMax; tuvSum ++) {
        for (int n = 0; n<nMax-tuvSum; n++) {
            for (int t = 0; t<subMax+1; t++) {
                for (int u = 0; u<subMax+1; u++) {
                    for (int v = 0; v<subMax+1; v++) {
                        // Check if this combination is available now:
                        if (t+u+v != tuvSum || t+u+v == 0) {
                            continue;
                        }
                        // The largest element of t, u, v is the one that can absorb subtraction! There are several ways to the elements of R, but the one here should work.
                        int largestElement = std::max(t,std::max(u, v));
                        if (largestElement == t) {
                            R(n)(t, u, v) = PC(0)*R(n+1)(t-1, u, v);
                            if (t>=2)
                                R(n)(t,u,v) += (t-1)*R(n+1)(t-2, u,v);
                        }
                        else if (largestElement == u) {
                            R(n)(t, u, v) = PC(1)*R(n+1)(t, u-1, v);
                            if(u>=2)
                                R(n)(t,u,v) += (u-1)*R(n+1)(t, u-2, v);
                        }
                        else {
                            R(n)(t, u, v) = PC(2)*R(n+1)(t, u, v-1);
                            if (v>=2)
                                R(n)(t,u,v) += (v-1)*R(n+1)(t, u, v-2);

                        }
                    }
                }
            }
        }
    }
}
