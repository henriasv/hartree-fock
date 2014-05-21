#include <primitive.h>
#include <integrator.h>
#include <unittest++/UnitTest++.h>
#include <boysfunction.h>
#include <armadillo>

int main()
{
    return UnitTest::RunAllTests();
}

TEST(hermite_overlap_integral_simple)
{
    arma::vec A(3); arma::vec B(3);
    Integrator integrator;
    int i, j, k, l, m, n;
    double a, b, weight;
    a = 0.2; weight=1.0; i=0; k=0; m=0; A(0) = 1.2; A(1) = 2.3; A(2) = 3.4;
    b = 0.3; weight=1.0; j=0; l=0; n=0; B(0) = -1.3; B(1) = 1.4; B(2) = -2.4;
    Primitive primitiveA(weight, i, j, m, a, A);
    Primitive primitiveB(weight, j, l, n, b, B);
    CHECK_CLOSE(1.191723635809e-1, integrator.overlapIntegral(primitiveA, primitiveB), 1e-5);
}

TEST(hermite_overlap_integral_general)
{
    // PrimitiveA:
    arma::vec A(3); arma::vec B(3);
    Integrator integrator;
    int i, j, k, l, m, n;
    double a, b, weight;
    a = 0.2;
    weight = 1;
    i = m = 0;
    k = 1;
    A(0) = 1.2; A(1) = 2.3; A(2) = 3.4;
    // PrimitiveB:
    b = 0.3;
    weight = 1;
    j = 0;
    l = n = 1;
    B(0) = -1.3; B(1) = 1.4; B(2) = -2.4;
    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    CHECK_CLOSE(2.2273219e-1, integrator.overlapIntegral(primitiveA,
    primitiveB), 1e-5);
}

TEST(hermite_overlap_integral_general2)
{
    Integrator integrator;
    arma::vec A(3); arma::vec B(3);
    int i, j, k, l, m, n;
    double a, b, weight;
    // PrimitiveA:
    a = 0.2;
    weight = 1;
    i = m = 0;
    k = 2;
    A(0) = 1.2; A(1) = 2.3; A(2) = 3.4;
    // PrimitiveB:
    b = 0.3;
    weight = 1;
    j = l = 1;
    n = 0;
    B(0) = -1.3; B(1) = 1.4; B(2) = -2.4;
    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    CHECK_CLOSE(-7.329386373895e-02, integrator.overlapIntegral(primitiveA, primitiveB), 1e-5);
}

TEST(kinetic_integral1)
{
    Integrator integrator;
    arma::vec A(3); arma::vec B(3);
    int i, j, k, l, m, n;
    double a, b, weight;
    // PrimitiveA
    a = 0.2;
    weight = 1;
    i = k = m = 0;
    A(0) = 1.2; A(1) = 2.3; A(2) = 3.4;
    // PrimitiveB
    b = 0.3;
    weight = 1;
    j = l = n = 0;
    B(0) = -1.3; B(1) = 1.4; B(2) = -2.4;
    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    CHECK_CLOSE(-9.678702680582e-02, integrator.kineticIntegral(primitiveA,
    primitiveB), 1e-5);
}

TEST(kinetic_integral2)
{
    Integrator integrator;
    arma::vec A(3); arma::vec B(3);
    int i, j, k, l, m, n;
    double a, b, weight;

    a = 0.2;
    weight = 1;
    i = m = 0;
    k = 1;
    A(0) = 1.2; A(1) = 2.3; A(2) = 3.4;

    b = 0.3;
    weight = 1;
    j = 0;
    l = n = 1;
    B(0) = -1.3; B(1) = 1.4; B(2) = -2.4;
    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    CHECK_CLOSE(-8.688217105502e-02, integrator.kineticIntegral(primitiveA, primitiveB), 1e-5);
}

TEST(kinetic_integral3)
{
    Integrator integrator;
    arma::vec A(3); arma::vec B(3);
    int i, j, k, l, m, n;
    double a, b, weight;
    // PrimitiveA:
    a = 0.2;
    weight = 1;
    i = m = 0;
    k = 2;
    A(0) = 1.2; A(1) = 2.3; A(2) = 3.4;
    // PrimitiveB:
    b = 0.3;
    weight = 1;
    j = l = 1;
    n = 0;
    B(0) = -1.3; B(1) = 1.4; B(2) = -2.4;
    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    CHECK_CLOSE(-1.598401092187e-02, integrator.kineticIntegral(primitiveA, primitiveB), 1e-5);
}

TEST(boysfunction_construct)
{
    //BoysFunction boys;
    //boys.set(1.2, 4);
    Integrator integrator;
    arma::vec A(3); arma::vec B(3); arma::vec C(3);
    int i, j, k, l, m, n;
    double a, b, weight;
    a = 0.2;
    weight = 1;
    i = m = 0;
    k = 2;
    A(0) = 1.2; A(1) = 2.3; A(2) = 3.4;
    // PrimitiveB:
    b = 0.3;
    weight = 1;
    j = l = 1;
    n = 0;
    B(0) = -1.3; B(1) = 1.4; B(2) = -2.4;
    C(0) = 0.2; C(1) = 1.9; C(2) = 1.1;
    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    integrator.setupHermiteIntegrals(a, b, A, B, C, i, j, k, l, m, n);
}

TEST(GTOnuclear_attraction_integral)
{
    Integrator integrator;
    int i, j, k, l, m, n;
    i=j=k=l=m=n=0;
    arma::vec A(3); arma::vec B(3); arma::vec C(3);
    A(0) = 1.2; A(1) = 2.3; A(2) = 3.4;
    B(0) = -1.3; B(1) = 1.4; B(2) = -2.4;
    C(0) = 2.3; C(1) = 0.9; C(2) = 3.2;
    double a = 0.2; double b = 0.3;
    double weight = 1;
    Primitive primitiveA(weight, i, k, m, a, A);
    Primitive primitiveB(weight, j, l, n, b, B);
    CHECK_CLOSE(2.788948987251e-02, integrator.nuclearElectronIntegral(primitiveA, primitiveB, C), 1e-5);

    primitiveA = Primitive(1.0, 2, 0, 0, 0.2, A);
    primitiveB = Primitive(1.0, 1, 1, 0, 0.3, B);
    CHECK_CLOSE(4.176920693786e-03, integrator.nuclearElectronIntegral(primitiveA, primitiveB, C), 1e-5);

    primitiveA = Primitive(1,2,0,0, 0.2, A);
    primitiveB = Primitive(1,1,0,1, 0.3, B);
    CHECK_CLOSE(3.878852644576e-02, integrator.nuclearElectronIntegral(primitiveA, primitiveB, C), 1e-5);
}
