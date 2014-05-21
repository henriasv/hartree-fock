#include "boysfunction.h"

BoysFunction::BoysFunction()
{
    m_nMax = 0;
    Ftabulated.load("../../data/boys_tabulated.dat"); // This file must be in ../data/
}

void BoysFunction::set(double x, int nMax)
{
    m_nMax = nMax;
    F = arma::zeros(m_nMax+1);

    if (x<=50)
    {
        F(m_nMax) = tabulated(m_nMax, x);
    }
    else
    {
        F(m_nMax) = asymptotic(m_nMax, x);
    }

    double ex = exp(-x);

    for (int n = m_nMax; n>0; n--)
    {
        F(n-1) = (2*x*F(n)+ex)/(2*n-1);
    }
}

double BoysFunction::returnValue(int n)
{
    return F(n);
}

double BoysFunction::nMax()
{
    return m_nMax;
}

double BoysFunction::tabulated(int n, double x)
{
    int nxVals = Ftabulated.n_rows;
    double dx = 50.0/(nxVals -1);
    int xIndex = int ((x+0.5*dx)/dx);
    double xt = xIndex*dx;
    double Dx = x-xt;

    double value = 0;
    double factorial = 1;

    for (int k = 0; k<7; k++)
    {
        if (k != 0) {
            factorial *=k;
        }
        value += Ftabulated(xIndex, n+k)*pow(-Dx, k)/factorial;
    }
    return value;
}

double BoysFunction::asymptotic(int n, double x)
{
    return factorial2(n)*sqrt(PI/(pow(x, 2*n+1)))/(pow(2, n+1));
}

double BoysFunction::factorial2(int n)
{
    double value = 1;
    double i = 1;

    while (i<2*n-1)
    {
        i += 2;
        value *= i;
    }
    return value;
}
