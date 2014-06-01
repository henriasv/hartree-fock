#include "primitive.h"

Primitive::Primitive(double weight, int xExponent, int yExponent, int zExponent, double exponent, arma::vec nucleusPosition):
    m_weight(weight), m_xExponent(xExponent), m_yExponent(yExponent), m_zExponent(zExponent), m_exponent(exponent), m_nucleusPosition(nucleusPosition)
{
    normalizePrimitive();
}

double Primitive::exponent() const
{
    return m_exponent;
}

int Primitive::zExponent() const
{
    return m_zExponent;
}

int Primitive::yExponent() const
{
    return m_yExponent;
}

int Primitive::xExponent() const
{
    return m_xExponent;
}

double Primitive::weight() const
{
    return m_weight;
}

const arma::vec & Primitive::nucleusPosition() const
{
    return m_nucleusPosition;
}

void Primitive::normalizePrimitive()
{
    int i, j, k;
    i = m_xExponent;
    j = m_yExponent;
    k = m_zExponent;
    m_weight = m_weight*pow(2*m_exponent/PI, 0.75)*sqrt(pow(8*m_exponent, i+j+k)*
                                                       factorial(i)*factorial(j)*factorial(k)/
                                                       (factorial(2*i)*factorial(2*j)*factorial(2*k)));
}

int Primitive::factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


