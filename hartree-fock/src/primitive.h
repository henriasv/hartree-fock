#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#ifndef PI
#define PI 3.14159265358979323
#endif

#include <armadillo>

class Primitive
{
public:
    explicit Primitive(double weight, int xExponent, int yExponent, int zExponent, double exponent, arma::vec nucleusPosition);
    double exponent() const;
    int zExponent() const;
    int yExponent() const;
    int xExponent() const;
    double weight() const;
    const arma::vec& nucleusPosition() const;
    void normalizePrimitive();

private:
    double m_weight;
    int m_xExponent;
    int m_yExponent;
    int m_zExponent;
    double m_exponent;
    arma::vec m_nucleusPosition;
    int factorial(int n);

};

#endif // PRIMITIVE_H
