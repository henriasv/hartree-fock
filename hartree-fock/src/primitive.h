#ifndef PRIMITIVE_H
#define PRIMITIVE_H

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

private:
    double m_weight;
    int m_xExponent;
    int m_yExponent;
    int m_zExponent;
    double m_exponent;
    arma::vec m_nucleusPosition;

};

#endif // PRIMITIVE_H
