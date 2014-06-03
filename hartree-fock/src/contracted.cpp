#include "contracted.h"

Contracted::Contracted()
{
}

std::vector<Primitive> Contracted::primitives()
{
    return m_primitives;
}

void Contracted::addPrimitive(Primitive primitive)
{
    m_primitives.push_back(primitive);
}

double Contracted::evaluate(double x, double y, double z)
{
    arma::vec nucPos = m_primitives.at(0).nucleusPosition();
    double x_diff = x - nucPos(0);
    double y_diff = y - nucPos(1);
    double z_diff = z - nucPos(2);
    double r_squared = (x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
    double result = 0;
    for (Primitive& p : m_primitives) {
        result += p.weight() * pow(x_diff, p.xExponent()) * pow (y_diff, p.yExponent()) * pow(z_diff, p.zExponent())
                        * exp(-p.exponent() * r_squared);
    }
    return result;
}
