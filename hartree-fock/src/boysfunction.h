#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#ifndef PI
#define PI 3.1415926535897932384626433
#endif

#include <armadillo>

class BoysFunction
{
public:
    BoysFunction();
    void set(double x, int nMax);
    double returnValue(int n);
    double nMax();
    arma::vec F;
private:
    double tabulated(int n, double x);
    double asymptotic(int n, double x);
    double factorial2(int n);

    arma::mat Ftabulated;

    int m_nMax;
};

#endif // BOYSFUNCTION_H
