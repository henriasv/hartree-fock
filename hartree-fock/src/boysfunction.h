#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#ifndef PI
#define PI 3.1415926535897932384626433
#endif

#include <armadillo>

class BoysFunction
{
public:
    BoysFunction(int angMomMax);
    void setx(double x);
    double returnValue(int n);
private:
    double tabulated(int n, double x);
    double asymptotic(int n, double x);
    double factorial2(int n);

    arma::mat Ftabulated;
    arma::vec F;
    int nMax;
};

#endif // BOYSFUNCTION_H
