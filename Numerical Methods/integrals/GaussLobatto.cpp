#include "GaussLobatto.h"
#include "../pch.h"
#include <cmath>
#include <iostream>

const double alpha = std::sqrt(2.0 / 3.0);
const double beta = 1.0 / std::sqrt(5.0);
const double x1 = 0.94288241569547971906;
const double x2 = 0.64185334234578130578;
const double x3 = 0.23638319966214988028;

GaussLobattoIntegration::GaussLobattoIntegration()
{
}

GaussLobattoIntegration::~GaussLobattoIntegration()
{

}
double GaussLobattoIntegration::gaussLobattoIntegrate(std::function<double(double)> f, double a, double b, double absAccuracy, double relAccuracy) {
    // Compute absolute tolerance
    double m = (a + b) / 2;
    double h = (b - a) / 2;
    double y1 = f(a);
    double y3 = f(m - alpha * h);
    double y5 = f(m - beta * h);
    double y7 = f(m);
    double y9 = f(m + beta * h);
    double y11 = f(m + alpha * h);
    double y13 = f(b);

    double f1 = f(m - x1 * h);
    double f2 = f(m + x1 * h);
    double f3 = f(m - x2 * h);
    double f4 = f(m + x2 * h);
    double f5 = f(m - x3 * h);
    double f6 = f(m + x3 * h);

    double acc = h * (0.0158271919734801831 * (y1 + y13)
        + 0.0942738402188500455 * (f1 + f2)
        + 0.1550719873365853963 * (y3 + y11)
        + 0.1888215739601824544 * (f3 + f4)
        + 0.1997734052268585268 * (y5 + y9)
        + 0.2249264653333395270 * (f5 + f6)
        + 0.2426110719014077338 * y7);

    // Ensure that the accuracy is met
    if (std::abs(acc) < absAccuracy || std::abs(acc) < relAccuracy * (std::abs(a) + std::abs(b)))
    {
        return acc;
    }
    else
    {
        return this->gaussLobattoIntegrate(f, a, m, absAccuracy, relAccuracy) + this->gaussLobattoIntegrate(f, m, b, absAccuracy, relAccuracy);
    }
}