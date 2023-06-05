#include "pch.h"
#include "Norms.h"
#include <cmath>
#include <algorithm>
double Norms::euclideanNorm(const std::vector<double>& v) {
    double norm = 0.0;
    for (double x : v) {
        norm += x * x;
    }
    return sqrt(norm);
}
double Norms::manhattanNorm(const std::vector<double>& v) {
    double norm = 0.0;
    for (double x : v) {
        norm += std::abs(x);
    }
    return norm;
}
double Norms::infinityNorm(const std::vector<double>& v) {
    double norm = 0.0;
    for (double x : v) {
        norm = (norm > fabs(x)) ? (norm) : (fabs(x));
    }
    return norm;
}   
