#include "pch.h"
#include "Norms.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

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
double Norms::logSquaredDifference(const std::vector<double>& market_prices, const std::vector<double>& model_prices) {
    if (market_prices.size() != model_prices.size()) {
        throw std::invalid_argument("Market and model prices vectors must have the same size.");
    }

    double total_error = 0.0;

    for (std::size_t i = 0; i < market_prices.size(); ++i) {
        double market_price = market_prices[i];
        double model_price = model_prices[i];

        if (market_price <= 0 || model_price <= 0) {
            throw std::invalid_argument("Prices must be positive.");
        }

        double error = std::log(market_price) - std::log(model_price);
        total_error += error * error / market_price;
    }

    return total_error;
}
