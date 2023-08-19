/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2005 Klaus Spanderen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file gaussianquadratures.hpp
    \brief Integral of a 1-dimensional function using the Gauss quadratures
*/
#include "../pch.h"

#include "gaussianorthogonalpolynomial.hpp"
#include "gammadistribution.hpp"
#include "comparison.hpp"
#include "errors.hpp"
#include <cmath>

namespace QuantLib {

    double GaussianOrthogonalPolynomial::value(size_t n, double x) const {
        if (n > 1) {
            return  (x-alpha(n-1)) * value(n-1, x)
                       - beta(n-1) * value(n-2, x);
        }
        else if (n == 1) {
            return x-alpha(0);
        }

        return 1;
    }

    double GaussianOrthogonalPolynomial::weightedValue(size_t n, double x) const {
        return std::sqrt(w(x))*value(n, x);
    }


    GaussLaguerrePolynomial::GaussLaguerrePolynomial(double s)
    : s_(s) {
        QL_REQUIRE(s > -1.0, "s must be bigger than -1");
    }

    double GaussLaguerrePolynomial::mu_0() const {
        return std::exp(GammaFunction().logValue(s_+1));
    }

    double GaussLaguerrePolynomial::alpha(size_t i) const {
        return 2*i+1+s_;
    }

    double GaussLaguerrePolynomial::beta(size_t i) const {
        return i*(i+s_);
    }

    double GaussLaguerrePolynomial::w(double x) const {
        return std::pow(x, s_)*std::exp(-x);
    }


    GaussHermitePolynomial::GaussHermitePolynomial(double mu)
    : mu_(mu) {
        QL_REQUIRE(mu > -0.5, "mu must be bigger than -0.5");
    }

    double GaussHermitePolynomial::mu_0() const {
        return std::exp(GammaFunction().logValue(mu_+0.5));
    }

    double GaussHermitePolynomial::alpha(size_t) const {
        return 0.0;
    }

    double GaussHermitePolynomial::beta(size_t i) const {
        return (i % 2) != 0U ? double(i / 2.0 + mu_) : double(i / 2.0);
    }

    double GaussHermitePolynomial::w(double x) const {
        return std::pow(std::fabs(x), 2*mu_)*std::exp(-x*x);
    }

    GaussJacobiPolynomial::GaussJacobiPolynomial(double alpha, double beta)
    : alpha_(alpha), beta_ (beta) {
        QL_REQUIRE(alpha_+beta_ > -2.0,"alpha+beta must be bigger than -2");
        QL_REQUIRE(alpha_       > -1.0,"alpha must be bigger than -1");
        QL_REQUIRE(beta_        > -1.0,"beta  must be bigger than -1");
    }

    double GaussJacobiPolynomial::mu_0() const {
        return std::pow(2.0, alpha_+beta_+1)
            * std::exp( GammaFunction().logValue(alpha_+1)
                        +GammaFunction().logValue(beta_ +1)
                        -GammaFunction().logValue(alpha_+beta_+2));
    }

    double GaussJacobiPolynomial::alpha(size_t i) const {
        double num = beta_*beta_ - alpha_*alpha_;
        double denom = (2.0*i+alpha_+beta_)*(2.0*i+alpha_+beta_+2);

        if (close_enough(denom,0.0)) {
            if (!close_enough(num,0.0)) {
                QL_FAIL("can't compute a_k for jacobi integration\n");
            }
            else {
                // l'Hospital
                num  = 2*beta_;
                denom= 2*(2.0*i+alpha_+beta_+1);

                QL_ASSERT(!close_enough(denom,0.0), "can't compute a_k for jacobi integration\n");
            }
        }

        return num / denom;
    }

    double GaussJacobiPolynomial::beta(size_t i) const {
        double num = 4.0*i*(i+alpha_)*(i+beta_)*(i+alpha_+beta_);
        double denom = (2.0*i+alpha_+beta_)*(2.0*i+alpha_+beta_)
                   * ((2.0*i+alpha_+beta_)*(2.0*i+alpha_+beta_)-1);

        if (close_enough(denom,0.0)) {
            if (!close_enough(num,0.0)) {
                QL_FAIL("can't compute b_k for jacobi integration\n");
            } else {
                // l'Hospital
                num  = 4.0*i*(i+beta_)* (2.0*i+2*alpha_+beta_);
                denom= 2.0*(2.0*i+alpha_+beta_);
                denom*=denom-1;
                QL_ASSERT(!close_enough(denom,0.0), "can't compute b_k for jacobi integration\n");
            }
        }
        return num / denom;
    }

    double GaussJacobiPolynomial::w(double x) const {
        return std::pow(1-x, alpha_)*std::pow(1+x, beta_);
    }


    GaussLegendrePolynomial::GaussLegendrePolynomial()
    : GaussJacobiPolynomial(0.0, 0.0) {
    }

    GaussChebyshev2ndPolynomial::GaussChebyshev2ndPolynomial()
    : GaussJacobiPolynomial(0.5, 0.5) {
    }

    GaussChebyshevPolynomial::GaussChebyshevPolynomial()
    : GaussJacobiPolynomial(-0.5, -0.5) {
    }

    GaussGegenbauerPolynomial::GaussGegenbauerPolynomial(double lambda)
    : GaussJacobiPolynomial(lambda-0.5, lambda-0.5){
    }

    double GaussHyperbolicPolynomial::mu_0() const {
        return M_PI;
    }

    double GaussHyperbolicPolynomial::alpha(size_t) const {
        return 0.0;
    }

    double GaussHyperbolicPolynomial::beta(size_t i) const {
        return i != 0U ? M_PI_2 * M_PI_2 * i * i : M_PI;
    }

    double GaussHyperbolicPolynomial::w(double x) const {
        return 1/std::cosh(x);
    }

}

