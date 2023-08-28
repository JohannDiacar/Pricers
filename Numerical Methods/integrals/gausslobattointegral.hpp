/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Klaus Spanderen

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

/*! \file gausslobattointegral.hpp
    \brief integral of a one-dimensional function using the adaptive
    Gauss-Lobatto integral
*/

#ifndef quantlib_gauss_lobatto_integral_hpp
#define quantlib_gauss_lobatto_integral_hpp

#include "errors.hpp"
#include "null.hpp"
#include "integral.hpp"

namespace QuantLib {

    //! Integral of a one-dimensional function
    /*! Given a target accuracy \f$ \epsilon \f$, the integral of
        a function \f$ f \f$ between \f$ a \f$ and \f$ b \f$ is
        calculated by means of the Gauss-Lobatto formula
    */

    /*! References:
       This algorithm is a C++ implementation of the algorithm outlined in

       W. Gander and W. Gautschi, Adaptive Quadrature - Revisited.
       BIT, 40(1):84-101, March 2000. CS technical report:
       ftp.inf.ethz.ch/pub/publications/tech-reports/3xx/306.ps.gz

       The original MATLAB version can be downloaded here
       http://www.inf.ethz.ch/personal/gander/adaptlob.m
    */

    class GaussLobattoIntegral : public Integrator {
    public:
        GaussLobattoIntegral(size_t maxIterations,
                             double absAccuracy,
                             double relAccuracy = Null<double>(),
                             bool useConvergenceEstimate = true);

      
        double integrate(const ext::function<double(double)>& f, double a, double b) const override;
    protected:
        double adaptivGaussLobattoStep(const ext::function<double (double)>& f,
                                     double a, double b, double fa, double fb,
                                     double is) const;
        double calculateAbsTolerance(const ext::function<double (double)>& f,
                                   double a, double b) const;

        double relAccuracy_;
        const bool useConvergenceEstimate_;
        const static double alpha_, beta_, x1_, x2_, x3_;
    };
}
#endif
