/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Klaus Spanderen

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
#include "../pch.h"
#include "discreteintegrals.hpp"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum.hpp>

using namespace boost::accumulators;

namespace QuantLib {

    double DiscreteTrapezoidIntegral::operator()(
        const Array& x, const Array& f)    const {

        const size_t n = f.size();
        QL_REQUIRE(n == x.size(), "inconsistent size_t");

        accumulator_set<double, features<tag::sum> > acc;

        for (size_t i=0; i < n-1; ++i) {
            acc((x[i+1]-x[i])*(f[i]+f[i+1]));
        }

        return 0.5*sum(acc);
    }

    double DiscreteSimpsonIntegral::operator()(
        const Array& x, const Array& f)    const {

        const size_t n = f.size();
        QL_REQUIRE(n == x.size(), "inconsistent size_t");

        accumulator_set<double, features<tag::sum> > acc;

        for (size_t j=0; j < n-2; j+=2) {
            const double dxj   = x[j+1]-x[j];
            const double dxjp1 = x[j+2]-x[j+1];

            const double alpha = -dxjp1*(2*x[j]-3*x[j+1]+x[j+2]);
            const double dd = x[j+2]-x[j];
            const double k = dd/(6*dxjp1*dxj);
            const double beta = dd*dd;
            const double gamma = dxj*(x[j]-3*x[j+1]+2*x[j+2]);

            acc(k*alpha*f[j]+k*beta*f[j+1]+k*gamma*f[j+2]);
        }
        if ((n & 1) == 0U) {
            acc(0.5*(x[n-1]-x[n-2])*(f[n-1]+f[n-2]));
        }

        return sum(acc);
    }


    double DiscreteTrapezoidIntegrator::integrate(
        const ext::function<double (double)>& f, double a, double b) const {
            const Array x(maxEvaluations(), a, (b-a)/(maxEvaluations()-1));
            Array fv(x.size());
            std::transform(x.begin(), x.end(), fv.begin(), f);

            increaseNumberOfEvaluations(maxEvaluations());
            return DiscreteTrapezoidIntegral()(x, fv);
    }

    double DiscreteSimpsonIntegrator::integrate(
        const ext::function<double (double)>& f, double a, double b) const {
            const Array x(maxEvaluations(), a, (b-a)/(maxEvaluations()-1));
            Array fv(x.size());
            std::transform(x.begin(), x.end(), fv.begin(), f);

            increaseNumberOfEvaluations(maxEvaluations());
            return DiscreteSimpsonIntegral()(x, fv);
    }
}
