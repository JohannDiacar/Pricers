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

/*! \file discreteintegrals.hpp
    \brief integrals on non uniform grids
*/

#ifndef quantlib_discrete_integrals_hpp
#define quantlib_discrete_integrals_hpp

#include "array.hpp"
#include "integral.hpp"
#include "null.hpp"

namespace QuantLib {

    /*! References:
        Levy, D. Numerical Integration
        http://www2.math.umd.edu/~dlevy/classes/amsc466/lecture-notes/integration-chap.pdf
    */
    class DiscreteTrapezoidIntegral {
      public:
        double operator()(const Array& x, const Array& f) const;
    };

    class DiscreteSimpsonIntegral {
      public:
        double operator()(const Array& x, const Array& f) const;
    };

    class DiscreteTrapezoidIntegrator: public Integrator {
      public:
        explicit DiscreteTrapezoidIntegrator(size_t evaluations)
        : Integrator(Null<double>(), evaluations) {}

      protected:
        double integrate(const ext::function<double(double)>& f, double a, double b) const override;
    };

    class DiscreteSimpsonIntegrator: public Integrator {
      public:
        explicit DiscreteSimpsonIntegrator(size_t evaluations)
        : Integrator(Null<double>(), evaluations) {}

      protected:
        double integrate(const ext::function<double(double)>& f, double a, double b) const override;
    };
}
#endif
