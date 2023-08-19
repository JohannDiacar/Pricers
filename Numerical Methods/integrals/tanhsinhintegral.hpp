/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2021 Klaus Spanderen

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

/*! \file tanhsinhintegral.hpp
*/

#ifndef quantlib_tanh_sinh_integral_hpp
#define quantlib_tanh_sinh_integral_hpp

#include "types.hpp"
#include "null.hpp"
#include "integral.hpp"

#if BOOST_VERSION >= 106900
#define QL_BOOST_HAS_TANH_SINH

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <limits>

namespace QuantLib {
    using tanh_sinh = boost::math::quadrature::tanh_sinh<double>;

    /*! The tanh-sinh quadrature routine provided by boost is a rapidly convergent
        numerical integration scheme for holomorphic integrands. The tolerance
        is used against the error estimate for the L1 norm of the integral.
    */
    class TanhSinhIntegral : public Integrator {
      public:
        TanhSinhIntegral(
            double relTolerance = std::sqrt(std::numeric_limits<double>::epsilon()),
            size_t maxRefinements = 15,
            double minComplement = std::numeric_limits<double>::min() * 4
            )
      : Integrator(QL_MAX_double, Null<size_t>()),
        relTolerance_(relTolerance),
        tanh_sinh_(maxRefinements, minComplement) {}

      protected:
        double integrate(const ext::function<double(double)>& f, double a, double b)
        const override {
            double error;
            double value = tanh_sinh_.integrate(f, a, b, relTolerance_, &error);
            setAbsoluteError(error);

            return value;
        }

      private:
        const double relTolerance_;
        mutable tanh_sinh tanh_sinh_;
    };
}

#else

namespace QuantLib {

    class TanhSinhIntegral : public Integrator {
      public:
        TanhSinhIntegral(double relTolerance = 0, size_t maxRefinements = 0, double minComplement = 0)
        : Integrator(Null<double>(), Null<size_t>()) {
            QL_FAIL("boost version 1.69 or higher is required in order to use TanhSinhIntegral");
        }

      protected:
        double integrate(const ext::function<double(double)>& f, double a, double b)
        const override { return 0.0; }
    };

}

#endif

#endif
